/*
MIT License
Copyright (c) 2021 Germán Molina
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/
use crate::Float;
use calendar::Date;
use communication_protocols::error_handling::ErrorHandling;
use communication_protocols::simulation_model::SimulationModel;
use weather::Weather;

use crate::construction::discretize_construction;
use crate::heating_cooling::calc_cooling_heating_power;
use crate::surface::ThermalSurface;
use crate::zone::ThermalZone;
use simple_model::{Boundary, SimpleModel, SimulationState, SimulationStateHeader};

pub struct ThermalModel {
    /// All the Thermal Zones in the model
    pub zones: Vec<ThermalZone>,

    /// All the surfaces in the model
    pub surfaces: Vec<ThermalSurface>,

    /// All the Fenestrations in the model
    pub fenestrations: Vec<ThermalSurface>,

    /// The number of steps that this model needs
    /// to take in order to advance one step of the main
    /// simulation.
    pub dt_subdivisions: usize,

    /// The model's dt (i.e., main_dt / self.dt_subdivisions)
    pub dt: Float,
}

impl ErrorHandling for ThermalModel {
    fn module_name() -> &'static str {
        "Thermal model"
    }
}

impl SimulationModel for ThermalModel {
    type Type = Self;

    /// Creates a new ThermalModel from a SimpleModel.
    ///    
    /// # Inputs:
    /// * model: the SimpleModel that the model represents
    /// * state: the SimulationState attached to the SimpleModel
    /// * n: the number of timesteps per hour taken by the main simulation.
    fn new(
        model: &SimpleModel,
        state: &mut SimulationStateHeader,
        n: usize,
    ) -> Result<Self, String> {
        /* CREATE ALL ZONES, ONE PER SPACE */
        let mut thermal_zones: Vec<ThermalZone> = Vec::with_capacity(model.spaces.len());
        for (i, space) in model.spaces.iter().enumerate() {
            // Add the zone to the model... this pushes it to the sate
            // as well
            thermal_zones.push(ThermalZone::from_space(space, state, i));
        }

        /* FIND MODEL TIMESTEP */
        // choose the smallest timestep in all constructions
        let max_dx = 0.04; // 4cm
        let min_dt = 60.; // 60 seconds

        let mut n_subdivisions: usize = 1;
        let main_dt = 60. * 60. / n as Float;

        // Store the dts and n_nodes somwehere. Take note of the largest
        // number of subditivions required
        let mut all_n_elements: Vec<Vec<usize>> = Vec::with_capacity(model.constructions.len());
        for construction in &model.constructions {
            let (mut found_n_subdivisions, n_elements) =
                discretize_construction(construction, main_dt, max_dx, min_dt);
            found_n_subdivisions *= n_subdivisions;
            if found_n_subdivisions > n_subdivisions {
                n_subdivisions = found_n_subdivisions;
            }
            all_n_elements.push(n_elements);
        }

        // This is the model's dt now. When marching
        let mut dt = 60. * 60. / (n as Float * n_subdivisions as Float);

        // safety..?
        dt *= 0.5;
        n_subdivisions *= 2;

        /* CREATE SURFACES USING THE MINIMUM TIMESTEP */
        // The rationale here is the following: We find the minimum
        // timestep (or maximum timestep_subdivisions), and that will be the
        // dt for the whole model. Some constructions will have a larger
        // dt (due to their discretization scheme), but simulating them
        // with a smaller (i.e. the official) dt is of no harm.

        // For the Thermal Model
        let mut thermal_surfaces: Vec<ThermalSurface> = Vec::with_capacity(model.surfaces.len());

        for (i, surface) in model.surfaces.iter().enumerate() {
            let construction_index = *surface.construction.index().unwrap();

            let thermal_surface = match ThermalSurface::new_surface(
                state,
                surface,
                dt,
                &all_n_elements[construction_index],
            ) {
                Ok(v) => v,
                Err(e) => return Err(e),
            };

            thermal_surfaces.push(thermal_surface);

            // Match surface and zones
            if let Ok(b) = surface.front_boundary() {
                thermal_surfaces[i].set_front_boundary(b.clone());
            }
            if let Ok(b) = surface.back_boundary() {
                thermal_surfaces[i].set_back_boundary(b.clone());
            }
        }

        let mut thermal_fenestrations: Vec<ThermalSurface> =
            Vec::with_capacity(model.fenestrations.len());
        for (i, fenestration) in model.fenestrations.iter().enumerate() {
            let construction_index = *fenestration.construction.index().unwrap();

            let thermal_surface = match ThermalSurface::new_fenestration(
                state,
                fenestration,
                dt,
                &all_n_elements[construction_index],
            ) {
                Ok(v) => v,
                Err(e) => return Err(e),
            };

            thermal_fenestrations.push(thermal_surface);

            // Match surface and zones
            if let Ok(b) = fenestration.front_boundary() {
                thermal_fenestrations[i].set_front_boundary(b.clone());
            }
            if let Ok(b) = fenestration.back_boundary() {
                thermal_fenestrations[i].set_back_boundary(b.clone());
            }
        }

        Ok(ThermalModel {
            zones: thermal_zones,
            surfaces: thermal_surfaces,
            fenestrations: thermal_fenestrations,
            dt_subdivisions: n_subdivisions,
            dt,
        })
    }

    /// Advances one main_timestep through time. That is,
    /// it performs `self.dt_subdivisions` steps, advancing
    /// `self.dt` seconds in each of them.
    fn march(
        &self,
        mut date: Date,
        weather: &dyn Weather,
        model: &SimpleModel,
        state: &mut SimulationState,
    ) -> Result<(), String> {
        // Iterate through all the sub-subdivitions
        for _sub_index in 0..self.dt_subdivisions {
            // advance in time
            date.add_seconds(self.dt);
            let current_weather = weather.get_weather_data(date);
            
            let t_out = match current_weather.dry_bulb_temperature {
                Some(v) => v,
                None => return Err(
                    "Trying to march on Thermal Model, but dry bulb temperature was not provided"
                    .to_string(),
                ),
            };            

            let t_current = self.get_current_zones_temperatures(state);
            // let (a_before, b_before, c_before) = self.calculate_zones_abc(model, state);
            // let t_current = self.estimate_zones_future_temperatures(&t_current, &a_before, &b_before, &c_before, self.dt);

            /* UPDATE SURFACE'S TEMPERATURES */
            for i in 0..self.surfaces.len() {
                // get surface
                let s = &self.surfaces[i];

                // find t_in and t_out of surface.
                let t_front = match s.front_boundary() {
                    Some(b) => match b {
                        Boundary::Space(space) => t_current[*space.index().unwrap()],
                        Boundary::Ground => unimplemented!(),
                    },
                    None => t_out,
                };
                let t_back = match s.back_boundary() {
                    Some(b) => match b {
                        Boundary::Space(space) => t_current[*space.index().unwrap()], //self.zones[z_index].temperature(model, state),
                        Boundary::Ground => unimplemented!(),
                    },
                    None => t_out,
                };

                // Update temperatures
                let (q_front, q_back) = s.march(state, t_front, t_back, self.dt);
                model.surfaces[i].set_front_convective_heat_flow(state, q_front);
                model.surfaces[i].set_back_convective_heat_flow(state, q_back);
            } // end of iterating surface

            // What  if they are open???
            for i in 0..self.fenestrations.len() {
                // get surface
                let s = &self.fenestrations[i];

                // find t_in and t_out of surface.
                let t_front = match s.front_boundary() {
                    Some(b) => match b {
                        Boundary::Space(space) => t_current[*space.index().unwrap()],
                        Boundary::Ground => unimplemented!(),
                    },
                    None => t_out,
                };
                let t_back = match s.back_boundary() {
                    Some(b) => match b {
                        Boundary::Space(space) => t_current[*space.index().unwrap()],
                        Boundary::Ground => unimplemented!(),
                    },
                    None => t_out,
                };

                // Update temperatures
                let (q_front, q_back) = s.march(state, t_front, t_back, self.dt);
                model.fenestrations[i].set_front_convective_heat_flow(state, q_front);
                model.fenestrations[i].set_back_convective_heat_flow(state, q_back);

            } // end of iterating surface

            /* UPDATE ZONES' TEMPERATURE */
            // This is done analytically.
            let (a, b, c) = self.calculate_zones_abc(model, state);            
            let future_temperatures =
                self.estimate_zones_future_temperatures(&t_current, &a, &b, &c, self.dt);
            for (i, zone) in self.zones.iter().enumerate() {
                zone.reference_space
                    .set_dry_bulb_temperature(state, future_temperatures[i]);
            }
        } // End of 'in each sub-timestep-subdivision'

        Ok(())
    }
}

impl ThermalModel {
    /// Retrieves the dt_subdivisions (i.e. the
    /// number of substimesteps per timestep of this
    /// model)
    pub fn dt_subdivisions(&self) -> usize {
        self.dt_subdivisions
    }

    /// Retrieves a ThermalZone
    pub fn get_thermal_zone(&self, index: usize) -> Result<&ThermalZone, String> {
        if index >= self.zones.len() {
            return ThermalModel::internal_error(format!(
                "Ouf of bounds: Thermal Zone number {} does not exist",
                index
            ));
        }

        Ok(&self.zones[index])
    }

    /// Retrieves a ThermalSurface
    pub fn get_thermal_surface(&self, index: usize) -> Result<&ThermalSurface, String> {
        if index >= self.surfaces.len() {
            return ThermalModel::internal_error(format!(
                "Ouf of bounds: Thermal Surface number {} does not exist",
                index
            ));
        }

        Ok(&self.surfaces[index])
    }

    /// Retrieves a THermalFenestration
    pub fn get_thermal_fenestration(&self, index: usize) -> Result<&ThermalSurface, String> {
        if index >= self.fenestrations.len() {
            return ThermalModel::internal_error(format!(
                "Ouf of bounds: Thermal Surface number {} does not exist",
                index
            ));
        }

        Ok(&self.fenestrations[index])
    }

    /// This estimation assumes nothing changes during this time.
    /// This is self evidently wrong, as we know that, for example, the surface temperatures
    /// will change together with the zone air temperature. However, in short periods of time
    /// this can actually work.
    ///
    /// Everything starts from the following equation, representing a heat balance over
    /// the air and the contents of the Thermal zone.
    ///
    /// ```math
    /// C_{zone}\frac{dT(t)}{dt} = \displaystyle\sum_{i=loads}{Q_i} + \displaystyle\sum_{i=surf.}{h_iA_i(T_i-T)}+\displaystyle\sum_{i=otherzones}{\dot{m_i}C_p(T_i-T)}+\dot{m}_{inf}C_p(T_{out}-T)+\dot{m}_{supplied}C_p(T_{sup}-T)
    /// ```
    /// Which can be expanded into the following
    ///
    /// ```math
    /// C_{zone}\frac{dT(t)}{dt} = A - B T
    /// ```
    ///
    /// Where $`A`$ and $`B`$ are constant terms (at least they can be assumed to be during this brief period of time).
    ///
    /// ```math
    /// A = \displaystyle\sum_{i=loads}{Q_i} + \displaystyle\sum_{i=surf.}{h_iA_i T_i}+\displaystyle\sum_{i=otherzones}{\dot{m_i}C_pT_i}+\dot{m}_{inf}C_pT_{out}+\dot{m}_{supplied}C_pT_{sup}
    /// ```
    ///
    /// ```math
    /// B= \displaystyle\sum_{i=surf.}{h_iA_i}+\displaystyle\sum_{i=otherzones}{\dot{m_i}C_p}+\dot{m}_{inf}C_p +\dot{m}_{supplied}C_p
    /// ```
    ///
    /// And so, (solving the differential equation) the Temperature $`T`$ at a time $`t`$ into the future
    /// can be estimated based on the current Temperature of the zone ($`T_{current}`$) and the following
    /// equation:
    ///
    /// ```math
    ///  T(t) = \frac{A}{B} + \left( T_{current} - \frac{A}{B} \right)e^{-\frac{B}{C_{zone}}t}
    /// ```
    ///
    /// And the average temperature during the same periood is:
    /// ```math
    /// \frac{\displaystyle\int_{0}^t{T(t)dt}}{t} = \frac{A}{B}+\frac{C_{zone}\left(T_{current}-\frac{A}{B}\right)}{Bt}\left(1-e^{-\frac{Bt}{C_{zone}}} \right)
    /// ```
    fn calculate_zones_abc(
        &self,
        model: &SimpleModel,
        state: &SimulationState,
    ) -> (Vec<Float>, Vec<Float>, Vec<Float>) {
        let nzones = self.zones.len();
        // Initialize vectors containing a and b
        let mut a = vec![0.0; nzones];
        let mut b = vec![0.0; nzones];
        let mut c = vec![0.0; nzones];

        /* Qi */
        // Heating/Cooling
        for hvac in model.hvacs.iter() {
            for (target_space_index, heating_cooling) in calc_cooling_heating_power(hvac, state) {
                a[target_space_index] += heating_cooling;
            }
            // heating through air supply?
        }

        // Heating/Cooling
        for luminaire in model.luminaires.iter() {
            if let Ok(target_space) = luminaire.target_space() {
                let target_space_index = *target_space.index().unwrap();
                let consumption = luminaire
                    .power_consumption(state)
                    .expect("Luminaire has no Power Consumption state");
                a[target_space_index] += consumption;
            }
        }

        // Other
        for (i, zone) in self.zones.iter().enumerate() {
            /* INFILTRATION AND VENTILATION */
            // ... should we always asume that the Cp for
            // outside and inside are equal?  I will
            // assume them different... profiling can
            // tell us if this makes the program slow.
            let space = &model.spaces[i];
            // let t_zone = space.dry_bulb_temperature(state).expect("Zone has no Temperature!");
            let cp_zone = gas_properties::air::specific_heat();

            // infiltration from outside
            if let Some(t_inf_inwards) = space.infiltration_temperature(state) {
                let m_inf = space
                    .infiltration_volume(state)
                    .expect("Space has infiltration temperature but not volume");
                let cp_inf_inwards = gas_properties::air::specific_heat();
                a[i] += m_inf * cp_inf_inwards * t_inf_inwards;
                b[i] += m_inf * cp_zone;
            }

            // ventilation
            if let Some(t_vent_inwards) = space.ventilation_temperature(state) {
                let m_vent = space
                    .ventilation_volume(state)
                    .expect("Space has ventilation temperature but not volume");
                let cp_vent_inwards = gas_properties::air::specific_heat();
                a[i] += m_vent * cp_vent_inwards * t_vent_inwards;
                b[i] += m_vent * cp_zone;
            }

            // Mixing with other zones

            /* CAPACITANCE */
            c[i] = zone.mcp();
        }

        /* SURFACES */
        fn iterate_surfaces(
            surfaces: &[ThermalSurface],
            state: &SimulationState,
            a: &mut Vec<Float>,
            b: &mut Vec<Float>,
        ) {
            for surface in surfaces.iter() {
                let (rs_front, rs_back) = match &surface {
                    ThermalSurface::Fenestration(fen, data) => {
                        let front = fen.front_convection_coefficient(state).unwrap();
                        let back = fen.back_convection_coefficient(state).unwrap();
                        (front + data.r_front, back + data.r_back)
                    }
                    ThermalSurface::Surface(fen, data) => {
                        let front = fen.front_convection_coefficient(state).unwrap();
                        let back = fen.back_convection_coefficient(state).unwrap();
                        (front + data.r_front, back + data.r_back)
                    }
                };

                let ai = surface.area();
                // if front leads to a Zone
                if let Some(Boundary::Space(space)) = surface.front_boundary() {
                    let z_index = space.index().unwrap();
                    let hi = 1. / rs_front;
                    let temp = surface.front_temperature(state);
                    a[*z_index] += hi * ai * temp;
                    b[*z_index] += hi * ai;
                }

                // if back leads to a Zone
                if let Some(Boundary::Space(space)) = surface.back_boundary() {
                    let z_index = space.index().unwrap();
                    let hi = 1. / rs_back;
                    let temp = surface.back_temperature(state);
                    a[*z_index] += hi * ai * temp;
                    b[*z_index] += hi * ai;
                }
            }
        }

        iterate_surfaces(&self.surfaces, state, &mut a, &mut b);
        iterate_surfaces(&self.fenestrations, state, &mut a, &mut b);

        /* AIR MIXTURE WITH OTHER ZONES */
        // unimplemented();

        // RETURN
        (a, b, c)
    }

    /// Retrieves a vector of the current temperatures of all the Zones as
    /// registered in the Simulation State
    fn get_current_zones_temperatures(&self, state: &SimulationState) -> Vec<Float> {
        let nzones = self.zones.len();
        // Initialize return
        let mut ret: Vec<Float> = Vec::with_capacity(nzones);
        for zone in self.zones.iter() {
            let t_current = zone.reference_space.dry_bulb_temperature(state).unwrap();
            ret.push(t_current);
        }
        ret
    }

    /// Uses an analytical solution to estimate an average temperature for each Zone
    /// for the near future. Uses the coefficients $`A`$, $`B`$ and $`C`$
    /// calculated by `calculate_zones_abc` and the Zones' current temperatures
    /// `t_current` as calculated by `get_current_temperatures`.
    #[allow(dead_code)]
    fn estimate_zones_mean_future_temperatures(
        &self,
        t_current: &[Float],
        a: &[Float],
        b: &[Float],
        c: &[Float],
        future_time: Float,
    ) -> Vec<Float> {
        let nzones = self.zones.len();
        // Initialize return
        let mut ret: Vec<Float> = Vec::with_capacity(nzones);

        for i in 0..self.zones.len() {
            let current_temp = t_current[i];
            if b[i].abs() > 1e-9 {
                // is this an apropriate threshold?
                ret.push(
                    a[i] / b[i]
                        + (c[i] * (current_temp - a[i] / b[i]) / future_time / b[i])
                            * (1.0 - (-b[i] * future_time / c[i]).exp()),
                );
            } else {
                ret.push(current_temp);
            }
        }

        ret
    }

    /// Uses an analytical solution to estimate the future Zones temperature
    /// for the near future. Uses the coefficients $`A`$, $`B`$ and $`C`$
    /// calculated by `calculate_zones_abc` and the Zones' current temperatures
    /// `t_current` as calculated by `get_current_temperatures`.
    fn estimate_zones_future_temperatures(
        &self,
        t_current: &[Float],
        a: &[Float],
        b: &[Float],
        c: &[Float],
        future_time: Float,
    ) -> Vec<Float> {
        let nzones = self.zones.len();
        // Initialize return
        let mut ret: Vec<Float> = Vec::with_capacity(nzones);
        for i in 0..nzones {
            if b[i].abs() > 1e-9 {
                // is this an apropriate threshold?
                ret.push(
                    a[i] / b[i] + (t_current[i] - a[i] / b[i]) * (-b[i] * future_time / c[i]).exp(),
                );
            } else {
                // A space that is disconnected from everything... maintains its temperature
                ret.push(t_current[i]);
            }
        }

        ret
    }
}

/***********/
/* TESTING */
/***********/

#[cfg(test)]
mod testing {
    use super::*;
    // use crate::construction::*;

    use calendar::Date;
    use schedule::ScheduleConstant;
    use weather::SyntheticWeather;

    use gas_properties::air;
    use simple_model::{SimulationStateElement, HVAC};
    use simple_test_models::*;

    /// A single-zone test model with walls assumed to have
    /// no mass. It has a closed solution, which is nice.
    ///
    /// There is no sun.
    #[derive(Default)]
    struct SingleZoneTestModel {
        /// volume of the zone (m3)
        zone_volume: Float,

        /// Facade area (m2)
        surface_area: Float,

        /// the R-value of the facade
        facade_r: Float,

        /// Infiltration rate (m3/s)
        infiltration_rate: Float,

        /// Heating power (Watts)
        heating_power: Float,

        /// Lighting power (Watts)
        lighting_power: Float,

        /// Temperature outside of the zone
        temp_out: Float,

        /// Temperature at the beginning
        temp_start: Float,
    }

    impl SingleZoneTestModel {
        fn get_closed_solution(&self) -> Box<impl Fn(Float) -> Float> {
            // heat balance in the form
            // of C*dT/dt = A - B*T
            let rho = air::density(); //kg/m3
            let cp = air::specific_heat(); //J/kg.K
            let u = 1. / self.facade_r;

            let c = self.zone_volume * rho * cp;

            let mut a = self.heating_power
                + self.lighting_power
                + self.temp_out * u * self.surface_area
                + self.infiltration_rate * cp * self.temp_out;

            a /= c;

            let mut b = u * self.surface_area + self.infiltration_rate * cp;
            b /= c;

            let k1 = self.temp_start - a / b;

            let f = move |t: Float| -> Float { a / b + k1 * (-b * t).exp() };

            Box::new(f)
        }
    }

    #[test]
    fn test_calculate_zones_abc() {
        let (simple_model, mut state_header) = get_single_zone_test_building(
            // &mut state,
            &SingleZoneTestBuildingOptions {
                zone_volume: 40.,
                surface_area: 4.,
                material_is_massive: Some(false),
                ..Default::default()
            },
        );

        let n: usize = 1;
        let thermal_model = ThermalModel::new(&simple_model, &mut state_header, n).unwrap();
        let state = state_header.take_values().unwrap();
        // MAP THE STATE
        // model.map_simulation_state(&mut state).unwrap();

        // Test
        let (a, b, c) = thermal_model.calculate_zones_abc(&simple_model, &state);
        assert_eq!(a.len(), 1);
        assert_eq!(c.len(), 1);
        assert_eq!(b.len(), 1);
        assert_eq!(c[0], thermal_model.get_thermal_zone(0).unwrap().mcp());
        let rs_front = simple_model.surfaces[0]
            .front_convection_coefficient(&state)
            .unwrap();
        let hi = 1. / rs_front;
        let temp = thermal_model
            .get_thermal_surface(0)
            .unwrap()
            .front_temperature(&state);
        let area = thermal_model.get_thermal_surface(0).unwrap().area();
        assert_eq!(a[0], area * hi * temp);
        assert_eq!(b[0], area * hi);
    }

    #[test]
    fn test_very_simple_march() {
        let zone_volume = 40.;
        let surface_area = 4.;
        let (simple_model, mut state_header) = get_single_zone_test_building(
            // &mut state,
            &SingleZoneTestBuildingOptions {
                zone_volume,
                surface_area,
                material_is_massive: Some(false),
                ..Default::default()
            },
        );

        let n: usize = 60;
        let main_dt = 60. * 60. / n as Float;
        let thermal_model = ThermalModel::new(&simple_model, &mut state_header, n).unwrap();

        let mut state = state_header.take_values().unwrap();

        //println!("DT_SUBDIVISIONS = {}", model.dt_subdivisions);
        // MAP THE STATE
        // model.map_simulation_state(&mut state_).unwrap();

        /* START THE TEST */
        let construction = &simple_model.constructions[0];
        assert!(thermal_model.surfaces[0].is_massive());

        let rs_front = simple_model.surfaces[0]
            .front_convection_coefficient(&state)
            .unwrap();
        let rs_back = simple_model.surfaces[0]
            .back_convection_coefficient(&state)
            .unwrap();
        let r = construction.r_value().unwrap() + rs_front + rs_back;

        // Initial T of the zone
        let t_start = thermal_model.zones[0]
            .reference_space
            .dry_bulb_temperature(&state)
            .unwrap();

        let t_out: Float = 30.0; // T of surroundings

        // test model
        let tester = SingleZoneTestModel {
            zone_volume,
            surface_area,
            facade_r: r,
            temp_out: t_out,
            temp_start: t_start,
            ..SingleZoneTestModel::default()
        };
        let exp_fn = tester.get_closed_solution();

        let mut weather = SyntheticWeather::new();
        weather.dry_bulb_temperature = Box::new(ScheduleConstant::new(t_out));

        let mut date = Date {
            day: 1,
            hour: 0.0,
            month: 1,
        };

        // March:

        for i in 0..800 {
            let time = (i as Float) * main_dt;
            date.add_seconds(time);

            let found = thermal_model.zones[0]
                .reference_space
                .dry_bulb_temperature(&state)
                .unwrap();

            thermal_model
                .march(date, &weather, &simple_model, &mut state)
                .unwrap();

            // Get exact solution.
            let exp = exp_fn(time);

            //assert!((exp - found).abs() < 0.05);
            let max_error = 0.15;
            let diff = (exp - found).abs();
            println!("{},{}, {}", exp, found, diff);
            assert!(diff < max_error);
        }
    }
    /// END OF TEST_MODEL_MARCH

    #[test]
    fn test_march_with_window() {
        let surface_area = 4.;
        let window_area = 1.;
        let zone_volume = 40.;

        let (simple_model, mut state_header) = get_single_zone_test_building(
            // &mut state,
            &SingleZoneTestBuildingOptions {
                zone_volume,
                surface_area,
                window_area,
                material_is_massive: Some(false),
                ..Default::default()
            },
        );

        // Finished model the SimpleModel
        let n: usize = 6;
        let main_dt = 60. * 60. / n as Float;
        let thermal_model = ThermalModel::new(&simple_model, &mut state_header, n).unwrap();

        let mut state = state_header.take_values().unwrap();

        // MAP THE STATE
        // model.map_simulation_state(&mut state).unwrap();

        // START TESTING.
        let construction = &simple_model.constructions[0];
        // assert!(!model.surfaces[0].is_massive());

        let rs_front = simple_model.surfaces[0]
            .front_convection_coefficient(&state)
            .unwrap();
        let rs_back = simple_model.surfaces[0]
            .back_convection_coefficient(&state)
            .unwrap();
        let r = construction.r_value().unwrap() + rs_front + rs_back;

        // Initial T of the zone
        let t_start = thermal_model.zones[0]
            .reference_space
            .dry_bulb_temperature(&state)
            .unwrap();

        let t_out: Float = 30.0; // T of surroundings

        let mut weather = SyntheticWeather::new();
        weather.dry_bulb_temperature = Box::new(ScheduleConstant::new(t_out));

        let dt = main_dt;

        let mut date = Date {
            day: 1,
            hour: 0.0,
            month: 1,
        };

        // test model
        let tester = SingleZoneTestModel {
            zone_volume,
            surface_area, // the window is a hole on the wall... does not add area
            facade_r: r,
            temp_out: t_out,
            temp_start: t_start,
            ..SingleZoneTestModel::default()
        };
        let exp_fn = tester.get_closed_solution();

        // March:
        for i in 0..80 {
            let time = (i as Float) * dt;
            date.add_seconds(time);

            let found = thermal_model.zones[0]
                .reference_space
                .dry_bulb_temperature(&state)
                .unwrap();

            thermal_model
                .march(date, &weather, &simple_model, &mut state)
                .unwrap();

            // Get exact solution.
            let exp = exp_fn(time);
            let max_error = 0.15;
            println!("{}, {}", exp, found);
            assert!((exp - found).abs() < max_error);
        }
    }

    #[test]
    fn test_model_march_with_window_and_luminaire() {
        let surface_area = 4.;
        let zone_volume = 40.;
        let lighting_power = 100.;

        let (simple_model, mut state_header) = get_single_zone_test_building(
            // &mut state,
            &SingleZoneTestBuildingOptions {
                zone_volume,
                surface_area,
                lighting_power,
                material_is_massive: Some(false),
                ..Default::default()
            },
        );

        // Finished model the SimpleModel

        let n: usize = 20;
        let main_dt = 60. * 60. / n as Float;
        let thermal_model = ThermalModel::new(&simple_model, &mut state_header, n).unwrap();

        let mut state = state_header.take_values().unwrap();

        // MAP THE STATE
        // model.map_simulation_state(&mut state).unwrap();

        // turn the lights on
        let lum_state_i = simple_model.luminaires[0]
            .power_consumption_index()
            .unwrap();
        // state.update_value(lum_state_i, SimulationStateElement::LuminairePowerConsumption(0, lighting_power));
        state[lum_state_i] = lighting_power;

        // START TESTING.
        let construction = &simple_model.constructions[0];
        // assert!(!model.surfaces[0].is_massive());
        println!("IS MASSIVE??? {}", thermal_model.surfaces[0].is_massive());

        let rs_front = simple_model.surfaces[0]
            .front_convection_coefficient(&state)
            .unwrap();
        let rs_back = simple_model.surfaces[0]
            .back_convection_coefficient(&state)
            .unwrap();
        let r = construction.r_value().unwrap() + rs_front + rs_back;

        // Initial T of the zone
        let t_start = thermal_model.zones[0]
            .reference_space
            .dry_bulb_temperature(&state)
            .unwrap();

        let t_out: Float = 30.0; // T of surroundings

        // test model
        let tester = SingleZoneTestModel {
            zone_volume,
            surface_area, // the window is a hole on the wall... does not add area
            lighting_power,
            facade_r: r,
            temp_out: t_out,
            temp_start: t_start,
            ..SingleZoneTestModel::default()
        };
        let exp_fn = tester.get_closed_solution();

        let mut weather = SyntheticWeather::new();
        weather.dry_bulb_temperature = Box::new(ScheduleConstant::new(t_out));

        let dt = main_dt; // / model.dt_subdivisions() as Float;

        let mut date = Date {
            day: 1,
            hour: 0.0,
            month: 1,
        };

        // March:
        for i in 0..800 {
            let time = (i as Float) * dt;
            date.add_seconds(time);

            let found = thermal_model.zones[0]
                .reference_space
                .dry_bulb_temperature(&state)
                .unwrap();

            thermal_model
                .march(date, &weather, &simple_model, &mut state)
                .unwrap();

            // Get exact solution.
            let exp = exp_fn(time);

            let max_error = 0.4;
            println!("{}, {}", exp, found);
            assert!((exp - found).abs() < max_error);
        }
    }

    #[test]
    fn test_model_march_with_window_and_heater() {
        let surface_area = 4.;
        let zone_volume = 40.;
        let heating_power = 100.;

        let (simple_model, mut state_header) = get_single_zone_test_building(
            // &mut state,
            &SingleZoneTestBuildingOptions {
                zone_volume,
                surface_area,
                heating_power,
                material_is_massive: Some(false),
                ..Default::default()
            },
        );

        // Finished model the SimpleModel

        let n: usize = 20;
        let main_dt = 60. * 60. / n as Float;
        let thermal_model = ThermalModel::new(&simple_model, &mut state_header, n).unwrap();

        let mut state = state_header.take_values().unwrap();
        // MAP THE STATE
        // model.map_simulation_state(&mut state).unwrap();

        // turn the heater on
        if let HVAC::ElectricHeater(heater) = &simple_model.hvacs[0] {
            let hvac_state_i = heater.heating_cooling_consumption_index().unwrap();
            state[hvac_state_i] = heating_power;
        }

        // START TESTING.
        let construction = &simple_model.constructions[0];
        // assert!(!model.surfaces[0].is_massive());

        let rs_front = simple_model.surfaces[0]
            .front_convection_coefficient(&state)
            .unwrap();
        let rs_back = simple_model.surfaces[0]
            .back_convection_coefficient(&state)
            .unwrap();
        let r = construction.r_value().unwrap() + rs_front + rs_back;

        // Initial T of the zone
        let t_start = thermal_model.zones[0]
            .reference_space
            .dry_bulb_temperature(&state)
            .unwrap();
        let t_out: Float = 30.0; // T of surroundings

        // test model
        let tester = SingleZoneTestModel {
            zone_volume,
            surface_area, // the window is a hole on the wall... does not add area
            heating_power,
            facade_r: r,
            temp_out: t_out,
            temp_start: t_start,
            ..SingleZoneTestModel::default()
        };
        let exp_fn = tester.get_closed_solution();

        let mut weather = SyntheticWeather::new();
        weather.dry_bulb_temperature = Box::new(ScheduleConstant::new(t_out));

        let dt = main_dt; // / model.dt_subdivisions() as Float;

        let mut date = Date {
            day: 1,
            hour: 0.0,
            month: 1,
        };

        // March:
        for i in 0..800 {
            let time = (i as Float) * dt;
            date.add_seconds(time);

            let found = thermal_model.zones[0]
                .reference_space
                .dry_bulb_temperature(&state)
                .unwrap();

            thermal_model
                .march(date, &weather, &simple_model, &mut state)
                .unwrap();

            // Get exact solution.
            let exp = exp_fn(time);

            let max_error = 0.2;
            println!("{}, {}", exp, found);
            assert!((exp - found).abs() < max_error);
        }
    }

    #[test]
    fn test_model_march_with_window_heater_and_infiltration() {
        let surface_area = 4.;
        let zone_volume = 40.;
        let heating_power = 100.;
        let infiltration_rate = 0.1;
        let t_out: Float = 30.0; // T of surroundings

        let (simple_model, mut state_header) = get_single_zone_test_building(
            // &mut state,
            &SingleZoneTestBuildingOptions {
                zone_volume,
                surface_area,
                heating_power,
                infiltration_rate,
                material_is_massive: Some(false),
                ..Default::default()
            },
        );

        // Finished model the SimpleModel

        let n: usize = 20;
        let main_dt = 60. * 60. / n as Float;
        let thermal_model = ThermalModel::new(&simple_model, &mut state_header, n).unwrap();

        // Set infiltration
        let inf_vol_index = state_header.push(
            SimulationStateElement::SpaceInfiltrationVolume(0),
            infiltration_rate,
        );
        simple_model.spaces[0].set_infiltration_volume_index(inf_vol_index);
        let inf_temp_index = state_header.push(
            SimulationStateElement::SpaceInfiltrationTemperature(0),
            t_out,
        );
        simple_model.spaces[0].set_infiltration_temperature_index(inf_temp_index);

        // MAP THE STATE

        let mut state = state_header.take_values().unwrap();

        // turn the heater on
        if let HVAC::ElectricHeater(heater) = &simple_model.hvacs[0] {
            let hvac_state_i = heater.heating_cooling_consumption_index().unwrap();
            state[hvac_state_i] = heating_power;
        }

        // START TESTING.
        let construction = &simple_model.constructions[0];
        // assert!(!model.surfaces[0].is_massive());

        let rs_front = simple_model.surfaces[0]
            .front_convection_coefficient(&state)
            .unwrap();
        let rs_back = simple_model.surfaces[0]
            .back_convection_coefficient(&state)
            .unwrap();
        let r = construction.r_value().unwrap() + rs_front + rs_back;

        // Initial T of the zone
        let t_start = thermal_model.zones[0]
            .reference_space
            .dry_bulb_temperature(&state)
            .unwrap();

        // test model
        let tester = SingleZoneTestModel {
            zone_volume,
            surface_area, // the window is a hole on the wall... does not add area
            heating_power,
            facade_r: r,
            temp_out: t_out,
            temp_start: t_start,
            infiltration_rate,
            ..SingleZoneTestModel::default()
        };
        let exp_fn = tester.get_closed_solution();

        let mut weather = SyntheticWeather::new();
        weather.dry_bulb_temperature = Box::new(ScheduleConstant::new(t_out));

        let dt = main_dt; // / model.dt_subdivisions() as Float;

        let mut date = Date {
            day: 1,
            hour: 0.0,
            month: 1,
        };

        // March:
        for i in 0..800 {
            let time = (i as Float) * dt;
            date.add_seconds(time);

            let found = thermal_model.zones[0]
                .reference_space
                .dry_bulb_temperature(&state)
                .unwrap();

            thermal_model
                .march(date, &weather, &simple_model, &mut state)
                .unwrap();

            // Get exact solution.
            let exp = exp_fn(time);

            let max_error = 0.1;
            println!("{}, {}", exp, found);
            assert!((exp - found).abs() < max_error);
        }
    }

    #[test]
    fn test_model_march_solar_radiation(){

        let surface_area = 20. * 3.;
        let zone_volume = 600.;
        

        let (simple_model, mut state_header) = get_single_zone_test_building(
            // &mut state,
            &SingleZoneTestBuildingOptions {
                zone_volume,
                surface_area,                
                material_is_massive: Some(true),
                ..Default::default()
            },
        );

        // Finished model the SimpleModel

        let n: usize = 20;
        // let main_dt = 60. * 60. / n as Float;
        let thermal_model = ThermalModel::new(&simple_model, &mut state_header, n).unwrap();

        let mut state = state_header.take_values().unwrap();

        

        let outdoor_temp = [15.085000, 15.170000, 15.255000, 15.340000, 15.425000, 15.510000, 15.595000, 15.680000, 15.765000, 15.850000, 15.935000, 16.020000, 16.105000, 16.190000, 16.275000, 16.360000, 16.445000, 16.530000, 16.615000, 16.700000, 16.620000, 16.540000, 16.460000, 16.380000, 16.300000, 16.220000, 16.140000, 16.060000, 15.980000, 15.900000, 15.820000, 15.740000, 15.660000, 15.580000, 15.500000, 15.420000, 15.340000, 15.260000, 15.180000, 15.100000, 15.035000, 14.970000, 14.905000, 14.840000, 14.775000, 14.710000, 14.645000, 14.580000, 14.515000, 14.450000, 14.385000, 14.320000, 14.255000, 14.190000, 14.125000, 14.060000, 13.995000, 13.930000, 13.865000, 13.800000, 13.745000, 13.690000, 13.635000, 13.580000, 13.525000, 13.470000, 13.415000, 13.360000, 13.305000, 13.250000, 13.195000, 13.140000, 13.085000, 13.030000, 12.975000, 12.920000, 12.865000, 12.810000, 12.755000, 12.700000, 12.650000, 12.600000, 12.550000, 12.500000, 12.450000, 12.400000, 12.350000, 12.300000, 12.250000, 12.200000, 12.150000, 12.100000, 12.050000, 12.000000, 11.950000, 11.900000, 11.850000, 11.800000, 11.750000, 11.700000, 11.645000, 11.590000, 11.535000, 11.480000, 11.425000, 11.370000, 11.315000, 11.260000, 11.205000, 11.150000, 11.095000, 11.040000, 10.985000, 10.930000, 10.875000, 10.820000, 10.765000, 10.710000, 10.655000, 10.600000, 10.720000, 10.840000, 10.960000, 11.080000, 11.200000, 11.320000, 11.440000, 11.560000, 11.680000, 11.800000, 11.920000, 12.040000, 12.160000, 12.280000, 12.400000, 12.520000, 12.640000, 12.760000, 12.880000, 13.000000, 13.250000, 13.500000, 13.750000, 14.000000, 14.250000, 14.500000, 14.750000, 15.000000, 15.250000, 15.500000, 15.750000, 16.000000, 16.250000, 16.500000, 16.750000, 17.000000, 17.250000, 17.500000, 17.750000, 18.000000, 18.150000, 18.300000, 18.450000, 18.600000, 18.750000, 18.900000, 19.050000, 19.200000, 19.350000, 19.500000, 19.650000, 19.800000, 19.950000, 20.100000, 20.250000, 20.400000, 20.550000, 20.700000, 20.850000, 21.000000, 21.150000, 21.300000, 21.450000, 21.600000, 21.750000, 21.900000, 22.050000, 22.200000, 22.350000, 22.500000, 22.650000, 22.800000, 22.950000, 23.100000, 23.250000, 23.400000, 23.550000, 23.700000, 23.850000, 24.000000, 24.110000, 24.220000, 24.330000, 24.440000, 24.550000, 24.660000, 24.770000, 24.880000, 24.990000, 25.100000, 25.210000, 25.320000, 25.430000, 25.540000, 25.650000, 25.760000, 25.870000, 25.980000, 26.090000, 26.200000, 26.240000, 26.280000, 26.320000, 26.360000, 26.400000, 26.440000, 26.480000, 26.520000, 26.560000, 26.600000, 26.640000, 26.680000, 26.720000, 26.760000, 26.800000, 26.840000, 26.880000, 26.920000, 26.960000, 27.000000, 27.100000, 27.200000, 27.300000, 27.400000, 27.500000, 27.600000, 27.700000, 27.800000, 27.900000, 28.000000, 28.100000, 28.200000, 28.300000, 28.400000, 28.500000, 28.600000, 28.700000, 28.800000, 28.900000, 29.000000, 29.050000, 29.100000, 29.150000, 29.200000, 29.250000, 29.300000, 29.350000, 29.400000, 29.450000, 29.500000, 29.550000, 29.600000, 29.650000, 29.700000, 29.750000, 29.800000, 29.850000, 29.900000, 29.950000, 30.000000, 29.975000, 29.950000, 29.925000, 29.900000, 29.875000, 29.850000, 29.825000, 29.800000, 29.775000, 29.750000, 29.725000, 29.700000, 29.675000, 29.650000, 29.625000, 29.600000, 29.575000, 29.550000, 29.525000, 29.500000, 29.475000, 29.450000, 29.425000, 29.400000, 29.375000, 29.350000, 29.325000, 29.300000, 29.275000, 29.250000, 29.225000, 29.200000, 29.175000, 29.150000, 29.125000, 29.100000, 29.075000, 29.050000, 29.025000, 29.000000, 28.945000, 28.890000, 28.835000, 28.780000, 28.725000, 28.670000, 28.615000, 28.560000, 28.505000, 28.450000, 28.395000, 28.340000, 28.285000, 28.230000, 28.175000, 28.120000, 28.065000, 28.010000, 27.955000, 27.900000, 27.805000, 27.710000, 27.615000, 27.520000, 27.425000, 27.330000, 27.235000, 27.140000, 27.045000, 26.950000, 26.855000, 26.760000, 26.665000, 26.570000, 26.475000, 26.380000, 26.285000, 26.190000, 26.095000, 26.000000, 25.900000, 25.800000, 25.700000, 25.600000, 25.500000, 25.400000, 25.300000, 25.200000, 25.100000, 25.000000, 24.900000, 24.800000, 24.700000, 24.600000, 24.500000, 24.400000, 24.300000, 24.200000, 24.100000, 24.000000, 23.925000, 23.850000, 23.775000, 23.700000, 23.625000, 23.550000, 23.475000, 23.400000, 23.325000, 23.250000, 23.175000, 23.100000, 23.025000, 22.950000, 22.875000, 22.800000, 22.725000, 22.650000, 22.575000, 22.500000, 22.325000, 22.150000, 21.975000, 21.800000, 21.625000, 21.450000, 21.275000, 21.100000, 20.925000, 20.750000, 20.575000, 20.400000, 20.225000, 20.050000, 19.875000, 19.700000, 19.525000, 19.350000, 19.175000, 19.000000, 18.900000, 18.800000, 18.700000, 18.600000, 18.500000, 18.400000, 18.300000, 18.200000, 18.100000, 18.000000, 17.900000, 17.800000, 17.700000, 17.600000, 17.500000, 17.400000, 17.300000, 17.200000, 17.100000, 17.000000, 16.950000, 16.900000, 16.850000, 16.800000, 16.750000, 16.700000, 16.650000, 16.600000, 16.550000, 16.500000, 16.450000, 16.400000, 16.350000, 16.300000, 16.250000, 16.200000, 16.150000, 16.100000, 16.050000, 16.000000, 15.950000, 15.900000, 15.850000, 15.800000, 15.750000, 15.700000, 15.650000, 15.600000, 15.550000, 15.500000, 15.450000, 15.400000, 15.350000, 15.300000, 15.250000, 15.200000, 15.150000, 15.100000, 15.050000, 15.000000, 14.930000, 14.860000, 14.790000, 14.720000, 14.650000, 14.580000, 14.510000, 14.440000, 14.370000, 14.300000, 14.230000, 14.160000, 14.090000, 14.020000, 13.950000, 13.880000, 13.810000, 13.740000, 13.670000, 13.600000, 13.600000, 13.600000, 13.600000, 13.600000, 13.600000, 13.600000, 13.600000, 13.600000, 13.600000, 13.600000, 13.600000, 13.600000, 13.600000, 13.600000, 13.600000, 13.600000, 13.600000, 13.600000, 13.600000, 13.600000, 13.520000, 13.440000, 13.360000, 13.280000, 13.200000, 13.120000, 13.040000, 12.960000, 12.880000, 12.800000, 12.720000, 12.640000, 12.560000, 12.480000, 12.400000, 12.320000, 12.240000, 12.160000, 12.080000, 12.000000, 12.000000, 12.000000, 12.000000, 12.000000, 12.000000, 12.000000, 12.000000, 12.000000, 12.000000, 12.000000, 12.000000, 12.000000, 12.000000, 12.000000, 12.000000, 12.000000, 12.000000, 12.000000, 12.000000, 12.000000, 12.005000, 12.010000, 12.015000, 12.020000, 12.025000, 12.030000, 12.035000, 12.040000, 12.045000, 12.050000, 12.055000, 12.060000, 12.065000, 12.070000, 12.075000, 12.080000, 12.085000, 12.090000, 12.095000, 12.100000, 11.995000, 11.890000, 11.785000, 11.680000, 11.575000, 11.470000, 11.365000, 11.260000, 11.155000, 11.050000, 10.945000, 10.840000, 10.735000, 10.630000, 10.525000, 10.420000, 10.315000, 10.210000, 10.105000, 10.000000, 10.050000, 10.100000, 10.150000, 10.200000, 10.250000, 10.300000, 10.350000, 10.400000, 10.450000, 10.500000, 10.550000, 10.600000, 10.650000, 10.700000, 10.750000, 10.800000, 10.850000, 10.900000, 10.950000, 11.000000, 11.255000, 11.510000, 11.765000, 12.020000, 12.275000, 12.530000, 12.785000, 13.040000, 13.295000, 13.550000, 13.805000, 14.060000, 14.315000, 14.570000, 14.825000, 15.080000, 15.335000, 15.590000, 15.845000, 16.100000, 16.295000, 16.490000, 16.685000, 16.880000, 17.075000, 17.270000, 17.465000, 17.660000, 17.855000, 18.050000, 18.245000, 18.440000, 18.635000, 18.830000, 19.025000, 19.220000, 19.415000, 19.610000, 19.805000, 20.000000, 20.100000, 20.200000, 20.300000, 20.400000, 20.500000, 20.600000, 20.700000, 20.800000, 20.900000, 21.000000, 21.100000, 21.200000, 21.300000, 21.400000, 21.500000, 21.600000, 21.700000, 21.800000, 21.900000, 22.000000, 22.130000, 22.260000, 22.390000, 22.520000, 22.650000, 22.780000, 22.910000, 23.040000, 23.170000, 23.300000, 23.430000, 23.560000, 23.690000, 23.820000, 23.950000, 24.080000, 24.210000, 24.340000, 24.470000, 24.600000, 24.770000, 24.940000, 25.110000, 25.280000, 25.450000, 25.620000, 25.790000, 25.960000, 26.130000, 26.300000, 26.470000, 26.640000, 26.810000, 26.980000, 27.150000, 27.320000, 27.490000, 27.660000, 27.830000, 28.000000, 28.020000, 28.040000, 28.060000, 28.080000, 28.100000, 28.120000, 28.140000, 28.160000, 28.180000, 28.200000, 28.220000, 28.240000, 28.260000, 28.280000, 28.300000, 28.320000, 28.340000, 28.360000, 28.380000, 28.400000, 28.420000, 28.440000, 28.460000, 28.480000, 28.500000, 28.520000, 28.540000, 28.560000, 28.580000, 28.600000, 28.620000, 28.640000, 28.660000, 28.680000, 28.700000, 28.720000, 28.740000, 28.760000, 28.780000, 28.800000, 28.810000, 28.820000, 28.830000, 28.840000, 28.850000, 28.860000, 28.870000, 28.880000, 28.890000, 28.900000, 28.910000, 28.920000, 28.930000, 28.940000, 28.950000, 28.960000, 28.970000, 28.980000, 28.990000, 29.000000, 28.950000, 28.900000, 28.850000, 28.800000, 28.750000, 28.700000, 28.650000, 28.600000, 28.550000, 28.500000, 28.450000, 28.400000, 28.350000, 28.300000, 28.250000, 28.200000, 28.150000, 28.100000, 28.050000, 28.000000, 27.945000, 27.890000, 27.835000, 27.780000, 27.725000, 27.670000, 27.615000, 27.560000, 27.505000, 27.450000, 27.395000, 27.340000, 27.285000, 27.230000, 27.175000, 27.120000, 27.065000, 27.010000, 26.955000, 26.900000, 26.805000, 26.710000, 26.615000, 26.520000, 26.425000, 26.330000, 26.235000, 26.140000, 26.045000, 25.950000, 25.855000, 25.760000, 25.665000, 25.570000, 25.475000, 25.380000, 25.285000, 25.190000, 25.095000, 25.000000, 24.900000, 24.800000, 24.700000, 24.600000, 24.500000, 24.400000, 24.300000, 24.200000, 24.100000, 24.000000, 23.900000, 23.800000, 23.700000, 23.600000, 23.500000, 23.400000, 23.300000, 23.200000, 23.100000, 23.000000, 22.975000, 22.950000, 22.925000, 22.900000, 22.875000, 22.850000, 22.825000, 22.800000, 22.775000, 22.750000, 22.725000, 22.700000, 22.675000, 22.650000, 22.625000, 22.600000, 22.575000, 22.550000, 22.525000, 22.500000, 22.475000, 22.450000, 22.425000, 22.400000, 22.375000, 22.350000, 22.325000, 22.300000, 22.275000, 22.250000, 22.225000, 22.200000, 22.175000, 22.150000, 22.125000, 22.100000, 22.075000, 22.050000, 22.025000, 22.000000, 21.850000, 21.700000, 21.550000, 21.400000, 21.250000, 21.100000, 20.950000, 20.800000, 20.650000, 20.500000, 20.350000, 20.200000, 20.050000, 19.900000, 19.750000, 19.600000, 19.450000, 19.300000, 19.150000, 19.000000, 18.950000, 18.900000, 18.850000, 18.800000, 18.750000, 18.700000, 18.650000, 18.600000, 18.550000, 18.500000, 18.450000, 18.400000, 18.350000, 18.300000, 18.250000, 18.200000, 18.150000, 18.100000, 18.050000, 18.000000, 17.750000, 17.500000, 17.250000, 17.000000, 16.750000, 16.500000, 16.250000, 16.000000, 15.750000, 15.500000, 15.250000, 15.000000, 14.750000, 14.500000, 14.250000, 14.000000, 13.750000, 13.500000, 13.250000, 13.000000, 13.200000, 13.400000, 13.600000, 13.800000, 14.000000, 14.200000, 14.400000, 14.600000, 14.800000, 15.000000, 15.200000, 15.400000, 15.600000, 15.800000, 16.000000, 16.200000, 16.400000, 16.600000, 16.800000, 17.000000, 16.875000, 16.750000, 16.625000, 16.500000, 16.375000, 16.250000, 16.125000, 16.000000, 15.875000, 15.750000, 15.625000, 15.500000, 15.375000, 15.250000, 15.125000, 15.000000, 14.875000, 14.750000, 14.625000, 14.500000, 14.475000, 14.450000, 14.425000, 14.400000, 14.375000, 14.350000, 14.325000, 14.300000, 14.275000, 14.250000, 14.225000, 14.200000, 14.175000, 14.150000, 14.125000, 14.100000, 14.075000, 14.050000, 14.025000, 14.000000, 14.000000, 14.000000, 14.000000, 14.000000, 14.000000, 14.000000, 14.000000, 14.000000, 14.000000, 14.000000, 14.000000, 14.000000, 14.000000, 14.000000, 14.000000, 14.000000, 14.000000, 14.000000, 14.000000, 14.000000, 13.960000, 13.920000, 13.880000, 13.840000, 13.800000, 13.760000, 13.720000, 13.680000, 13.640000, 13.600000, 13.560000, 13.520000, 13.480000, 13.440000, 13.400000, 13.360000, 13.320000, 13.280000, 13.240000, 13.200000, 13.140000, 13.080000, 13.020000, 12.960000, 12.900000, 12.840000, 12.780000, 12.720000, 12.660000, 12.600000, 12.540000, 12.480000, 12.420000, 12.360000, 12.300000, 12.240000, 12.180000, 12.120000, 12.060000, 12.000000, 12.050000, 12.100000, 12.150000, 12.200000, 12.250000, 12.300000, 12.350000, 12.400000, 12.450000, 12.500000, 12.550000, 12.600000, 12.650000, 12.700000, 12.750000, 12.800000, 12.850000, 12.900000, 12.950000, 13.000000, 13.200000, 13.400000, 13.600000, 13.800000, 14.000000, 14.200000, 14.400000, 14.600000, 14.800000, 15.000000, 15.200000, 15.400000, 15.600000, 15.800000, 16.000000, 16.200000, 16.400000, 16.600000, 16.800000, 17.000000, 17.150000, 17.300000, 17.450000, 17.600000, 17.750000, 17.900000, 18.050000, 18.200000, 18.350000, 18.500000, 18.650000, 18.800000, 18.950000, 19.100000, 19.250000, 19.400000, 19.550000, 19.700000, 19.850000, 20.000000, 20.150000, 20.300000, 20.450000, 20.600000, 20.750000, 20.900000, 21.050000, 21.200000, 21.350000, 21.500000, 21.650000, 21.800000, 21.950000, 22.100000, 22.250000, 22.400000, 22.550000, 22.700000, 22.850000, 23.000000, 23.130000, 23.260000, 23.390000, 23.520000, 23.650000, 23.780000, 23.910000, 24.040000, 24.170000, 24.300000, 24.430000, 24.560000, 24.690000, 24.820000, 24.950000, 25.080000, 25.210000, 25.340000, 25.470000, 25.600000, 25.670000, 25.740000, 25.810000, 25.880000, 25.950000, 26.020000, 26.090000, 26.160000, 26.230000, 26.300000, 26.370000, 26.440000, 26.510000, 26.580000, 26.650000, 26.720000, 26.790000, 26.860000, 26.930000, 27.000000, 27.050000, 27.100000, 27.150000, 27.200000, 27.250000, 27.300000, 27.350000, 27.400000, 27.450000, 27.500000, 27.550000, 27.600000, 27.650000, 27.700000, 27.750000, 27.800000, 27.850000, 27.900000, 27.950000, 28.000000, 27.970000, 27.940000, 27.910000, 27.880000, 27.850000, 27.820000, 27.790000, 27.760000, 27.730000, 27.700000, 27.670000, 27.640000, 27.610000, 27.580000, 27.550000, 27.520000, 27.490000, 27.460000, 27.430000, 27.400000, 27.380000, 27.360000, 27.340000, 27.320000, 27.300000, 27.280000, 27.260000, 27.240000, 27.220000, 27.200000, 27.180000, 27.160000, 27.140000, 27.120000, 27.100000, 27.080000, 27.060000, 27.040000, 27.020000, 27.000000, 27.000000, 27.000000, 27.000000, 27.000000, 27.000000, 27.000000, 27.000000, 27.000000, 27.000000, 27.000000, 27.000000, 27.000000, 27.000000, 27.000000, 27.000000, 27.000000, 27.000000, 27.000000, 27.000000, 27.000000, 26.945000, 26.890000, 26.835000, 26.780000, 26.725000, 26.670000, 26.615000, 26.560000, 26.505000, 26.450000, 26.395000, 26.340000, 26.285000, 26.230000, 26.175000, 26.120000, 26.065000, 26.010000, 25.955000, 25.900000, 25.805000, 25.710000, 25.615000, 25.520000, 25.425000, 25.330000, 25.235000, 25.140000, 25.045000, 24.950000, 24.855000, 24.760000, 24.665000, 24.570000, 24.475000, 24.380000, 24.285000, 24.190000, 24.095000, 24.000000, 23.950000, 23.900000, 23.850000, 23.800000, 23.750000, 23.700000, 23.650000, 23.600000, 23.550000, 23.500000, 23.450000, 23.400000, 23.350000, 23.300000, 23.250000, 23.200000, 23.150000, 23.100000, 23.050000, 23.000000, 22.950000, 22.900000, 22.850000, 22.800000, 22.750000, 22.700000, 22.650000, 22.600000, 22.550000, 22.500000, 22.450000, 22.400000, 22.350000, 22.300000, 22.250000, 22.200000, 22.150000, 22.100000, 22.050000, 22.000000, 21.900000, 21.800000, 21.700000, 21.600000, 21.500000, 21.400000, 21.300000, 21.200000, 21.100000, 21.000000, 20.900000, 20.800000, 20.700000, 20.600000, 20.500000, 20.400000, 20.300000, 20.200000, 20.100000, 20.000000, 19.900000, 19.800000, 19.700000, 19.600000, 19.500000, 19.400000, 19.300000, 19.200000, 19.100000, 19.000000, 18.900000, 18.800000, 18.700000, 18.600000, 18.500000, 18.400000, 18.300000, 18.200000, 18.100000, 18.000000, 17.920000, 17.840000, 17.760000, 17.680000, 17.600000, 17.520000, 17.440000, 17.360000, 17.280000, 17.200000, 17.120000, 17.040000, 16.960000, 16.880000, 16.800000, 16.720000, 16.640000, 16.560000, 16.480000, 16.400000, 16.330000, 16.260000, 16.190000, 16.120000, 16.050000, 15.980000, 15.910000, 15.840000, 15.770000, 15.700000, 15.630000, 15.560000, 15.490000, 15.420000, 15.350000, 15.280000, 15.210000, 15.140000, 15.070000, 15.000000];
        let thermal_heat_gain = [-2207.935358, -2179.483538, -2152.297408, -2126.166830, -2100.944370, -2076.520058, -2052.808116, -2029.739434, -2007.256955, -1985.312660, -1963.865514, -1942.879995, -1922.325029, -1902.173184, -1882.400051, -1862.983765, -1843.904616, -1825.144748, -1806.687904, -1788.519219, -1807.829005, -1825.193741, -1841.292016, -1856.369533, -1870.595847, -1884.095351, -1896.963431, -1909.275613, -1921.093166, -1932.466754, -1943.438921, -1954.045852, -1964.318659, -1974.284337, -1983.966501, -1993.385954, -2002.561135, -2011.508484, -2020.242728, -2028.777123, -2032.573439, -2036.417875, -2040.247668, -2044.046346, -2047.804310, -2051.515850, -2055.177600, -2058.787688, -2062.345229, -2065.850006, -2069.302259, -2072.702549, -2076.051653, -2079.350503, -2082.600127, -2085.801620, -2088.956118, -2092.064771, -2095.128739, -2098.149172, -2098.989253, -2099.922834, -2100.906842, -2101.927174, -2102.974242, -2104.041142, -2105.122692, -2106.214889, -2107.314572, -2108.419204, -2109.526724, -2110.635443, -2111.743963, -2112.851124, -2113.955956, -2115.057646, -2116.155512, -2117.248982, -2118.337574, -2119.420884, -2119.876689, -2120.386257, -2120.929666, -2121.500033, -2122.092511, -2122.703483, -2123.330128, -2123.970180, -2124.621779, -2125.283368, -2125.953630, -2126.631435, -2127.315805, -2128.005887, -2128.700935, -2129.400289, -2130.103362, -2130.809634, -2131.518637, -2132.229953, -2132.677044, -2133.080713, -2133.455515, -2133.806176, -2134.135951, -2134.447221, -2134.741809, -2135.021158, -2135.286439, -2135.538624, -2135.778529, -2136.006852, -2136.224194, -2206.857516, -2235.041280, -2259.377024, -2282.468244, -2304.830931, -2326.772610, -2348.478051, -2335.549854, -2324.857385, -2313.204410, -2298.649224, -2284.878696, -2272.037759, -2260.075714, -2248.930884, -2238.540284, -2228.843470, -2238.930905, -2279.101455, -2300.205617, -2321.164609, -2342.814638, -2365.049901, -2387.751612, -2410.803761, -2434.098402, -2457.895503, -2483.138296, -2487.879628, -2491.649785, -2495.514770, -2499.463032, -2503.447626, -2507.408420, -2511.281845, -2515.004938, -2518.517062, -2515.013858, -2510.002550, -2503.918155, -2496.868859, -2488.917743, -2505.922668, -2501.152420, -2493.661467, -2484.656450, -2474.310555, -2480.187133, -2483.584851, -2485.003486, -2484.620310, -2482.559594, -2478.913375, -2473.752737, -2467.134378, -2459.104728, -2449.702682, -2434.770577, -2434.635426, -2419.021065, -2400.763047, -2380.732557, -2359.136011, -2336.118605, -2311.790352, -2286.239337, -2259.538887, -2236.999345, -2213.336040, -2188.628095, -2162.930276, -2136.290471, -2108.751645, -2080.353045, -2051.130992, -2021.119449, -2003.915263, -1988.426791, -1974.093020, -1960.671510, -1947.994125, -1935.940057, -1924.418549, -1913.359163, -1902.705844, -1892.413096, -1882.443384, -1895.159420, -1894.391667, -1892.554243, -1890.272893, -1887.681466, -1884.861673, -1881.867240, -1878.735590, -1875.493914, -1872.162628, -1866.509753, -1860.391490, -1853.947961, -1847.237694, -1840.302993, -1833.176209, -1825.883083, -1821.033283, -1814.658665, -1808.120947, -1817.520950, -1825.829467, -1833.358734, -1840.232262, -1846.537436, -1852.340139, -1857.692358, -1862.636512, -1867.208095, -1871.437421, -1876.812014, -1882.123802, -1887.304617, -1892.337567, -1897.214754, -1901.933529, -1906.494624, -1910.901176, -1915.158174, -1919.272156, -1907.700471, -1896.900242, -1886.602769, -1876.714051, -1867.173185, -1857.939425, -1848.985392, -1840.293149, -1831.851697, -1823.655210, -1814.998816, -1806.462842, -1798.088229, -1789.892394, -1781.887917, -1774.083793, -1766.485929, -1759.097360, -1751.918416, -1744.946922, -1748.079322, -1750.715715, -1753.057234, -1755.177343, -1757.123612, -1758.928144, -1760.613099, -1762.193904, -1763.681265, -1765.082501, -1766.665033, -1768.232712, -1769.776823, -1771.295388, -1772.787005, -1774.250434, -1775.684431, -1777.087685, -1778.458808, -1779.796322, -1801.676984, -1822.352518, -1842.170010, -1861.265955, -1879.735211, -1897.647839, -1915.057863, -1933.236769, -1950.389156, -1967.171168, -1980.515807, -1992.939569, -2004.628489, -2015.660395, -2026.091999, -2035.967143, -2045.321147, -2054.183317, -2062.578513, -2070.528186, -2076.312933, -2081.719088, -2086.751626, -2091.420219, -2095.734411, -2099.703363, -2103.335757, -2106.639793, -2109.623204, -2112.293287, -2116.229942, -2120.216883, -2124.179762, -2128.098003, -2131.959685, -2135.757813, -2139.488401, -2143.149405, -2146.740083, -2150.260589, -2162.413876, -2174.007960, -2182.696573, -2210.325320, -2240.969169, -2273.532896, -2307.639247, -2343.025091, -2379.494073, -2416.892109, -2451.407690, -2485.752718, -2520.032585, -2554.226273, -2588.301225, -2622.220074, -2655.943696, -2689.432608, -2722.647617, -2755.550095, -2795.561713, -2834.572600, -2872.749524, -2910.123556, -2946.705912, -2982.495891, -2996.646629, -3026.411780, -3056.797433, -3086.840934, -3108.518000, -3127.714896, -3144.857981, -3160.031792, -3173.285478, -3184.648156, -3194.136805, -3201.760536, -3207.523088, -3211.424382, -3184.940467, -3182.412977, -3180.344869, -3177.162515, -3172.688388, -3166.829251, -3159.538719, -3150.799168, -3140.613290, -3098.037682, -3079.194729, -3060.679216, -3041.755753, -3021.877452, -3000.873166, -2978.632022, -2955.065778, -2930.097453, -2903.654429, -2875.663662, -2838.233031, -2799.520832, -2759.279451, -2717.348998, -2673.565699, -2598.257279, -2544.757136, -2491.008289, -2434.832521, -2373.571725, -2343.138216, -2317.872872, -2297.297333, -2279.711828, -2264.672643, -2251.830474, -2240.647354, -2194.648522, -2177.991974, -2164.437437, -2174.907879, -2185.217075, -2195.626251, -2206.152703, -2216.798426, -2227.558633, -2238.425579, -2249.390326, -2260.443618, -2271.576328, -2282.779682, -2294.045373, -2305.365602, -2316.733091, -2328.141070, -2339.583255, -2351.053819, -2362.547368, -2374.058902, -2385.583793, -2381.089201, -2377.689137, -2375.019476, -2372.955248, -2371.408047, -2370.311431, -2369.613207, -2369.271052, -2369.249817, -2369.519762, -2370.055348, -2370.834373, -2371.837339, -2373.046973, -2374.447864, -2376.026168, -2377.769382, -2379.666158, -2381.706149, -2383.879884, -2374.438260, -2365.839835, -2357.832799, -2350.329115, -2343.265539, -2336.593812, -2330.275463, -2324.278848, -2318.577329, -2313.148071, -2307.971218, -2303.029299, -2298.306796, -2293.789808, -2289.465802, -2285.323402, -2281.352235, -2277.542796, -2273.886341, -2270.374793, -2267.879473, -2265.499847, -2263.234435, -2261.078875, -2259.028684, -2257.079433, -2255.226821, -2253.466703, -2251.795111, -2250.208250, -2248.702499, -2247.274403, -2245.920669, -2244.638158, -2243.423874, -2242.274962, -2241.188696, -2240.162473, -2239.193810, -2238.280332, -2240.840642, -2243.190503, -2245.414650, -2247.540347, -2249.586189, -2251.565618, -2253.488774, -2255.363544, -2257.196198, -2258.991814, -2260.754556, -2262.487882, -2264.194688, -2265.877417, -2267.538145, -2269.178645, -2270.800439, -2272.404839, -2273.992979, -2275.565843, -2263.069617, -2251.497857, -2240.534525, -2230.075274, -2220.047138, -2210.395891, -2201.079390, -2192.063798, -2183.321276, -2174.828474, -2166.565499, -2158.515186, -2150.662565, -2142.994457, -2135.499167, -2128.166248, -2120.986307, -2113.950853, -2107.052177, -2100.283243, -2109.838949, -2118.453319, -2126.468995, -2134.000294, -2141.126169, -2147.904550, -2154.379879, -2160.587364, -2166.555588, -2172.308204, -2177.865095, -2183.243187, -2188.457049, -2193.519337, -2198.441136, -2203.232222, -2207.901274, -2212.456038, -2216.903465, -2221.249820, -2209.458434, -2198.633143, -2188.422959, -2178.714238, -2169.428324, -2160.507369, -2151.906874, -2143.591480, -2135.532383, -2127.705661, -2120.091129, -2112.671535, -2105.431967, -2098.359414, -2091.442433, -2084.670881, -2078.035719, -2071.528836, -2065.142923, -2058.871361, -2051.208549, -2043.720964, -2036.379574, -2029.171547, -2022.086877, -2015.117347, -2008.255978, -2001.496703, -1994.834159, -1988.263548, -1981.780538, -1975.381183, -1969.061871, -1962.819277, -1956.650326, -1950.552165, -1944.522136, -1938.557757, -1932.656704, -1926.816796, -1944.347915, -1960.499476, -1975.731556, -1990.205263, -2004.033100, -2017.298590, -2030.066550, -2042.388925, -2054.308366, -2065.860563, -2077.075844, -2087.980308, -2098.596655, -2172.893509, -2209.956872, -2243.304796, -2275.220838, -2306.212340, -2336.593844, -2366.559820, -2363.394348, -2362.101785, -2361.175923, -2355.796170, -2350.875565, -2346.656908, -2343.112828, -2340.198972, -2337.866128, -2336.065323, -2353.955622, -2401.462651, -2430.112833, -2458.433445, -2487.309806, -2516.639049, -2546.306223, -2576.199406, -2606.214826, -2636.624306, -2654.694482, -2653.635184, -2652.060349, -2651.003600, -2650.382840, -2650.100152, -2650.056622, -2650.158189, -2650.317311, -2650.453138, -2643.827481, -2635.837911, -2626.917742, -2617.161336, -2606.619833, -2620.775990, -2613.641242, -2603.836266, -2592.602226, -2580.103743, -2577.641205, -2573.278188, -2567.335403, -2559.926430, -2551.130011, -2541.003599, -2529.590820, -2516.925870, -2503.036294, -2487.944847, -2483.013826, -2462.010929, -2438.112252, -2412.273507, -2384.737755, -2355.675386, -2325.214796, -2293.458108, -2260.489742, -2226.381528, -2210.013488, -2191.356467, -2170.889962, -2148.808060, -2125.256695, -2100.352253, -2074.191498, -2046.857310, -2018.422284, -2002.943524, -1992.612780, -1983.892512, -1976.326369, -1986.377127, -1983.449748, -1980.111167, -1976.902367, -1973.858136, -1970.983127, -1968.269504, -1960.742619, -1953.734462, -1947.107505, -1940.802831, -1934.775952, -1928.991401, -1923.419868, -1918.036521, -1912.819957, -1907.751477, -1899.357105, -1890.478912, -1881.311545, -1871.923807, -1862.362255, -1852.660409, -1842.843520, -1832.931261, -1822.939368, -1812.880698, -1791.936317, -1772.581545, -1753.221859, -1734.169890, -1715.388299, -1696.840839, -1678.499830, -1660.343609, -1642.354927, -1624.519872, -1607.703749, -1591.144147, -1574.781016, -1558.585758, -1542.537704, -1526.621493, -1510.825671, -1495.141903, -1479.564478, -1464.090010, -1484.329900, -1502.562303, -1519.401652, -1535.108446, -1549.867791, -1563.820644, -1577.079941, -1589.739676, -1601.880358, -1613.572433, -1624.226276, -1634.438238, -1644.292687, -1653.847865, -1663.152496, -1672.247679, -1681.167925, -1689.941823, -1698.592585, -1707.138558, -1707.764189, -1708.446888, -1709.152978, -1709.870518, -1710.590274, -1711.304188, -1712.004746, -1712.684735, -1713.337159, -1713.955228, -1713.690964, -1713.205524, -1712.533370, -1711.683459, -1710.660142, -1709.465233, -1708.099069, -1706.561099, -1704.850237, -1702.965080, -1711.916296, -1720.390307, -1728.473972, -1736.201292, -1743.595599, -1750.673919, -1757.449238, -1763.931796, -1770.129877, -1776.050324, -1782.192642, -1788.198081, -1794.051598, -1799.751463, -1805.297786, -1810.691560, -1815.934212, -1821.027358, -1825.972683, -1830.771861, -1847.612066, -1856.294788, -1870.242266, -1883.985145, -1897.351267, -1910.347922, -1922.990024, -1935.294483, -1947.278072, -1958.956673, -1968.754391, -1977.980194, -1986.735892, -1995.067318, -2003.009230, -2010.589531, -2017.831512, -2024.755142, -2031.377885, -2037.715242, -2047.108337, -2056.127994, -2079.493798, -2107.732150, -2120.563590, -2149.322754, -2180.431888, -2213.042140, -2246.883104, -2281.761133, -2315.789153, -2350.074042, -2384.584674, -2419.248068, -2453.994648, -2488.759093, -2523.480199, -2558.100421, -2592.565352, -2626.823208, -2669.092159, -2710.436658, -2751.013209, -2790.844817, -2829.934596, -2868.273687, -2885.216814, -2917.637889, -2950.667811, -2983.374604, -3007.261285, -3028.661651, -3048.037946, -3065.491842, -3081.085327, -3094.857769, -3106.834732, -3117.032803, -3125.462460, -3132.129898, -3135.988893, -3138.050296, -3104.980744, -3096.616531, -3088.898929, -3080.228380, -3070.412342, -3059.346740, -3046.976576, -3033.276877, -3017.479787, -3000.127576, -2981.160906, -2960.605313, -2938.448115, -2914.671701, -2857.083317, -2823.533427, -2790.495621, -2756.431324, -2708.235064, -2659.647693, -2610.148717, -2559.479911, -2507.417910, -2453.750903, -2398.264807, -2312.582113, -2247.479741, -2178.742644, -2141.654260, -2109.832259, -2083.167039, -2059.752328, -2039.104770, -2020.855124, -2004.474273, -1953.930403, -1932.639256, -1914.450895, -1901.497506, -1889.906984, -1879.485193, -1870.083478, -1861.586273, -1853.900320, -1846.948325, -1840.664936, -1834.994040, -1829.886865, -1825.300597, -1821.197332, -1817.543281, -1814.308139, -1811.464576, -1808.987837, -1806.855397, -1805.046686, -1803.542858, -1802.326587, -1822.789090, -1841.853412, -1860.008653, -1877.434543, -1894.254117, -1910.556872, -1926.410821, -1941.869317, -1956.975236, -1971.763696, -1986.263902, -2000.500459, -2014.494318, -2028.263487, -2041.823572, -2055.188193, -2068.369316, -2081.377512, -2094.222171, -2106.911675, -2099.228153, -2092.792328, -2087.165624, -2082.195826, -2077.777447, -2073.832755, -2070.301851, -2067.137045, -2064.299417, -2061.756582, -2059.481159, -2057.449698, -2055.641889, -2054.039976, -2052.628310, -2051.392995, -2050.321618, -2049.403025, -2048.627143, -2047.984835, -2085.926914, -2121.241479, -2154.786332, -2186.864321, -2217.686700, -2247.410111, -2276.155976, -2304.021498, -2331.086414, -2357.417398, -2383.071076, -2408.096157, -2432.534998, -2456.424775, -2479.798380, -2502.685120, -2525.111273, -2547.100531, -2568.674363, -2589.852314, -2523.958159, -2463.804495, -2407.248643, -2353.630780, -2302.492689, -2253.496493, -2206.381622, -2160.940496, -2117.003627, -2074.429936, -2033.100148, -1992.912112, -1953.777389, -1915.618695, -1878.367948, -1841.964752, -1806.355186, -1771.490847, -1737.328067, -1703.827265, -1736.965239, -1766.525770, -1793.772818, -1819.168508, -1843.030139, -1865.589346, -1887.022869, -1907.469964, -1927.043063, -1945.834705, -1963.922243, -1981.371182, -1998.237601, -2014.569959, -2030.410478, -2045.796207, -2060.759866, -2075.330518, -2089.534109, -2103.393916, -2095.474631, -2088.611993, -2082.379298, -2076.639077, -2071.298344, -2066.290284, -2061.564692, -2057.082594, -2052.812969, -2048.730634, -2044.814815, -2041.048142, -2037.415921, -2033.905606, -2030.506382, -2027.208863, -2024.004841, -2020.887097, -2017.849247, -2014.885613, -2006.198116, -1997.913833, -1989.918166, -1982.168830, -1974.635461, -1967.294907, -1960.128743, -1953.121852, -1946.261545, -1939.536987, -1932.938790, -1926.458731, -1920.089537, -1913.824723, -1907.658468, -1901.585516, -1895.601093, -1889.700845, -1883.880785, -1878.137246, -1881.915286, -1885.224625, -1888.234212, -1891.003139, -1893.572153, -1895.971120, -1898.222914, -1900.345631, -1902.353935, -1904.259938, -1906.073800, -1907.804152, -1909.458409, -1911.042993, -1912.563519, -1914.024924, -1915.431583, -1916.787387, -1918.095820, -1919.360013, -1924.847185, -1930.032259, -1935.001863, -1939.788071, -1944.413811, -1948.896502, -1953.249964, -1957.485509, -1961.612613, -1965.639356, -1969.572731, -1973.418859, -1977.183151, -1980.870431, -2063.994097, -2095.244781, -2122.245151, -2147.905129, -2172.760590, -2197.122388, -2197.030670, -2198.230880, -2200.300937, -2197.478007, -2194.281442, -2191.604019, -2189.478625, -2187.883655, -2186.786077, -2186.149103, -2205.583424, -2254.884759, -2284.935208, -2314.680808, -2344.943151, -2375.622029, -2406.604752, -2437.781550, -2469.050663, -2500.319286, -2501.748889, -2532.068541, -2541.471715, -2549.693810, -2557.848549, -2565.950355, -2573.969550, -2581.858491, -2589.563269, -2597.028781, -2597.817225, -2597.134155, -2595.397677, -2592.708703, -2589.125261, -2584.680642, -2604.686705, -2603.324460, -2599.235729, -2593.637521, -2596.666481, -2597.770241, -2597.251949, -2595.228730, -2591.782127, -2586.971898, -2580.843637, -2573.433287, -2564.769989, -2554.877984, -2538.008397, -2519.014612, -2514.456011, -2494.844955, -2472.602015, -2448.591387, -2423.013534, -2396.010875, -2367.692380, -2338.146137, -2307.446143, -2275.656331, -2242.833135, -2209.027220, -2174.284701, -2138.648036, -2102.156691, -2064.847659, -2026.755864, -2002.198217, -1981.473759, -1962.358343, -1944.473313, -1944.460113, -1931.388953, -1918.016212, -1904.879805, -1892.017697, -1879.436116, -1867.127865, -1860.452237, -1853.732618, -1847.045105, -1840.409204, -1833.834805, -1827.326521, -1820.885901, -1814.512647, -1808.205327, -1801.961811, -1796.401175, -1790.968360, -1785.617703, -1780.327496, -1775.081102, -1769.865078, -1764.668176, -1759.480770, -1754.294499, -1749.102036, -1749.747479, -1749.669238, -1749.078354, -1748.054617, -1746.652983, -1744.913851, -1742.868411, -1740.541674, -1737.954347, -1735.124056, -1736.331889, -1737.900205, -1739.575809, -1741.256898, -1742.873743, -1744.376512, -1741.098929, -1742.168811, -1743.469385, -1744.803512, -1754.196370, -1763.198976, -1771.901573, -1780.340593, -1788.547384, -1796.552093, -1804.385534, -1812.080114, -1819.670245, -1827.192407, -1832.639349, -1837.743295, -1842.659059, -1847.471508, -1852.250726, -1857.056675, -1911.993421, -1926.490522, -1938.318687, -1949.190733, -1974.391811, -1998.033947, -2020.591378, -2042.285170, -2063.271314, -2083.665181, -2103.555006, -2123.009929, -2142.085142, -2160.825356, -2174.721394, -2187.487327, -2199.410309, -2210.606328, -2221.159129, -2231.132647, -2240.577669, -2249.535700, -2258.041401, -2266.124216, -2278.473545, -2290.473073, -2302.134411, -2313.471995, -2362.502738, -2379.598707, -2393.936909, -2407.069489, -2419.271924, -2430.713900, -2438.773770, -2445.735498, -2451.827099, -2457.154637, -2461.797364, -2465.817384, -2469.265082, -2472.182403, -2474.604973, -2476.563549, -2470.910576, -2465.202182, -2459.342192, -2453.308736, -2447.090487, -2440.681968, -2434.081156, -2427.288180, -2420.304557, -2413.132711, -2404.126337, -2394.631731, -2401.304629, -2393.348522, -2383.909573, -2373.724725, -2362.942122, -2351.657515, -2339.938094, -2327.834140, -2332.727066, -2336.377027, -2342.670023, -2350.557734, -2358.196709, -2365.411362, -2374.656057, -2380.077876, -2384.528490, -2388.149313, -2391.364648, -2394.024479, -2396.137098, -2397.713480, -2398.761933, -2399.288833, -2399.299113, -2398.796610, -2397.784300, -2396.264478, -2401.577025, -2405.801482, -2409.127830, -2411.623611, -2413.336106, -2414.300408, -2414.543673, -2414.087542, -2398.648095, -2393.263972, -2386.840368, -2379.798709, -2372.146981, -2363.876349, -2354.984881, -2345.473483, -2335.344026, -2324.598419, -2313.238116, -2301.263860, -2275.505212, -2249.858441, -2224.084653, -2198.096231, -2171.829843, -2145.236329, -2118.275450, -2090.912985, -2063.119053, -2034.867171, -2011.908697, -1989.142069, -1967.141482, -1945.422094, -1923.874255, -1902.421599, -1881.000536, -1859.555256, -1838.034705, -1816.390692, -1798.070285, -1779.475646, -1760.581958, -1741.356954, -1721.774386, -1702.902490, -1685.081936, -1667.462517, -1649.988090, -1632.641712, -1622.086387, -1612.797366, -1604.400073, -1596.739520, -1589.706209, -1583.217415, -1577.206184, -1559.431693, -1553.023109, -1547.782118, -1568.383696, -1588.636014, -1608.723950, -1628.695691, -1648.580582, -1668.397737, -1688.160217, -1707.877242, -1727.555466, -1747.199768, -1766.813755, -1786.400107, -1805.960817, -1825.497358, -1845.010809, -1864.501942, -1883.971296, -1903.419228, -1922.845952, -1942.251569, -1945.864893, -1949.720900, -1953.735163, -1957.876413, -1962.122676, -1966.457522, -1970.868097, -1975.344018, -1979.876688, -1984.458849, -1989.084281, -1993.747584, -1998.444023, -2003.169405, -2007.919992, -2012.692430, -2017.483688, -2022.291019, -2027.111919, -2031.944096, -2041.191139, -2050.578599, -2060.063421, -2069.629111, -2079.263727, -2088.958086, -2098.704812, -2108.497800, -2118.331882, -2128.202602, -2138.106070, -2148.038850, -2157.997877, -2167.980401, -2177.983933, -2188.006213, -2198.045175, -2208.098925, -2218.165719, -2228.243949, -2227.377190, -2226.805850, -2226.432975, -2226.224824, -2226.157477, -2226.212908, -2226.376920, -2226.637967, -2226.986431, -2227.414150, -2227.914095, -2228.480133, -2229.106864, -2229.789490, -2230.523715, -2231.305669, -2232.131847, -2232.999059, -2233.904385, -2234.845149];
        let incident_solar_radiation = [0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,22.10646184,26.52037076,30.33829044,33.79710130,37.04109062,40.15278299,43.17905067,46.14675988,49.07144879,51.19563413,52.14290011,53.29763049,54.60341139,56.02013467,57.51850948,59.07665290,60.67789384,68.25726950,83.81988936,90.92365164,97.80406506,104.46817846,110.92196916,117.17054514,123.21830720,129.06907951,134.83617675,148.72356897,153.93420568,158.95443680,163.78603711,168.43066209,172.88987384,177.16516255,181.25796418,185.16967522,188.90166504,190.44027874,191.85942187,193.15613074,194.32795866,195.37291154,203.91891758,204.61891625,205.16582491,205.56146725,205.80761819,205.90602062,205.85839916,205.66647165,205.33195858,204.85659096,204.24211687,203.49030686,202.60295838,201.58189951,200.42899200,197.93774987,200.19598775,197.44138075,194.61992859,191.73396061,188.78573340,185.77744128,182.71122544,179.58918199,176.41336904,173.18581289,169.90851357,166.58344960,163.21258231,159.79785958,156.34121917,152.84459163,149.30990292,145.73907670,145.97194049,146.13855212,146.30793020,146.48000509,146.65469839,146.83192360,147.01158667,147.19358657,147.37781580,147.56416092,147.75250304,152.40605388,152.52694955,152.64886990,152.77173962,152.89547830,153.02000109,153.14521931,153.27104110,153.39737209,153.52411601,153.02832911,152.53673955,152.04936794,151.56622827,151.08732863,150.61267198,150.14225691,150.38910672,150.12390348,149.85243731,149.57511197,149.29231584,149.00442363,148.71179812,148.41479194,148.11374943,147.80900841,147.50090217,147.18976138,146.87591617,146.95961376,147.03040108,147.08880358,147.13541550,147.17090804,147.19603827,147.21165879,147.21872809,147.21832134,147.21164144,147.20002965,147.18497536,147.16812382,147.15128070,147.13641183,147.12563622,147.12121008,147.12549989,147.14094233,147.16999019,147.02496009,146.89979494,146.79619967,146.71555349,146.65883821,146.62658429,146.61884169,146.63517927,146.67471283,146.73615793,146.81790071,146.91807853,147.03466207,147.16553244,147.30854825,147.46160046,147.62265410,147.78977781,147.96116267,148.13513218,148.38022739,148.62913154,148.88080937,149.13430685,149.38874676,149.64332326,149.89729593,150.14998368,150.40075870,150.64904066,150.89429116,151.13600865,151.37372358,151.60699407,151.83540186,152.05854859,152.27605248,152.81376150,153.12561597,153.43617948,152.91801629,152.38479491,151.83645122,151.27294047,150.69423544,150.10032469,149.49121114,148.86691070,148.22745109,147.57287076,146.90321796,146.21854987,145.51893179,144.80443651,144.07514361,143.33113892,142.57251400,141.79936565,141.01179550,140.20990960,139.81237053,139.41645201,139.02231086,138.63009901,138.23996296,137.85204337,137.46647450,137.08338372,136.70289098,136.32510826,135.95013896,135.57807730,134.54194697,138.77681406,143.00254432,147.21735673,151.41944153,155.60695532,159.77801560,163.93069461,167.07959490,170.16187426,173.17546450,176.11829608,178.98829450,181.78337611,184.50144338,187.14037933,189.69804124,192.17225337,194.56079848,196.86140822,199.07175191,201.18942358,203.21192704,205.13665849,201.38946838,203.03487284,204.59261842,206.06146465,205.32773213,204.41030375,203.30761192,202.01814515,200.54044395,198.87309566,197.01472779,194.96399956,192.71959123,190.28019049,178.83168243,175.94809029,172.89364914,169.67095124,166.28317387,162.73417381,159.02860269,155.17204922,151.17121533,138.61328325,134.44108107,129.98624010,125.50139638,120.86380121,116.06902341,111.11194598,105.98662062,100.68608844,95.20215877,89.52513638,83.64348630,77.54342589,71.20843781,64.61871138,57.75055459,42.39370004,34.94249822,27.16532184,18.83913944,9.33338590,8.87520948,8.50088377,8.60399710,8.79584231,9.08482304,9.45721446,9.82720772,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,20.12198201,24.64408842,28.52139290,32.00008240,35.23848406,38.32905220,41.32483783,44.25648352,47.14187417,49.72669210,50.61461959,51.71679678,52.97478834,54.34729050,55.80428294,57.32342693,58.88776130,66.50903168,82.08106092,89.21951225,96.13165442,102.82505117,109.30611375,115.58031672,121.65237145,127.52636595,133.31978678,147.20221714,152.43559962,157.47840283,162.33249901,166.99962651,171.48141808,175.77942416,179.89513211,183.82998226,187.58538127,189.14353932,190.58108254,191.89469696,193.08163820,194.13965918,202.69813075,203.41122734,203.96706341,204.36762343,204.61482349,204.71053133,204.65658298,204.45479645,204.10698320,203.61495764,202.98054497,202.20558780,201.29195158,200.24152909,199.05624419,201.04841486,197.90500558,194.68383790,191.38855729,188.02260006,184.58922010,181.09151225,177.53243285,173.91481774,170.24139807,166.51481429,162.73762832,158.91233445,155.04136883,151.12711792,147.17192590,143.17810125,139.14792257,135.08364369,134.99661105,135.87624128,136.74783585,137.61124265,143.23755894,144.03650452,144.82536476,145.60400097,146.37227533,147.13005114,147.87719318,148.61356804,149.33904447,150.05349373,150.75679003,151.44881093,152.12943782,152.79855641,153.45605725,154.10183636,154.73579584,154.38705664,154.03973462,153.69384558,153.34940257,153.00641659,152.66489744,152.32485455,151.98629784,151.64923866,151.31369076,150.97967131,150.93991232,150.74608655,150.54728808,150.34386108,150.13615227,149.92451334,149.70930354,149.49089242,149.26966284,149.28657198,149.29148724,149.28498570,149.26771753,149.24041464,149.20390001,149.15909796,149.10704495,149.04890096,148.98596095,148.91966600,148.85161340,148.78356462,148.71744992,148.65536792,148.59957811,148.55248399,148.51660470,148.49453309,148.48887893,148.32608792,148.18619074,148.07103934,147.98214200,147.92058531,147.88697406,147.88139666,147.90342057,147.95211845,148.02612146,148.12369291,148.24281378,148.38127119,148.53674257,148.70686990,148.88932119,149.08183802,149.28226966,149.48859554,149.69893773,149.68604213,149.66803471,149.64333147,149.61051308,149.56831236,149.51560127,149.45137813,149.37475533,149.28494778,149.18126231,149.06308790,148.92988692,148.78118724,148.61657515,148.43568907,148.23821393,148.02387624,147.79243969,147.54370119,147.27748751,147.12527847,146.96617442,146.80007584,146.62689134,146.44653598,146.25892990,146.06399688,145.86166312,145.65185600,145.43450299,145.20953055,143.07222932,142.87400092,142.67213031,142.46670439,142.25781078,142.04553729,141.82997144,141.61120006,141.38930879,140.73981903,140.08621339,139.42866332,138.76733937,138.10241094,137.43404582,136.76240993,136.08766694,135.40997791,134.72950098,134.04639100,133.36079918,136.60406629,140.50608318,139.73047962,143.66375240,147.58798139,151.50174036,155.40359714,159.29211145,162.69935672,166.04058858,169.31356345,172.51602101,175.64567841,178.70022359,181.67730754,184.57453557,187.38945710,190.11955407,192.76222754,195.31478222,197.77440873,200.13816301,202.40294253,204.56545878,201.07623835,202.94265876,204.71830128,206.40163610,205.74006320,204.90760618,203.90285985,202.72448708,201.37121844,199.84185140,198.13524919,196.25033885,194.18610868,191.94160450,189.51592469,186.90821342,175.06617520,172.03723143,168.84724006,165.49857862,161.99416195,158.33753240,154.53297100,150.58563572,146.29273562,141.84788683,137.22578030,132.43753057,127.48013583,122.35012647,108.18370097,102.65397505,96.93806748,91.02834308,84.91568696,78.58915357,72.03553469,65.23883921,58.17969443,50.83472538,43.17608624,27.32210907,18.94151892,9.33210842,8.87376829,8.43318564,8.52394725,8.70129845,8.97403444,9.33084933,9.69406244,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,24.68201157,28.48486283,31.83544914,34.92455262,37.86050076,40.70343756,43.48691894,46.22983442,48.94282581,49.88341458,50.82685118,51.95180072,53.21011325,54.56714647,55.99735095,57.48148018,65.06922396,80.57549238,87.68000193,94.55849942,101.21881318,107.66756469,113.91039598,119.95215211,125.79702810,131.44868749,137.00932689,150.60140512,155.62129569,160.45400659,165.10129223,169.56479680,173.84607829,177.94662838,181.86788894,185.61126583,187.26823999,188.80465770,190.21750243,191.50428429,192.66297373,193.69194590,202.07596574,202.73798281,203.24718599,203.60533261,203.81414994,203.87534880,203.79063491,203.56171828,203.19032090,202.67818320,202.02706930,201.23877144,200.31511352,199.25795410,196.39605689,193.46771098,195.13615391,191.89977138,188.58862408,185.20620286,181.75581757,178.24062004,174.66362428,171.02772420,167.33570924,163.59027812,159.79405095,155.94957999,152.05935908,148.12583199,144.15139985,140.13842766,136.08925009,136.07441525,136.55420466,137.03225056,137.50843002,142.76600259,143.18042016,143.59170150,143.99975036,144.40446687,144.80574800,145.20348799,145.59757882,145.98791065,146.37437222,146.75685136,147.13523543,147.50941178,147.87926829,148.24469381,148.60557880,148.96181579,149.48634872,149.99324812,150.48253245,150.95423107,151.40838523,151.84504912,152.26429114,152.66619523,153.05086238,153.41841240,153.76898578,154.10274589,154.41988140,154.72060905,155.00517671,155.27386689,155.52700067,155.76494216,155.98810351,156.19695053,157.55759199,158.84857248,160.07145733,161.22796167,162.31997107,163.34956455,163.06039542,164.17251566,165.24070397,166.26767950,167.25645441,168.21036222,169.13308499,170.02867723,170.90158389,171.75664920,172.59911208,173.43458377,174.26900304,175.10856501,175.40893769,175.73305220,176.08666671,176.47521968,176.90366992,177.37635066,191.32199252,191.97285186,192.65048791,193.35562649,194.08845315,194.84867792,195.63561379,196.44826036,197.28538531,198.14559841,199.02741496,199.92930727,200.84974433,201.78722085,201.52770332,201.27873816,201.03850994,200.80538827,200.57791538,200.35479254,200.13486618,199.91711422,199.70063296,199.48462484,199.26838704,199.05130112,198.83282354,198.61247715,208.51821851,208.07789804,207.63417278,207.18688686,206.73590921,206.28113029,205.09391183,203.89909175,202.69675042,201.48698509,200.26990774,199.04564325,197.81432767,196.57610673,195.33113447,194.07957203,192.82158650,191.55734997,190.28703857,189.01083165,187.72891103,186.44146026,185.14866401,183.85070746,182.54777572,181.24005331,179.48795163,177.73281576,180.38783740,178.39644627,176.40538232,174.41496140,172.42549460,170.43728798,168.45064238,166.46585320,164.48321014,162.50299700,161.47607291,160.90987760,160.25862727,159.52119205,159.37456116,158.13946713,156.82016008,155.41581508,154.03542501,152.61631400,151.15813634,149.66055085,148.12322065,146.54581272,144.92799762,143.26944897,141.56984300,139.82885804,138.04617387,136.22147107,134.35443029,132.44473138,130.49205243,128.49606876,126.45645168,124.37286719,118.36616374,115.95663398,113.19958350,110.40586202,107.57462275,104.70496476,101.79592903,98.84649428,95.85557279,92.82200617,89.74456137,86.62192696,83.45271003,80.23543407,76.96853831,73.65037929,70.27923558,66.85331714,63.37078118,59.82975719,56.22838516,52.56487208,50.43324707,48.13477265,45.94487135,43.72999767,41.48781984,39.21570047,36.91066347,34.56936477,32.18807387,29.76267872,27.28873722,24.76161643,22.17679516,19.53046844,16.82071437,14.35092375,12.16644009,9.98793881,7.81685340,5.65526000,5.37462581,5.09464566,4.81541184,4.53702066,4.25952888,3.98281622,3.70619720,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000];
        let exp_zone_air_temp = [25.443431, 25.414320, 25.384897, 25.355170, 25.325147, 25.294836, 25.264247, 25.233386, 25.202262, 25.170883, 25.139255, 25.107386, 25.075283, 25.042955, 25.010408, 24.977652, 24.944693, 24.911541, 24.878205, 24.844693, 24.811016, 24.777187, 24.743214, 24.709108, 24.674881, 24.640542, 24.606102, 24.571574, 24.536969, 24.502298, 24.467573, 24.432805, 24.398006, 24.363186, 24.328355, 24.293524, 24.258702, 24.223897, 24.189117, 24.154369, 24.119657, 24.084987, 24.050364, 24.015789, 23.981268, 23.946802, 23.912391, 23.878038, 23.843742, 23.809503, 23.775318, 23.741188, 23.707110, 23.673081, 23.639100, 23.605162, 23.571265, 23.537405, 23.503580, 23.469784, 23.436016, 23.402270, 23.368544, 23.334833, 23.301133, 23.267442, 23.233755, 23.200069, 23.166380, 23.132685, 23.098981, 23.065263, 23.031530, 22.997777, 22.964002, 22.930203, 22.896376, 22.862519, 22.828629, 22.794705, 22.760744, 22.726743, 22.692702, 22.658618, 22.624490, 22.590315, 22.556093, 22.521822, 22.487502, 22.453130, 22.418705, 22.384228, 22.349696, 22.315110, 22.280467, 22.245769, 22.211013, 22.176201, 22.141330, 22.106401, 22.071414, 22.036368, 22.001263, 21.966100, 21.930877, 21.895595, 21.860254, 21.824854, 21.789395, 21.753877, 21.718300, 21.682665, 21.646972, 21.611220, 21.575411, 21.539543, 21.503619, 21.467636, 21.431597, 21.395501, 21.359348, 21.323137, 21.286870, 21.250548, 21.214173, 21.177747, 21.141274, 21.104759, 21.068206, 21.031621, 20.995012, 20.958385, 20.921751, 20.885117, 20.848496, 20.811897, 20.775334, 20.738818, 20.702363, 20.665983, 20.629694, 20.593512, 20.557455, 20.521539, 20.485783, 20.450207, 20.414832, 20.379679, 20.344772, 20.310134, 20.275790, 20.241767, 20.208092, 20.174792, 20.141896, 20.109433, 20.077433, 20.045926, 20.014942, 19.984513, 19.954680, 19.925469, 19.896906, 19.869015, 19.841825, 19.815363, 19.789655, 19.764728, 19.740606, 19.717312, 19.694870, 19.673299, 19.652619, 19.632848, 19.614001, 19.596094, 19.579137, 19.563142, 19.548118, 19.534072, 19.521009, 19.508931, 19.497840, 19.487738, 19.478621, 19.470499, 19.463358, 19.457183, 19.451961, 19.447682, 19.444353, 19.441963, 19.440496, 19.439936, 19.440267, 19.441474, 19.443539, 19.446445, 19.450175, 19.454709, 19.460028, 19.466114, 19.472948, 19.480510, 19.488782, 19.497746, 19.507385, 19.517680, 19.528616, 19.540178, 19.552351, 19.565120, 19.578473, 19.592399, 19.606884, 19.621920, 19.637495, 19.653602, 19.670230, 19.687373, 19.705023, 19.723173, 19.741816, 19.760946, 19.780556, 19.800640, 19.821191, 19.842205, 19.863675, 19.885596, 19.907961, 19.930765, 19.954001, 19.977665, 20.001762, 20.026283, 20.051217, 20.076552, 20.102282, 20.128399, 20.154895, 20.181763, 20.208995, 20.236583, 20.264518, 20.292794, 20.321403, 20.350335, 20.379584, 20.409140, 20.438996, 20.469145, 20.499578, 20.530287, 20.561266, 20.592508, 20.624005, 20.655752, 20.687742, 20.719970, 20.752431, 20.785119, 20.818031, 20.851162, 20.884508, 20.918067, 20.951834, 20.985808, 21.019985, 21.054364, 21.088943, 21.123720, 21.158693, 21.193861, 21.229222, 21.264775, 21.300518, 21.336451, 21.372572, 21.408880, 21.445373, 21.482052, 21.518913, 21.555957, 21.593180, 21.630581, 21.668157, 21.705907, 21.743829, 21.781921, 21.820179, 21.858603, 21.897188, 21.935933, 21.974835, 22.013890, 22.053095, 22.092447, 22.131940, 22.171572, 22.211336, 22.251226, 22.291238, 22.331365, 22.371603, 22.411944, 22.452383, 22.492912, 22.533525, 22.574214, 22.614971, 22.655788, 22.696658, 22.737571, 22.778519, 22.819494, 22.860487, 22.901490, 22.942492, 22.983485, 23.024461, 23.065411, 23.106325, 23.147195, 23.188012, 23.228768, 23.269455, 23.310065, 23.350589, 23.391020, 23.431351, 23.471573, 23.511681, 23.551667, 23.591526, 23.631252, 23.670839, 23.710283, 23.749581, 23.788729, 23.827726, 23.866570, 23.905261, 23.943798, 23.982184, 24.020420, 24.058509, 24.096453, 24.134256, 24.171923, 24.209458, 24.246865, 24.284150, 24.321317, 24.358371, 24.395317, 24.432159, 24.468902, 24.505549, 24.542105, 24.578569, 24.614945, 24.651233, 24.687434, 24.723549, 24.759576, 24.795514, 24.831360, 24.867110, 24.902760, 24.938303, 24.973732, 25.009040, 25.044217, 25.079253, 25.114137, 25.148855, 25.183395, 25.217741, 25.251879, 25.285792, 25.319463, 25.352874, 25.386006, 25.418840, 25.451356, 25.483533, 25.515351, 25.546788, 25.577822, 25.608433, 25.638598, 25.668296, 25.697491, 25.726166, 25.754303, 25.781884, 25.808888, 25.835291, 25.861071, 25.886207, 25.910675, 25.934454, 25.957525, 25.979866, 26.001459, 26.022287, 26.042333, 26.061583, 26.080024, 26.097644, 26.114433, 26.130382, 26.145483, 26.159732, 26.173123, 26.185652, 26.197318, 26.208120, 26.218056, 26.227129, 26.235339, 26.242676, 26.249151, 26.254776, 26.259560, 26.263509, 26.266614, 26.268878, 26.270308, 26.270911, 26.270694, 26.269664, 26.267827, 26.265190, 26.261761, 26.257546, 26.252553, 26.246792, 26.240269, 26.232995, 26.224979, 26.216229, 26.206755, 26.196568, 26.185676, 26.174090, 26.161821, 26.148877, 26.135271, 26.121011, 26.106110, 26.090576, 26.074422, 26.057658, 26.040294, 26.022342, 26.003812, 25.984716, 25.965066, 25.944870, 25.924142, 25.902891, 25.881130, 25.858869, 25.836120, 25.812895, 25.789191, 25.765024, 25.740410, 25.715364, 25.689896, 25.664019, 25.637743, 25.611077, 25.584034, 25.556622, 25.528851, 25.500732, 25.472274, 25.443486, 25.414377, 25.384957, 25.355234, 25.325216, 25.294911, 25.264329, 25.233476, 25.202360, 25.170989, 25.139369, 25.107508, 25.075413, 25.043089, 25.010543, 24.977781, 24.944809, 24.911632, 24.878255, 24.844684, 24.810923, 24.776976, 24.742849, 24.708543, 24.674065, 24.639417, 24.604602, 24.569624, 24.534486, 24.499191, 24.463742, 24.428140, 24.392390, 24.356494, 24.320454, 24.284274, 24.247957, 24.211505, 24.174922, 24.138211, 24.101377, 24.064423, 24.027352, 23.990170, 23.952881, 23.915489, 23.877999, 23.840415, 23.802743, 23.764988, 23.727155, 23.689248, 23.651272, 23.613232, 23.575133, 23.536979, 23.498773, 23.460520, 23.422222, 23.383884, 23.345507, 23.307093, 23.268645, 23.230164, 23.191651, 23.153107, 23.114532, 23.075927, 23.037291, 22.998625, 22.959927, 22.921197, 22.882435, 22.843639, 22.804810, 22.765945, 22.727046, 22.688112, 22.649143, 22.610138, 22.571100, 22.532027, 22.492921, 22.453784, 22.414617, 22.375422, 22.336199, 22.296953, 22.257685, 22.218398, 22.179094, 22.139776, 22.100447, 22.061111, 22.021771, 21.982430, 21.943091, 21.903757, 21.864434, 21.825123, 21.785828, 21.746554, 21.707303, 21.668079, 21.628886, 21.589728, 21.550607, 21.511527, 21.472492, 21.433504, 21.394568, 21.355685, 21.316859, 21.278091, 21.239384, 21.200739, 21.162158, 21.123641, 21.085189, 21.046800, 21.008475, 20.970213, 20.932012, 20.893871, 20.855789, 20.817764, 20.779798, 20.741888, 20.704036, 20.666243, 20.628510, 20.590839, 20.553234, 20.515697, 20.478234, 20.440849, 20.403549, 20.366338, 20.329224, 20.292215, 20.255317, 20.218541, 20.181894, 20.145390, 20.109041, 20.072859, 20.036859, 20.001058, 19.965473, 19.930120, 19.895020, 19.860193, 19.825661, 19.791446, 19.757573, 19.724065, 19.690949, 19.658252, 19.626000, 19.594220, 19.562943, 19.532195, 19.502019, 19.472438, 19.443473, 19.415148, 19.387490, 19.360523, 19.334273, 19.308764, 19.284019, 19.260061, 19.236910, 19.214585, 19.193106, 19.172487, 19.152745, 19.133892, 19.115941, 19.098902, 19.082783, 19.067590, 19.053331, 19.040007, 19.027623, 19.016178, 19.005672, 18.996103, 18.987467, 18.979769, 18.972995, 18.967128, 18.962155, 18.958062, 18.954854, 18.952516, 18.951032, 18.950384, 18.950555, 18.951525, 18.953278, 18.955793, 18.959051, 18.963034, 18.967720, 18.973092, 18.979131, 18.985819, 18.993138, 19.001072, 19.009606, 19.018725, 19.028414, 19.038661, 19.049454, 19.060782, 19.072635, 19.085002, 19.097876, 19.111249, 19.125114, 19.139464, 19.154295, 19.169599, 19.185373, 19.201612, 19.218311, 19.235466, 19.253074, 19.271131, 19.289634, 19.308579, 19.327964, 19.347785, 19.368039, 19.388725, 19.409840, 19.431381, 19.453347, 19.475736, 19.498545, 19.521774, 19.545421, 19.569497, 19.593999, 19.618921, 19.644259, 19.670012, 19.696180, 19.722762, 19.749760, 19.777172, 19.804998, 19.833239, 19.861892, 19.890957, 19.920433, 19.950316, 19.980605, 20.011296, 20.042387, 20.073872, 20.105748, 20.138008, 20.170647, 20.203657, 20.237032, 20.270764, 20.304844, 20.339266, 20.374020, 20.409097, 20.444489, 20.480186, 20.516180, 20.552461, 20.589020, 20.625847, 20.662935, 20.700273, 20.737852, 20.775663, 20.813697, 20.851946, 20.890401, 20.929055, 20.967900, 21.006929, 21.046133, 21.085506, 21.125040, 21.164728, 21.204564, 21.244540, 21.284649, 21.324886, 21.365243, 21.405714, 21.446292, 21.486971, 21.527745, 21.568607, 21.609551, 21.650571, 21.691660, 21.732813, 21.774023, 21.815284, 21.856592, 21.897938, 21.939319, 21.980727, 22.022157, 22.063602, 22.105058, 22.146517, 22.187973, 22.229420, 22.270852, 22.312260, 22.353641, 22.394987, 22.436290, 22.477545, 22.518743, 22.559876, 22.600937, 22.641919, 22.682812, 22.723610, 22.764305, 22.804889, 22.845354, 22.885693, 22.925898, 22.965964, 23.005883, 23.045651, 23.085262, 23.124711, 23.163994, 23.203107, 23.242047, 23.280814, 23.319409, 23.357832, 23.396085, 23.434171, 23.472093, 23.509854, 23.547459, 23.584911, 23.622216, 23.659378, 23.696403, 23.733296, 23.770062, 23.806705, 23.843230, 23.879641, 23.915942, 23.952135, 23.988223, 24.024207, 24.060088, 24.095866, 24.131541, 24.167111, 24.202572, 24.237922, 24.273156, 24.308267, 24.343250, 24.378095, 24.412795, 24.447340, 24.481717, 24.515915, 24.549922, 24.583722, 24.617302, 24.650647, 24.683741, 24.716567, 24.749106, 24.781341, 24.813252, 24.844821, 24.876027, 24.906850, 24.937269, 24.967262, 24.996795, 25.025851, 25.054412, 25.082461, 25.109976, 25.136935, 25.163315, 25.189095, 25.214252, 25.238764, 25.262611, 25.285772, 25.308229, 25.329964, 25.350962, 25.371207, 25.390686, 25.409388, 25.427303, 25.444423, 25.460742, 25.476254, 25.490958, 25.504852, 25.517936, 25.530213, 25.541685, 25.552358, 25.562238, 25.571331, 25.579647, 25.587194, 25.593970, 25.599994, 25.605286, 25.609862, 25.613737, 25.616923, 25.619435, 25.621270, 25.622442, 25.622969, 25.622866, 25.622148, 25.620826, 25.618914, 25.616422, 25.613363, 25.609747, 25.605583, 25.600882, 25.595651, 25.589899, 25.583634, 25.576863, 25.569591, 25.561825, 25.553570, 25.544831, 25.535614, 25.525922, 25.515759, 25.505131, 25.494040, 25.482492, 25.470491, 25.458041, 25.445146, 25.431811, 25.418042, 25.403844, 25.389221, 25.374180, 25.358726, 25.342866, 25.326605, 25.309949, 25.292906, 25.275481, 25.257680, 25.239509, 25.220973, 25.202078, 25.182828, 25.163227, 25.143278, 25.122983, 25.102343, 25.081360, 25.060034, 25.038363, 25.016346, 24.993980, 24.971262, 24.948190, 24.924757, 24.900960, 24.876794, 24.852240, 24.827297, 24.801966, 24.776247, 24.750138, 24.723636, 24.696742, 24.669456, 24.641781, 24.613720, 24.585281, 24.556470, 24.527297, 24.497772, 24.467909, 24.437721, 24.407224, 24.376434, 24.345369, 24.314048, 24.282489, 24.250713, 24.218739, 24.186587, 24.154277, 24.121826, 24.089253, 24.056574, 24.023806, 23.990963, 23.958059, 23.925106, 23.892115, 23.859095, 23.826053, 23.792996, 23.759928, 23.726854, 23.693775, 23.660695, 23.627612, 23.594528, 23.561442, 23.528352, 23.495257, 23.462154, 23.429041, 23.395916, 23.362777, 23.329620, 23.296444, 23.263245, 23.230021, 23.196770, 23.163491, 23.130184, 23.096845, 23.063476, 23.030074, 22.996640, 22.963174, 22.929675, 22.896144, 22.862581, 22.828987, 22.795363, 22.761709, 22.728028, 22.694321, 22.660589, 22.626835, 22.593060, 22.559267, 22.525458, 22.491635, 22.457801, 22.423959, 22.390111, 22.356261, 22.322411, 22.288564, 22.254722, 22.220890, 22.187068, 22.153261, 22.119470, 22.085697, 22.051945, 22.018217, 21.984512, 21.950834, 21.917183, 21.883561, 21.849968, 21.816406, 21.782874, 21.749374, 21.715905, 21.682467, 21.649061, 21.615685, 21.582340, 21.549026, 21.515740, 21.482483, 21.449253, 21.416050, 21.382871, 21.349716, 21.316584, 21.283471, 21.250378, 21.217302, 21.184242, 21.151198, 21.118169, 21.085156, 21.052159, 21.019180, 20.986223, 20.953290, 20.920388, 20.887520, 20.854693, 20.821914, 20.789191, 20.756532, 20.723946, 20.691442, 20.659031, 20.626722, 20.594528, 20.562460, 20.530531, 20.498755, 20.467147, 20.435723, 20.404499, 20.373494, 20.342724, 20.312212, 20.281976, 20.252038, 20.222421, 20.193148, 20.164256, 20.135766, 20.107696, 20.080069, 20.052911, 20.026247, 20.000103, 19.974505, 19.949478, 19.925047, 19.901237, 19.878070, 19.855571, 19.833760, 19.812658, 19.792286, 19.772661, 19.753802, 19.735725, 19.718443, 19.701971, 19.686319, 19.671499, 19.657519, 19.644385, 19.632104, 19.620679, 19.610112, 19.600405, 19.591557, 19.583566, 19.576440, 19.570166, 19.564730, 19.560120, 19.556328, 19.553356, 19.551195, 19.549830, 19.549246, 19.549430, 19.550365, 19.552035, 19.554426, 19.557519, 19.561298, 19.565745, 19.570844, 19.576576, 19.582927, 19.589879, 19.597419, 19.605531, 19.614202, 19.623420, 19.633172, 19.643449, 19.654239, 19.665536, 19.677329, 19.689612, 19.702379, 19.715622, 19.729338, 19.743521, 19.758167, 19.773271, 19.788832, 19.804845, 19.821308, 19.838218, 19.855574, 19.873373, 19.891615, 19.910296, 19.929417, 19.948975, 19.968969, 19.989398, 20.010261, 20.031555, 20.053280, 20.075433, 20.098013, 20.121018, 20.144445, 20.168291, 20.192568, 20.217269, 20.242385, 20.267911, 20.293843, 20.320178, 20.346913, 20.374044, 20.401569, 20.429485, 20.457788, 20.486476, 20.515545, 20.544992, 20.574815, 20.605009, 20.635572, 20.666500, 20.697789, 20.729437, 20.761440, 20.793794, 20.826495, 20.859540, 20.892926, 20.926649, 20.960706, 20.995094, 21.029810, 21.064852, 21.100217, 21.135901, 21.171904, 21.208221, 21.244849, 21.281786, 21.319027, 21.356569, 21.394407, 21.432537, 21.470952, 21.509649, 21.548620, 21.587859, 21.627358, 21.667111, 21.707110, 21.747345, 21.787808, 21.828491, 21.869383, 21.910476, 21.951760, 21.993224, 22.034861, 22.076659, 22.118610, 22.160704, 22.202933, 22.245286, 22.287755, 22.330330, 22.373002, 22.415763, 22.458601, 22.501509, 22.544476, 22.587493, 22.630550, 22.673637, 22.716745, 22.759862, 22.802980, 22.846088, 22.889177, 22.932236, 22.975256, 23.018227, 23.061140, 23.103985, 23.146753, 23.189435, 23.232023, 23.274508, 23.316882, 23.359136, 23.401264, 23.443256, 23.485106, 23.526805, 23.568347, 23.609724, 23.650927, 23.691951, 23.732788, 23.773429, 23.813869, 23.854099, 23.894113, 23.933902, 23.973460, 24.012780, 24.051854, 24.090675, 24.129236, 24.167530, 24.205548, 24.243285, 24.280732, 24.317883, 24.354730, 24.391265, 24.427481, 24.463370, 24.498924, 24.534136, 24.568998, 24.603499, 24.637633, 24.671391, 24.704762, 24.737738, 24.770310, 24.802467, 24.834200, 24.865497, 24.896350, 24.926748, 24.956680, 24.986123, 25.015071, 25.043518, 25.071459, 25.098884, 25.125783, 25.152146, 25.177965, 25.203230, 25.227932, 25.252062, 25.275614, 25.298578, 25.320948, 25.342717, 25.363879, 25.384428, 25.404360, 25.423668, 25.442349, 25.460400, 25.477817, 25.494596, 25.510737, 25.526237, 25.541094, 25.555306, 25.568874, 25.581796, 25.594072, 25.605703, 25.616688, 25.627029, 25.636727, 25.645785, 25.654206, 25.661992, 25.669146, 25.675673, 25.681564, 25.686832, 25.691487, 25.695541, 25.699000, 25.701869, 25.704155, 25.705861, 25.706979, 25.707516, 25.707481, 25.706883, 25.705729, 25.704025, 25.701777, 25.698990, 25.695670, 25.691822, 25.687452, 25.682565, 25.677167, 25.671262, 25.664856, 25.657953, 25.650559, 25.642678, 25.634314, 25.625474, 25.616161, 25.606380, 25.596136, 25.585434, 25.574279, 25.562674, 25.550625, 25.538137, 25.525212, 25.511857, 25.498075, 25.483872, 25.469250, 25.454215, 25.438771, 25.422921, 25.406672, 25.390026, 25.372988, 25.355561, 25.337751, 25.319562, 25.300996, 25.282058, 25.262751, 25.243081, 25.223050, 25.202663, 25.181924, 25.160836, 25.139404, 25.117631, 25.095522, 25.073079, 25.050306, 25.027208, 25.003789, 24.980050, 24.955998, 24.931634, 24.906951, 24.881955];
        // Set initial temperature
        simple_model.spaces[0].set_dry_bulb_temperature(&mut state, exp_zone_air_temp[0]);

        let mut date = Date{
            month: 1, day: 1, hour: 0.0,
        };
        for i in 0..outdoor_temp.len() {
            // Get zone's temp
            let found_temp = simple_model.spaces[0].dry_bulb_temperature(&state).unwrap();
            let exp_temp = exp_zone_air_temp[i];
            if i > 450{
                // ignore warmup period
                println!("{},{}", found_temp, exp_temp);
                assert!( (found_temp-exp_temp).abs() < 0.4, "found_temp = {}, exp_temp = {} ,  error = {}", found_temp, exp_temp, (found_temp-exp_temp).abs() )
            }
            
            
            // Set outdoor temp
            let mut weather = SyntheticWeather::new();
            weather.dry_bulb_temperature = Box::new(ScheduleConstant::new(outdoor_temp[i]));

            let surface = &simple_model.surfaces[0];

            // Set Solar Radiation
            surface.set_back_incident_solar_irradiance(&mut state, incident_solar_radiation[i]);

            // Set IR radiation
            surface.set_back_ir_irradiance(&mut state, thermal_heat_gain[i]/surface_area/0.9);

            // March
            thermal_model.march(date, &weather, &simple_model, &mut state).unwrap();

            
            

            // Advance
            date.add_hours(1./n as Float);
            // assert!(false)

            
        }
        
        
    }
}
