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
use calendar::date::Date;
use communication_protocols::error_handling::ErrorHandling;
use communication_protocols::simulation_model::SimulationModel;
use weather::Weather;

use crate::surface::ThermalSurface;
use simple_model::{
    SimpleModel, 
    SimulationState, SimulationStateHeader,    
    Boundary,
};
// use simple_model::simulation_state_element::SimulationStateElement;
// use simple_model::hvac::*;
// use simple_model::hvac::ideal_heater_cooler::IdealHeaterCooler;
// use simple_model::hvac::electric_heater::ElectricHeater;
// use simple_model::boundary::Boundary;

use crate::zone::ThermalZone;
use crate::heating_cooling::calc_cooling_heating_power;
use crate::construction::discretize_construction;

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
    fn new(model: &SimpleModel, state: &mut SimulationStateHeader, n: usize) -> Result<Self, String> {
        

        /* CREATE ALL ZONES, ONE PER SPACE */        
        let mut thermal_zones: Vec<ThermalZone> = Vec::with_capacity(model.spaces.len());
        for (i,space) in model.spaces.iter().enumerate() {
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
        let dt = 60. * 60. / (n as Float * n_subdivisions as Float);

        if n * n_subdivisions < 6 {
            eprintln!("Number of timesteps per hour (n) is too small in Finite Difference Thermal  Module... try to use 6 or more.");
        }

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
            if let Ok(b) = surface.front_boundary(){
                thermal_surfaces[i].set_front_boundary(b.clone());
            }
            if let Ok(b)=surface.back_boundary(){
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
            if let Ok(b) = fenestration.front_boundary(){
                thermal_fenestrations[i].set_front_boundary(b.clone());
            }
            if let Ok(b) = fenestration.back_boundary(){
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
        for _ in 0..self.dt_subdivisions {
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
                    Some(b)=> match b {
                            Boundary::Space(space) => t_current[*space.index().unwrap()], 
                            Boundary::Ground => unimplemented!(),
                    },                                
                    None => t_out    
                };
                let t_back = match s.back_boundary() {
                    Some(b)=> match b {
                        Boundary::Space(space) => t_current[*space.index().unwrap()], //self.zones[z_index].temperature(model, state),
                        Boundary::Ground => unimplemented!(),
                    },
                    None => t_out,                    
                };

                // Update temperatures
                let (_q_front, _q_back) = s.march( state, t_front, t_back, self.dt);
            } // end of iterating surface

            // What  if they are open???
            for i in 0..self.fenestrations.len() {
                // get surface
                let s = &self.fenestrations[i];

                // find t_in and t_out of surface.
                let t_front = match s.front_boundary() {
                    Some(b)=>match b{
                        Boundary::Space(space) => t_current[*space.index().unwrap()],
                        Boundary::Ground => unimplemented!(),
                    },
                    None => t_out                    
                };
                let t_back = match s.back_boundary() {
                    Some(b)=>match b{
                        Boundary::Space(space) => t_current[*space.index().unwrap()],
                        Boundary::Ground => unimplemented!(),
                    },
                    None => t_out                    
                };

                // Update temperatures
                s.march( state, t_front, t_back, self.dt);
            } // end of iterating surface

            /* UPDATE ZONES' TEMPERATURE */
            // This is done analytically.
            let (a, b, c) = self.calculate_zones_abc(model, state);
            let future_temperatures =
                self.estimate_zones_future_temperatures(&t_current, &a, &b, &c, self.dt);
            for (i, zone) in self.zones.iter().enumerate() {
                zone.reference_space.set_dry_bulb_temperature(state, future_temperatures[i]);
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
        for hvac in model.hvacs.iter(){
            for (target_space_index, heating_cooling) in calc_cooling_heating_power(hvac, state){
                a[target_space_index] += heating_cooling;
            }
            // heating through air supply?
        }

        // Heating/Cooling
        for luminaire in model.luminaires.iter(){
            if let Ok(target_space) = luminaire.target_space() {                
                let target_space_index = *target_space.index().unwrap();
                let consumption = luminaire.power_consumption(state).expect("Luminaire has no Power Consumption state");                
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
            if let Some(t_inf_inwards) = space.infiltration_temperature(state){
                let m_inf = space.infiltration_volume(state).expect("Space has infiltration temperature but not volume");            
                let cp_inf_inwards = gas_properties::air::specific_heat();
                a[i] += m_inf * cp_inf_inwards * t_inf_inwards;
                b[i] += m_inf * cp_zone;
            }

            
            // ventilation
            if let Some(t_vent_inwards) = space.ventilation_temperature(state){
                let m_vent = space.ventilation_volume(state).expect("Space has ventilation temperature but not volume");            
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
                let ai = surface.area();
                // if front leads to a Zone
                if let Some(Boundary::Space(space)) = surface.front_boundary() {
                    let z_index = space.index().unwrap();
                    let hi = 1. / surface.rs_front();
                    let temp = surface.front_temperature(state);
                    a[*z_index] += hi * ai * temp;
                    b[*z_index] += hi * ai;
                }

                // if back leads to a Zone
                if let Some(Boundary::Space(space)) = surface.back_boundary() {
                    let z_index = space.index().unwrap();
                    let hi = 1. / surface.rs_back();
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
    fn get_current_zones_temperatures(
        &self,        
        state: &SimulationState,
    ) -> Vec<Float> {
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

            ret.push(
                a[i] / b[i]
                    + (c[i] * (current_temp - a[i] / b[i]) / future_time / b[i])
                        * (1.0 - (-b[i] * future_time / c[i]).exp()),
            );
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
            ret.push(
                a[i] / b[i] + (t_current[i] - a[i] / b[i]) * (-b[i] * future_time / c[i]).exp(),
            );
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

    use calendar::date::Date;
    use schedule::constant::ScheduleConstant;
    use weather::synthetic_weather::SyntheticWeather;
    
    use simple_test_models::*;
    use gas_properties::air;
    use simple_model::{
        HVAC, 
        SimulationStateElement
    };
    
    /// A single-zone test model with walls assumed to have
    /// no mass. It has a closed solution, which is nice.
    /// 
    /// There is no sun.
    #[derive(Default)]
    struct SingleZoneTestModel{
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
        fn get_closed_solution(&self)->Box<impl Fn(Float)->Float>{
            // heat balance in the form
            // of C*dT/dt = A - B*T
            let rho = air::density(); //kg/m3
            let cp = air::specific_heat(); //J/kg.K
            let u = 1./self.facade_r;

            let c = self.zone_volume * rho * cp;

            let mut a =     
                    self.heating_power 
                    + self.lighting_power 
                    + self.temp_out * u * self.surface_area + self.infiltration_rate * cp * self.temp_out ;

            a/=c;

            let mut b = u * self.surface_area + self.infiltration_rate*cp;
            b /= c;

            let k1 = self.temp_start - a/b;

            let f = move |t: Float|->Float{
                a/b + k1*(-b*t).exp()
            };

            Box::new(f)            
        }
    }



    #[test]
    fn test_calculate_zones_abc() {
        
        let  (simple_model, mut state_header) = get_single_zone_test_building(
            // &mut state,
            &SingleZoneTestBuildingOptions {
                zone_volume: 40.,
                surface_area: 4.,
                material_is_massive: Some(false),
                ..Default::default()
            },
        );

        let n: usize = 1;
        let thermal_model = ThermalModel::new(& simple_model, &mut state_header, n).unwrap();
        let state = state_header.take_values().unwrap();
        // MAP THE STATE
        // model.map_simulation_state(&mut state).unwrap();

        // Test
        let (a, b, c) = thermal_model.calculate_zones_abc(&simple_model, &state);
        assert_eq!(a.len(), 1);
        assert_eq!(c.len(), 1);
        assert_eq!(b.len(), 1);
        assert_eq!(c[0], thermal_model.get_thermal_zone(0).unwrap().mcp());
        let hi = 1. / thermal_model.get_thermal_surface(0).unwrap().rs_front();
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
        // assert!(model.surfaces[0].is_massive());

        let r = construction.r_value().unwrap()
            + thermal_model.surfaces[0].rs_front()
            + thermal_model.surfaces[0].rs_back();


        // Initial T of the zone
        let t_start = thermal_model.zones[0].reference_space.dry_bulb_temperature(&state).unwrap(); 

        let t_out: Float = 30.0; // T of surroundings

        // test model
        let tester = SingleZoneTestModel{
            zone_volume,
            surface_area,
            facade_r: r,
            temp_out: t_out,
            temp_start: t_start,
            .. SingleZoneTestModel::default()
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

            let found = thermal_model.zones[0].reference_space.dry_bulb_temperature(&state).unwrap();

            thermal_model.march(date, &weather, &simple_model, &mut state).unwrap();

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

        
        let  (simple_model, mut state_header) = get_single_zone_test_building(
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

        let r = construction.r_value().unwrap()
            + thermal_model.surfaces[0].rs_front()
            + thermal_model.surfaces[0].rs_back();
        
        // Initial T of the zone
        let t_start = thermal_model.zones[0].reference_space.dry_bulb_temperature(&state).unwrap(); 

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
         let tester = SingleZoneTestModel{
            zone_volume,
            surface_area, // the window is a hole on the wall... does not add area
            facade_r: r,
            temp_out: t_out,
            temp_start: t_start,
            .. SingleZoneTestModel::default()
        };
        let exp_fn = tester.get_closed_solution();

        // March:
        for i in 0..80 {
            let time = (i as Float) * dt;
            date.add_seconds(time);

            let found = thermal_model.zones[0].reference_space.dry_bulb_temperature(&state).unwrap();

            thermal_model.march(date, &weather, &simple_model, &mut state).unwrap();

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
        let lum_state_i = simple_model.luminaires[0].power_consumption_index().unwrap();
        // state.update_value(lum_state_i, SimulationStateElement::LuminairePowerConsumption(0, lighting_power));
        state[lum_state_i] = lighting_power;

        // START TESTING.
        let construction = &simple_model.constructions[0];
        // assert!(!model.surfaces[0].is_massive());

        let r = construction.r_value().unwrap()
            + thermal_model.surfaces[0].rs_front()
            + thermal_model.surfaces[0].rs_back();



        // Initial T of the zone
        let t_start = thermal_model.zones[0].reference_space.dry_bulb_temperature(&state).unwrap(); 

        let t_out: Float = 30.0; // T of surroundings

        // test model
        let tester = SingleZoneTestModel{
            zone_volume,
            surface_area, // the window is a hole on the wall... does not add area
            lighting_power,
            facade_r: r,
            temp_out: t_out,
            temp_start: t_start,
            .. SingleZoneTestModel::default()
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

            let found = thermal_model.zones[0].reference_space.dry_bulb_temperature(&state).unwrap();

            thermal_model.march(date, &weather, &simple_model, &mut state).unwrap();

            // Get exact solution.
            let exp = exp_fn(time);

            let max_error = 0.55;            
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
        if let HVAC::ElectricHeater(heater) = &simple_model.hvacs[0]{
            let hvac_state_i = heater.heating_cooling_consumption_index().unwrap();
            state[hvac_state_i] = heating_power;
        }
        

        // START TESTING.
        let construction = &simple_model.constructions[0];
        // assert!(!model.surfaces[0].is_massive());

        let r = construction.r_value().unwrap()
            + thermal_model.surfaces[0].rs_front()
            + thermal_model.surfaces[0].rs_back();



        // Initial T of the zone
        let t_start = thermal_model.zones[0].reference_space.dry_bulb_temperature(&state).unwrap(); 
        let t_out: Float = 30.0; // T of surroundings

        // test model
        let tester = SingleZoneTestModel{
            zone_volume,
            surface_area, // the window is a hole on the wall... does not add area
            heating_power,
            facade_r: r,
            temp_out: t_out,
            temp_start: t_start,
            .. SingleZoneTestModel::default()
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

            let found = thermal_model.zones[0].reference_space.dry_bulb_temperature(&state).unwrap();

            thermal_model.march(date, &weather, &simple_model, &mut state).unwrap();

            // Get exact solution.
            let exp = exp_fn(time);

            let max_error = 0.55;            
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
        let inf_vol_index = state_header.push(SimulationStateElement::SpaceInfiltrationVolume(0), infiltration_rate);
        simple_model.spaces[0].set_infiltration_volume_index(inf_vol_index);
        let inf_temp_index = state_header.push(SimulationStateElement::SpaceInfiltrationTemperature(0), t_out);
        simple_model.spaces[0].set_infiltration_temperature_index(inf_temp_index);

        // MAP THE STATE
        
        let mut state = state_header.take_values().unwrap();

        // turn the heater on
        if let HVAC::ElectricHeater(heater) = &simple_model.hvacs[0]{
            let hvac_state_i = heater.heating_cooling_consumption_index().unwrap();
            state[hvac_state_i] = heating_power;
        }
    
        // START TESTING.
        let construction = &simple_model.constructions[0];
        // assert!(!model.surfaces[0].is_massive());

        let r = construction.r_value().unwrap()
            + thermal_model.surfaces[0].rs_front()
            + thermal_model.surfaces[0].rs_back();



        // Initial T of the zone
        let t_start = thermal_model.zones[0].reference_space.dry_bulb_temperature(&state).unwrap(); 
        

        // test model
        let tester = SingleZoneTestModel{
            zone_volume,
            surface_area, // the window is a hole on the wall... does not add area
            heating_power,
            facade_r: r,
            temp_out: t_out,
            temp_start: t_start,
            infiltration_rate,
            .. SingleZoneTestModel::default()
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

            let found = thermal_model.zones[0].reference_space.dry_bulb_temperature(&state).unwrap();

            thermal_model.march(date, &weather, &simple_model, &mut state).unwrap();

            // Get exact solution.
            let exp = exp_fn(time);

            let max_error = 0.55;            
            println!("{}, {}", exp, found);
            assert!((exp - found).abs() < max_error);
            
        }
    }
    
}
