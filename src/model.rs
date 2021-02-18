use crate::surface::*;
use crate::zone::ThermalZone;
use building_model::boundary::Boundary;
use building_model::building::Building;
use building_model::object_trait::ObjectTrait;
use calendar::date::Date;
use communication_protocols::error_handling::ErrorHandling;
use communication_protocols::simulation_model::SimulationModel;
use simulation_state::simulation_state::SimulationState;
use weather::Weather;

use crate::construction::discretize_construction;

pub struct ThermalModel {
    zones: Vec<ThermalZone>,
    surfaces: Vec<ThermalSurface>,
    fenestrations: Vec<ThermalSurface>,
    dt_subdivisions: usize,
    dt: f64,
}

impl ErrorHandling for ThermalModel {
    fn module_name() -> &'static str {
        "Finite Difference Thermal Model"
    }
}

impl SimulationModel for ThermalModel {
    type Type = Self;

    fn new(building: &Building, state: &mut SimulationState, n: usize) -> Result<Self, String> {
        let main_dt = 60. * 60. / n as f64;

        /* CREATE ALL ZONES, ONE PER SPACE */
        let spaces = building.get_spaces();
        let mut thermal_zones: Vec<ThermalZone> = Vec::with_capacity(spaces.len());
        for i in 0..spaces.len() {
            // Add the zone to the model... this pushes it to the sate
            // as well
            thermal_zones.push(ThermalZone::from_space(&spaces[i], state));
        }

        /* FIND MODEL TIMESTEP */
        // choose the smallest timestep in all constructions
        let max_dx = 0.04; // 4cm
        let min_dt = 60.; // 1 minute

        let mut n_subdivisions: usize = 1;
        let constructions = building.get_constructions();

        // Store the dts and n_nodes somwehere
        let mut all_n_elements: Vec<Vec<usize>> = Vec::with_capacity(constructions.len());
        for construction in constructions {
            let (found_n_subdivisions, n_elements) =
                discretize_construction(building, construction, main_dt, max_dx, min_dt);
            if found_n_subdivisions > n_subdivisions {
                n_subdivisions = found_n_subdivisions;
            }
            all_n_elements.push(n_elements);
        }

        let dt = main_dt / (n_subdivisions as f64);

        /* CREATE SURFACES USING THE MINIMUM TIMESTEP */
        // The rationale here is the following: We find the minimum
        // timestep (or maximum timestep_subdivisions), and that will be the
        // dt for the whole model. Some constructions will have a larger
        // dt (due to their discretization scheme), but simulating them
        // with a smaller (i.e. the official) dt is of no harm.
        let surfaces = building.get_surfaces();

        // For the Thermal Model
        let mut thermal_surfaces: Vec<ThermalSurface> = Vec::with_capacity(surfaces.len());

        for i in 0..surfaces.len() {
            let surface = &surfaces[i];
            let construction_index = surface.get_construction_index().unwrap();

            let thermal_surface = match ThermalSurface::from_surface(
                building,
                state,
                surface,
                dt,
                &all_n_elements[construction_index],
                surface.index(),
            ) {
                Ok(v) => v,
                Err(e) => return Err(e),
            };

            thermal_surfaces.push(thermal_surface);

            // Match surface and zones
            thermal_surfaces[i].set_front_boundary(*surface.front_boundary());
            thermal_surfaces[i].set_back_boundary(*surface.back_boundary());
        }

        let fenestrations = building.get_fenestrations();
        let mut thermal_fenestrations: Vec<ThermalSurface> =
            Vec::with_capacity(fenestrations.len());
        for i in 0..fenestrations.len() {
            let construction_index = fenestrations[i].get_construction_index().unwrap();

            let i = thermal_fenestrations.len();
            let thermal_surface = match ThermalSurface::from_fenestration(
                building,
                state,
                &fenestrations[i],
                dt,
                &all_n_elements[construction_index],
                fenestrations[i].index(),
            ) {
                Ok(v) => v,
                Err(e) => return Err(e),
            };

            thermal_fenestrations.push(thermal_surface);

            // Match surface and zones
            thermal_fenestrations[i].set_front_boundary(*fenestrations[i].front_boundary());
            thermal_fenestrations[i].set_back_boundary(*fenestrations[i].back_boundary());
        }

        let ret = ThermalModel {
            zones: thermal_zones,
            surfaces: thermal_surfaces,
            fenestrations: thermal_fenestrations,
            dt_subdivisions: n_subdivisions,
            dt: dt,
        };
        return Ok(ret);
    }

    /* ********************************* */
    /* ********************************* */
    /* ********************************* */
    /* ********************************* */
    /* ********************************* */

    fn march(
        &self,
        date: Date,
        weather: &dyn Weather,
        building: &Building,
        state: &mut SimulationState,
    ) -> Result<(), String> {
        let mut date = date.clone();

        // Iterate through all the sub-subdivitions
        for _ in 0..self.dt_subdivisions {
            // advance in time
            date.add_seconds(self.dt);
            let current_weather = weather.get_weather_data(date);

            let t_out = match current_weather.dry_bulb_temperature {
                Some(v) => v,
                None => {
                    return Err(format!(
                    "Trying to march on Thermal Model, but dry bulb temperature was not provided"
                ))
                }
            };

            // update surface temperatures
            let mut heat_storage: Vec<f64> = vec![0.0; self.zones.len()];

            for i in 0..self.surfaces.len() {
                // get surface
                let s = &self.surfaces[i];

                // find t_in and t_out of surface.
                let t_front = match s.front_boundary() {
                    Boundary::Space(z_index) => self.zones[z_index].temperature(building, state),
                    Boundary::Ground => unimplemented!(),
                    Boundary::None => t_out,
                };
                let t_back = match s.back_boundary() {
                    Boundary::Space(z_index) => self.zones[z_index].temperature(building, state),
                    Boundary::Ground => unimplemented!(),
                    Boundary::None => t_out,
                };

                // Update temperatures
                let (q_front, q_back) = s.march(building, state, t_front, t_back);

                // Distribute heat flows.
                match s.front_boundary() {
                    Boundary::Space(z_index) => {
                        heat_storage[z_index] += q_front * s.area() * self.dt
                    }
                    Boundary::Ground | Boundary::None => {}
                };
                match s.back_boundary() {
                    Boundary::Space(z_index) => {
                        heat_storage[z_index] += q_back * s.area() * self.dt
                    }
                    Boundary::Ground | Boundary::None => {}
                };
            } // end of iterating surface

            for i in 0..self.fenestrations.len() {
                // get surface
                let s = &self.fenestrations[i];

                // find t_in and t_out of surface.
                let t_front = match s.front_boundary() {
                    Boundary::Space(z_index) => self.zones[z_index].temperature(building, state),
                    Boundary::Ground => unimplemented!(),
                    Boundary::None => t_out,
                };
                let t_back = match s.back_boundary() {
                    Boundary::Space(z_index) => self.zones[z_index].temperature(building, state),
                    Boundary::Ground => unimplemented!(),
                    Boundary::None => t_out,
                };

                // Update temperatures
                let (q_front, q_back) = s.march(building, state, t_front, t_back);

                // Distribute heat flows.
                match s.front_boundary() {
                    Boundary::Space(z_index) => {
                        heat_storage[z_index] += q_front * s.area() * self.dt
                    }
                    Boundary::Ground | Boundary::None => {}
                };
                match s.back_boundary() {
                    Boundary::Space(z_index) => {
                        heat_storage[z_index] += q_back * s.area() * self.dt
                    }
                    Boundary::Ground | Boundary::None => {}
                };
            } // end of iterating surface

            // update air flows
            // assume constant during timestep... get a vector with Qs

            // ZONE:
            for i in 0..self.zones.len() {
                // calculate air-flow heat transfer

                // calculate infiltration

                // calculate Zone heating/cooling
                let heating_cooling_power =
                    self.zones[i].calc_heating_cooling_power(building, state);
                heat_storage[i] += heating_cooling_power * self.dt;

                // Calculate people

                // Calculate lighting
                heat_storage[i] += self.zones[i].calc_lighting_power(building, state) * self.dt;

                // update all zones temperatures
                self.zones[i].consume_heat(heat_storage[i], building, state);
            } // end of iterating zones
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
}

/***********/
/* TESTING */
/***********/

#[cfg(test)]
mod testing {
    use super::*;
    use geometry3d::loop3d::Loop3D;
    use geometry3d::point3d::Point3D;
    use geometry3d::polygon3d::Polygon3D;

    use building_model::building::Building;
    use building_model::fenestration::*;
    use building_model::material::MaterialProperties;
    use building_model::substance::SubstanceProperties;

    use crate::construction::*;
    use building_model::boundary::Boundary;

    use calendar::date::Date;
    use schedule::constant::ScheduleConstant;
    use weather::synthetic_weather::SyntheticWeather;

    #[test]
    fn test_model_march() {
        let mut building = Building::new("the building".to_string());
        let mut state = SimulationState::new();

        // Add the space
        let zone_volume = 40.;
        let space_index = building.add_space("Some space".to_string());
        building.set_space_volume(space_index, zone_volume).unwrap();

        // Add substance
        let poly_index = building.add_substance("polyurethane".to_string());
        building
            .set_substance_properties(
                poly_index,
                SubstanceProperties {
                    thermal_conductivity: 0.0252,  // W/m.K
                    specific_heat_capacity: 2400., // J/kg.K
                    density: 17.5,                 // kg/m3... reverse engineered from paper
                },
            )
            .unwrap();

        // add material
        let mat_index = building.add_material("20mm Poly".to_string());
        building
            .set_material_properties(
                mat_index,
                MaterialProperties {
                    thickness: 20. / 1000.,
                },
            )
            .unwrap();
        building
            .set_material_substance(mat_index, poly_index)
            .unwrap();

        // Add construction
        let c_index = building.add_construction("The construction".to_string());
        building
            .add_material_to_construction(c_index, mat_index)
            .unwrap();

        // Create surface geometry
        // Geometry
        let mut the_loop = Loop3D::new();
        let l = 1. as f64;
        the_loop.push(Point3D::new(-l, -l, 0.)).unwrap();
        the_loop.push(Point3D::new(l, -l, 0.)).unwrap();
        the_loop.push(Point3D::new(l, l, 0.)).unwrap();
        the_loop.push(Point3D::new(-l, l, 0.)).unwrap();
        the_loop.close().unwrap();

        let p = Polygon3D::new(the_loop).unwrap();

        // Add surface
        let surface_index = building.add_surface("Surface".to_string());
        building
            .set_surface_construction(surface_index, c_index)
            .unwrap();
        building.set_surface_polygon(surface_index, p).unwrap();

        building
            .set_surface_front_boundary(surface_index, Boundary::Space(space_index))
            .unwrap();

        if let Ok(surf) = building.get_surface(surface_index) {
            match surf.front_boundary() {
                Boundary::Space(s) => {
                    assert_eq!(*s, space_index)
                }
                _ => assert!(false),
            }
        } else {
            assert!(false);
        }

        let n: usize = 12;
        let main_dt = 60. * 60. / n as f64;
        let model = ThermalModel::new(&mut building, &mut state, n).unwrap();

        // MAP THE STATE
        building.map_simulation_state(&mut state).unwrap();

        /* START THE TEST */
        let construction = building.get_construction(c_index).unwrap();
        assert!(!model.surfaces[0].is_massive());

        let r = r_value(&building, construction).unwrap()
            + model.surfaces[0].rs_i()
            + model.surfaces[0].rs_o();
        let u = 1. / r;
        let area = model.surfaces[0].area();

        let t_o = model.zones[0].temperature(&building, &state); // Initial T of the zone

        let t_s: f64 = 30.0; // T of surroundings

        let mut weather = SyntheticWeather::new();
        weather.dry_bulb_temperature = Box::new(ScheduleConstant::new(t_s));

        let dt = main_dt / model.dt_subdivisions() as f64;

        let mut date = Date {
            day: 1,
            hour: 0.0,
            month: 1,
        };

        // March:
        for i in 0..3000 {
            let time = (i as f64) * dt;
            date.add_seconds(time);

            let found = model.zones[0].temperature(&building, &state);
            let zone_mass = model.zones[0].mcp();

            model.march(date, &weather, &building, &mut state).unwrap();

            // Get exact solution.
            let exp = t_s + (t_o - t_s) * (-time * u * area / zone_mass).exp();
            assert!((exp - found).abs() < 0.05);
        }
    }
    /// END OF TEST_MODEL_MARCH

    #[test]
    fn test_model_march_with_window() {
        let mut state: SimulationState = SimulationState::new();
        let mut building = Building::new("The Building".to_string());

        // Add the space
        let zone_volume = 40.;
        let space_index = building.add_space("Some space".to_string());
        building.set_space_volume(space_index, zone_volume).unwrap();

        // Add substance
        let poly_index = building.add_substance("polyurethane".to_string());
        building
            .set_substance_properties(
                poly_index,
                SubstanceProperties {
                    thermal_conductivity: 0.0252,  // W/m.K
                    specific_heat_capacity: 2400., // J/kg.K
                    density: 17.5,                 // kg/m3... reverse engineered from paper
                },
            )
            .unwrap();

        // add material
        let mat_index = building.add_material("20mm Poly".to_string());
        building
            .set_material_properties(
                mat_index,
                MaterialProperties {
                    thickness: 20. / 1000.,
                },
            )
            .unwrap();
        building
            .set_material_substance(mat_index, poly_index)
            .unwrap();

        // Add construction
        let c_index = building.add_construction("The construction".to_string());
        building
            .add_material_to_construction(c_index, mat_index)
            .unwrap();

        // Create surface geometry
        // Geometry
        let mut the_loop = Loop3D::new();
        let l = 1. as f64;
        the_loop.push(Point3D::new(-l, -l, 0.)).unwrap();
        the_loop.push(Point3D::new(l, -l, 0.)).unwrap();
        the_loop.push(Point3D::new(l, l, 0.)).unwrap();
        the_loop.push(Point3D::new(-l, l, 0.)).unwrap();
        the_loop.close().unwrap();

        let mut p = Polygon3D::new(the_loop).unwrap();

        let mut the_inner_loop = Loop3D::new();
        let l = 0.5 as f64;
        the_inner_loop.push(Point3D::new(-l, -l, 0.)).unwrap();
        the_inner_loop.push(Point3D::new(l, -l, 0.)).unwrap();
        the_inner_loop.push(Point3D::new(l, l, 0.)).unwrap();
        the_inner_loop.push(Point3D::new(-l, l, 0.)).unwrap();
        the_inner_loop.close().unwrap();
        p.cut_hole(the_inner_loop.clone()).unwrap();

        // Add surface
        let surface_index = building.add_surface("Surface".to_string());
        building
            .set_surface_construction(surface_index, c_index)
            .unwrap();
        building.set_surface_polygon(surface_index, p).unwrap();
        building
            .set_surface_front_boundary(surface_index, Boundary::Space(space_index))
            .unwrap();

        // Add window.
        let window_polygon = Polygon3D::new(the_inner_loop).unwrap();
        let window_index = building.add_fenestration(
            &mut state,
            "Window One".to_string(),
            FenestrationPositions::Binary,
            FenestrationType::Window,
        );
        building
            .set_fenestration_construction(window_index, c_index)
            .unwrap();
        building
            .set_fenestration_polygon(window_index, window_polygon)
            .unwrap();
        building
            .set_fenestration_front_boundary(surface_index, Boundary::Space(space_index))
            .unwrap();

        if let Ok(surf) = building.get_surface(surface_index) {
            match surf.front_boundary() {
                Boundary::Space(s) => {
                    assert_eq!(*s, space_index)
                }
                _ => assert!(false),
            }
        } else {
            assert!(false);
        }

        // Finished building the Building

        let n: usize = 12;
        let main_dt = 60. * 60. / n as f64;
        let model = ThermalModel::new(&mut building, &mut state, n).unwrap();

        // MAP THE STATE
        building.map_simulation_state(&mut state).unwrap();

        // START TESTING.
        let construction = building.get_construction(c_index).unwrap();

        assert!(!model.surfaces[0].is_massive());

        let r = r_value(&building, construction).unwrap()
            + model.surfaces[0].rs_i()
            + model.surfaces[0].rs_o();
        let u = 1. / r;
        let area = 4.0; //the area of both the window and the wall together

        let t_o = model.zones[0].temperature(&building, &state); // Initial T of the zone

        let t_s: f64 = 30.0; // T of surroundings

        let mut weather = SyntheticWeather::new();
        weather.dry_bulb_temperature = Box::new(ScheduleConstant::new(t_s));

        let dt = main_dt / model.dt_subdivisions() as f64;

        let mut date = Date {
            day: 1,
            hour: 0.0,
            month: 1,
        };

        // March:
        for i in 0..3000 {
            let time = (i as f64) * dt;
            date.add_seconds(time);

            let found = model.zones[0].temperature(&building, &state);
            let zone_mass = model.zones[0].mcp();

            model.march(date, &weather, &building, &mut state).unwrap();

            // Get exact solution.
            let exp = t_s + (t_o - t_s) * (-time * u * area / zone_mass).exp();
            assert!((exp - found).abs() < 0.05);
        }
    }

    use building_model::heating_cooling::HeatingCoolingKind;
    use simulation_state::simulation_state_element::SimulationStateElement;
    #[test]
    fn test_model_march_with_window_and_heater() {
        let mut state: SimulationState = SimulationState::new();
        let mut building = Building::new("The Building".to_string());

        // Add the space
        let zone_volume = 40.;
        let space_index = building.add_space("Some space".to_string());
        building.set_space_volume(space_index, zone_volume).unwrap();

        let heating_power = 500.;
        building
            .add_heating_cooling_to_space(
                &mut state,
                space_index,
                HeatingCoolingKind::IdealHeaterCooler,
            )
            .unwrap();
        building
            .set_space_max_heating_power(space_index, heating_power)
            .unwrap();
        let heating_state_index = building
            .get_space(space_index)
            .unwrap()
            .get_heating_cooling()
            .unwrap()
            .state_index();
        state[heating_state_index] =
            SimulationStateElement::SpaceHeatingCoolingPowerConsumption(space_index, heating_power);

        // Add substance
        let poly_index = building.add_substance("polyurethane".to_string());
        building
            .set_substance_properties(
                poly_index,
                SubstanceProperties {
                    thermal_conductivity: 0.0252,  // W/m.K
                    specific_heat_capacity: 2400., // J/kg.K
                    density: 17.5,                 // kg/m3... reverse engineered from paper
                },
            )
            .unwrap();

        // add material
        let mat_index = building.add_material("20mm Poly".to_string());
        building
            .set_material_properties(
                mat_index,
                MaterialProperties {
                    thickness: 20. / 1000.,
                },
            )
            .unwrap();
        building
            .set_material_substance(mat_index, poly_index)
            .unwrap();

        // Add construction
        let c_index = building.add_construction("The construction".to_string());
        building
            .add_material_to_construction(c_index, mat_index)
            .unwrap();

        // Create surface geometry
        // Geometry
        let mut the_loop = Loop3D::new();
        let l = 1. as f64;
        the_loop.push(Point3D::new(-l, -l, 0.)).unwrap();
        the_loop.push(Point3D::new(l, -l, 0.)).unwrap();
        the_loop.push(Point3D::new(l, l, 0.)).unwrap();
        the_loop.push(Point3D::new(-l, l, 0.)).unwrap();
        the_loop.close().unwrap();

        let mut p = Polygon3D::new(the_loop).unwrap();

        let mut the_inner_loop = Loop3D::new();
        let l = 0.5 as f64;
        the_inner_loop.push(Point3D::new(-l, -l, 0.)).unwrap();
        the_inner_loop.push(Point3D::new(l, -l, 0.)).unwrap();
        the_inner_loop.push(Point3D::new(l, l, 0.)).unwrap();
        the_inner_loop.push(Point3D::new(-l, l, 0.)).unwrap();
        the_inner_loop.close().unwrap();
        p.cut_hole(the_inner_loop.clone()).unwrap();

        // Add surface
        let surface_index = building.add_surface("Surface".to_string());
        building
            .set_surface_construction(surface_index, c_index)
            .unwrap();
        building.set_surface_polygon(surface_index, p).unwrap();
        building
            .set_surface_front_boundary(surface_index, Boundary::Space(space_index))
            .unwrap();

        // Add window.
        let window_polygon = Polygon3D::new(the_inner_loop).unwrap();
        let window_index = building.add_fenestration(
            &mut state,
            "Window One".to_string(),
            FenestrationPositions::Binary,
            FenestrationType::Window,
        );
        building
            .set_fenestration_construction(window_index, c_index)
            .unwrap();
        building
            .set_fenestration_polygon(window_index, window_polygon)
            .unwrap();
        building
            .set_fenestration_front_boundary(surface_index, Boundary::Space(space_index))
            .unwrap();

        if let Ok(surf) = building.get_surface(surface_index) {
            match surf.front_boundary() {
                Boundary::Space(s) => {
                    assert_eq!(*s, space_index)
                }
                _ => assert!(false),
            }
        } else {
            assert!(false);
        }

        // Finished building the Building

        let n: usize = 12;
        let main_dt = 60. * 60. / n as f64;
        let model = ThermalModel::new(&mut building, &mut state, n).unwrap();

        // MAP THE STATE
        building.map_simulation_state(&mut state).unwrap();

        // START TESTING.
        let construction = building.get_construction(c_index).unwrap();
        assert!(!model.surfaces[0].is_massive());

        let r = r_value(&building, construction).unwrap()
            + model.surfaces[0].rs_i()
            + model.surfaces[0].rs_o();
        let u = 1. / r;
        let area = 4.0; //the area of both the window and the wall together

        let t_o = model.zones[0].temperature(&building, &state); // Initial T of the zone

        let t_s: f64 = 30.0; // T of surroundings

        let mut weather = SyntheticWeather::new();
        weather.dry_bulb_temperature = Box::new(ScheduleConstant::new(t_s));

        let dt = main_dt / model.dt_subdivisions() as f64;

        let mut date = Date {
            day: 1,
            hour: 0.0,
            month: 1,
        };

        // March:
        for i in 0..3000 {
            let time = (i as f64) * dt;
            date.add_seconds(time);

            let found = model.zones[0].temperature(&building, &state);
            let zone_mass = model.zones[0].mcp();

            model.march(date, &weather, &building, &mut state).unwrap();

            // Get exact solution.
            let exp = t_s
                + heating_power / (u * area)
                + (t_o - t_s - heating_power / (u * area)) * (-time * (u * area) / zone_mass).exp();
            //println!("{} vs {}", exp, found);
            assert!((exp - found).abs() < 1.0);
        }
    }
}
