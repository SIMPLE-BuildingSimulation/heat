use crate::zone::ThermalZone;
use crate::surface::*;
use building_model::building::Building;
use building_model::building_state::BuildingState;

use building_model::boundary::Boundary;
use simulation_model_trait::SimulationModel;
use weather::current_weather::CurrentWeather;


use crate::construction::discretize_construction;

pub struct ThermalModel {
    zones: Vec<ThermalZone>,
    surfaces: Vec<ThermalSurface>,
    dt_subdivisions: usize,
    dt: f64,
}

impl SimulationModel for ThermalModel{

    fn new( building: &Building, state: &mut BuildingState, n: usize )-> Self {
        
        let main_dt = 60.*60./n as f64;
        
        /* CREATE ALL ZONES, ONE PER SPACE */
        let spaces = building.get_spaces();
        let mut thermal_zones : Vec<ThermalZone> = Vec::with_capacity(spaces.len());
        for space in spaces {
            
            // Add the zone to the model... this pushes it to the sate
            // as well
            thermal_zones.push(
                ThermalZone::from_space(space, state)
            );

            
        }
        
        /* FIND MODEL TIMESTEP */
        // choose the smallest timestep in all constructions
        let max_dx = 0.04; // 4cm
        let min_dt = 60.; // 1 minute

        let mut n_subdivisions : usize = 1;
        let constructions = building.get_constructions();

        // Store the dts and n_nodes somwehere        
        let mut all_n_elements : Vec<Vec<usize>> = Vec::with_capacity(constructions.len());
        for construction in constructions{
            let(found_n_subdivisions,n_elements)=discretize_construction(building, construction, main_dt, max_dx, min_dt);            
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
        let mut thermal_surfaces : Vec<ThermalSurface> = Vec::with_capacity(surfaces.len());

        for surface in surfaces {

            let construction_index = surface.get_construction_index().unwrap();

            let i = thermal_surfaces.len();
            thermal_surfaces.push(
                ThermalSurface::new(
                    building,
                    state,
                    surface,
                    dt,
                    &all_n_elements[construction_index],
                    i
                )
            );

            // Match surface and zones        
            thermal_surfaces[i].set_front_boundary(*surface.front_boundary());
            thermal_surfaces[i].set_back_boundary(*surface.back_boundary());
            
        }



        return ThermalModel{
            zones: thermal_zones,
            surfaces: thermal_surfaces,
            dt_subdivisions: n_subdivisions,
            dt: dt
        };
    }

    fn march(&self, state: &mut BuildingState, current_weather: CurrentWeather )->Result<(),String>{

        
        let t_out = match current_weather.dry_bulb_temperature{
            Some(v)=>v,
            None => return Err(format!("Trying to march on Thermal Model, but dry bulb temperature was not provided"))
        };

        // Iterate through all the sub-subdivitions
        for _ in 0..self.dt_subdivisions{

            // update surface temperatures
            let mut heat_storage : Vec<f64> = vec![0.0; self.zones.len()];

            for i in 0..self.surfaces.len(){            
                
                // get surface
                let s = &self.surfaces[i];
                
                
                // find t_in and t_out of surface.
                let t_front = match s.front_boundary(){
                    Boundary::Space(z_index)=>self.zones[z_index].temperature(state),
                    Boundary::Ground => unimplemented!(),
                    Boundary::None => t_out
                };
                let t_back = match s.back_boundary(){
                    Boundary::Space(z_index)=>self.zones[z_index].temperature(state),
                    Boundary::Ground => unimplemented!(),
                    Boundary::None => t_out
                };            

                // Update temperatures
                let (q_front, q_back) = s.march(state, t_front,t_back);

                // Distribute heat flows.
                match s.front_boundary(){
                    Boundary::Space(z_index)=> heat_storage[z_index] += q_front*s.area() * self.dt ,
                    Boundary::Ground | Boundary::None => {}                
                };
                match s.back_boundary(){
                    Boundary::Space(z_index)=> heat_storage[z_index] += q_back*s.area() * self.dt ,
                    Boundary::Ground | Boundary::None => {}
                };

            }// end of iterating surface             


            // update air flows
                // assume constant during timestep... get a vector with Qs

            // ZONE:
            for i in 0..self.zones.len(){            
                // calculate air-flow heat transfer

                // calculate infiltration
                
                // calculate Zone heating/cooling

                // Calculate people
                
                // Calculate lighting
                        
                // update all zones temperatures
                //self.zones[i].consume_heat();
                self.zones[i].consume_heat(heat_storage[i], state);
            }// end of iterating zones

        } // End of 'in each sub-timestep-subdivision'
        
        Ok(())
    }
}

impl ThermalModel {

    
    pub fn dt_subdivisions(&self)->usize{
        self.dt_subdivisions
    }



}










/***********/
/* TESTING */
/***********/




#[cfg(test)]
mod testing{
    use super::*;
    use geometry3d::polygon3d::Polygon3D;
    use geometry3d::point3d::Point3D;
    use geometry3d::loop3d::Loop3D;


    use building_model::building::Building;
    use building_model::substance::{SubstanceProperties};
    use building_model::material::{MaterialProperties};
    

    use crate::construction::*;
    use building_model::boundary::Boundary;
    
    use weather::Weather;
    use weather::synthetic_weather::SyntheticWeather;
    use schedule::constant::ScheduleConstant;
    use calendar::date::Date;

    #[test]
    fn test_model_march(){
        
        let mut building = Building::new("the building".to_string());
        
        // Add the space
        let zone_volume = 40.;
        let space_index = building.add_space("Some space".to_string());
        building.set_space_volume(space_index,zone_volume).unwrap();

        // Add substance
        let poly_index = building.add_substance("polyurethane".to_string());
        building.set_substance_properties(poly_index, SubstanceProperties{
            thermal_conductivity: 0.0252, // W/m.K            
            specific_heat_capacity: 2400., // J/kg.K
            density: 17.5, // kg/m3... reverse engineered from paper
        }).unwrap();

        // add material
        let mat_index = building.add_material("20mm Poly".to_string());
        building.set_material_properties(mat_index, MaterialProperties{
            thickness: 20./1000.
        }).unwrap();
        building.set_material_substance(mat_index,poly_index).unwrap();

        // Add construction
        let c_index = building.add_construction("The construction".to_string());
        building.add_material_to_construction(c_index, mat_index).unwrap();

        // Create surface geometry
        // Geometry
        let mut the_loop = Loop3D::new();
        let l = 1. as f64;
        the_loop.push( Point3D::new(-l, -l, 0.)).unwrap();
        the_loop.push( Point3D::new(l, -l, 0.)).unwrap();
        the_loop.push( Point3D::new(l, l, 0.)).unwrap();
        the_loop.push( Point3D::new(-l, l, 0.)).unwrap();
        the_loop.close().unwrap();
        
        let p = Polygon3D::new(the_loop).unwrap();

        // Add surface
        let surface_index = building.add_surface("Surface".to_string());
        building.set_surface_construction(surface_index,c_index).unwrap();
        building.set_surface_polygon(surface_index, p).unwrap();
        
        building.set_surface_front_boundary(surface_index, Boundary::Space(space_index)).unwrap();

        if let Ok(surf) = building.get_surface(surface_index){
            match surf.front_boundary(){
                Boundary::Space(s)=>{
                    assert_eq!(*s,space_index)
                },
                _ => assert!(false) 
            }
        }else{
            assert!(false);
        }

        // Finished building the Building
        let mut state : BuildingState = Vec::new();

        let n : usize = 12;
        let main_dt = 60. * 60. / n as f64;
        let model = ThermalModel::new(&building, &mut state, n);
        let construction = building.get_construction(c_index).unwrap();

        // START TESTING.
        assert!(!model.surfaces[0].is_massive());

        let r = r_value(&building, construction).unwrap() + model.surfaces[0].rs_i() + model.surfaces[0].rs_o();
        let u = 1./r;
        let area = model.surfaces[0].area();
        
        let t_o = model.zones[0].temperature(&state); // Initial T of the zone
        

        let t_s : f64 = 30.0; // T of surroundings

        let mut weather = SyntheticWeather::new();
        weather.dry_bulb_temperature = Box::new(ScheduleConstant::new(t_s));

        let dt = main_dt/model.dt_subdivisions() as f64;
        

        // March:
        for i in 0..3000 {
            let time = (i as f64)*dt;
            let weather_data = weather.get_weather_data(Date{
                month: 1,
                day: 1,
                hour: 0.0
            });

            
            let found = model.zones[0].temperature(&state);
            let zone_mass = model.zones[0].mcp();
            
            model.march(&mut state,weather_data).unwrap();
            
            // Get exact solution.            
            let exp = t_s + (t_o - t_s)*(-time * u * area / zone_mass ).exp();            
            assert!((exp - found).abs() < 0.05);                                    
        }


    
    }

    
}
