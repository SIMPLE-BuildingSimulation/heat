use crate::zone::ThermalZone;
use crate::surface::*;
use building_model::building::Building;
use building_model::object_trait::ObjectTrait;
use building_model::boundary::Boundary;

pub struct ThermalModel {
    zones: Vec<ThermalZone>,
    surfaces: Vec<ThermalSurface>,
    dt_subdivisions: usize,
}


impl ThermalModel {

    pub fn new(building: &Building, main_dt: f64 )-> Self {
        
        
        /* CREATE ALL ZONES, ONE PER SPACE */
        let spaces = building.get_spaces();
        let mut thermal_zones : Vec<ThermalZone> = Vec::with_capacity(spaces.len());
        for space in spaces {
            thermal_zones.push(
                ThermalZone::new(
                    format!("ThermalZone::{}",space.name()),
                    space.volume().unwrap(),
                    thermal_zones.len(),
                )
            );
        }
        
        /* FIND MODEL TIMESTEP */
        // choose the smallest timestep in all constructions
        let max_dx = 0.04; // 4cm
        let min_dt = 60. * 2.; // 2 minutes 

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
                    surface,
                    main_dt / (n_subdivisions as f64),
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
            dt_subdivisions: n_subdivisions
        };
    }

    pub fn dt_subdivisions(&self)->usize{
        self.dt_subdivisions
    }

    pub fn march(&mut self, t_out: f64, dt: f64 /*, building::BuildingState */){
        
        // update state


        // update solar component


        // update surface temperatures
        for i in 0..self.surfaces.len(){            
            
            // get surface
            let s = &mut self.surfaces[i];
            
            
            // find t_in and t_out of surface.
            let t_front = match s.front_boundary(){
                Boundary::Space(z_index)=>self.zones[z_index].temperature(),
                Boundary::Ground => unimplemented!(),
                Boundary::None => t_out
            };
            let t_back = match s.back_boundary(){
                Boundary::Space(z_index)=>self.zones[z_index].temperature(),
                Boundary::Ground => unimplemented!(),
                Boundary::None => t_out
            };            

            // Update temperatures
            let (q_front, q_back) = s.march(t_front,t_back);

            // Distribute heat flows.
            match s.front_boundary(){
                Boundary::Space(z_index)=>self.zones[z_index].accumulate_heat(q_front*s.area()*dt),
                Boundary::Ground | Boundary::None => {}                
            };
            match s.back_boundary(){
                Boundary::Space(z_index)=>self.zones[z_index].accumulate_heat(q_back*s.area()*dt),
                Boundary::Ground | Boundary::None => {}
            };

        }                       

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
            self.zones[i].consume_heat();
        }

                
    }

}




#[cfg(test)]
mod testing{
    use super::*;
    use geometry3d::polygon3d::Polygon3D;
    use geometry3d::point3d::Point3D;
    use geometry3d::loop3d::Loop3D;


    use building_model::building::Building;
    use building_model::substance::{Substance, SubstanceProperties};
    use building_model::material::{Material,MaterialProperties};
    use building_model::construction::Construction;
    use building_model::surface::Surface;

    use crate::construction::*;
    use building_model::boundary::Boundary;


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

        let main_dt = 300.;
        let mut model = ThermalModel::new(&building, main_dt);
        let construction = building.get_construction(c_index).unwrap();

        // START TESTING.
        assert!(!model.surfaces[0].is_massive());

        let r = r_value(&building, construction).unwrap() + model.surfaces[0].rs_i() + model.surfaces[0].rs_o();
        let u = 1./r;
        let area = model.surfaces[0].area();
        
        let t_o = model.zones[0].temperature(); // Initial T of the zone
        let t_s = 30.0; // T of surroundings

        let dt = main_dt/model.dt_subdivisions() as f64;
        // March:
        for i in 0..3000 {
            let time = (i as f64)*dt;

            
            let found = model.zones[0].temperature();
            let zone_mass = model.zones[0].mcp();
            
            model.march(t_s,dt);
            
            // Get exact solution.            
            let exp = t_s + (t_o - t_s)*(-time * u * area / zone_mass ).exp();            
            assert!((exp - found).abs() < 0.05);                                    
        }


    
    }

    
}
