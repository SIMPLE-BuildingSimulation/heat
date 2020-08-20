use crate::zone::ThermalZone;
use crate::surface::*;
use std::rc::Rc;


pub struct ThermalModel {
    zones: Vec<ThermalZone>,
    surfaces: Vec<ThermalSurface>,
}


impl ThermalModel {

    pub fn new_empty_model()-> Self {
        
        return ThermalModel{
            zones: vec![],
            surfaces: vec![],
        };

        // Create all zones

        // find minimum timestep and set it

        // Create all surfaces, using minimum tstep

        // Match zones and surfaces        
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
                Some(z_index)=>self.zones[z_index].temperature(),
                None => t_out
            };
            let t_back = match s.back_boundary(){
                Some(z_index)=>self.zones[z_index].temperature(),
                None => t_out
            };

            // Update temperatures
            let (q_front, q_back) = s.march(t_front,t_back);

            // Distribute heat flows.
            match s.front_boundary(){
                Some(z_index)=>self.zones[z_index].accumulate_heat(q_front*s.area()*dt),
                None => {}
            };
            match s.back_boundary(){
                Some(z_index)=>self.zones[z_index].accumulate_heat(q_back*s.area()*dt),
                None => {}
            };

        }                       

        // update air flows
            // assume constant during timestep... get a vector with Qs

        // ZONE:
            // calculate air-flow heat transfer

            // calculate infiltration
            
            // calculate Zone heating/cooling

            // Calculate people
            
            // Calculate lighting
        
            // update all zones temperatures
        
        for i in 0..self.zones.len(){            
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


    use building::substance::Substance;
    use building::material::Material;
    use building::construction::Construction;
    use building::surface::Surface;




    #[test]
    fn test_model_march(){

        let zone_volume = 40.;

        let mut model = ThermalModel::new_empty_model();

        
        // Add zone
        model.zones.push(ThermalZone::new("Zone 0".to_string(),zone_volume, 0));

        

        // Geometry
        let mut the_loop = Loop3D::new();
        let l = 1. as f64;
        the_loop.push( Point3D::new(-l, -l, 0.)).unwrap();
        the_loop.push( Point3D::new(l, -l, 0.)).unwrap();
        the_loop.push( Point3D::new(l, l, 0.)).unwrap();
        the_loop.push( Point3D::new(-l, l, 0.)).unwrap();
        the_loop.close().unwrap();
        
        let p = Polygon3D::new(the_loop).unwrap();


        // Materials

        let polyurethane = Substance::new(
            "polyurethane".to_string(),
            0.0252, // W/m.K            
            2400., // J/kg.K
            17.5, // kg/m3... reverse engineered from paper            
        );

        let m1 = Material::new(Rc::clone(&polyurethane),20./1000.);
        let c = Construction::new("wall 1".to_string(),vec![Rc::clone(&m1)]);
        let surface = Surface::new(p,Rc::clone(&c));
        
        // Build the thermal surface
        let main_dt = 300.0;
        let max_dx = m1.thickness()/4.0;
        let min_dt = 120.;
        let (dt,nodes)=find_dt_and_n_nodes(c.as_ref(), main_dt, 1, max_dx, min_dt);        
        

        // Add surface
        model.surfaces.push(ThermalSurface::new(Rc::clone(&surface),dt,&nodes,0)); 
        model.surfaces[0].set_front_boundary(0);
        model.zones[0].push_surface(0);

        // START TESTING.
        assert!(!model.surfaces[0].is_massive());

        let r = c.r_value()+model.surfaces[0].rs_i()+model.surfaces[0].rs_o();
        let u = 1./r;
        let area = model.surfaces[0].area();
        
        let t_o = model.zones[0].temperature(); // Initial T of the zone
        let t_s = 30.0; // T of surroundings

        


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
