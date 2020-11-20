use matrix::Matrix;

use building_model::building::Building;
use building_model::object_trait::ObjectTrait;
use building_model::substance::*;
use building_model::material::*;
use building_model::surface::Surface;

use building_model::construction::Construction;
use convection::*;

use crate::construction::*;


/// This is a Surface from the point of 
/// view of our thermal solver.
pub struct ThermalSurface {
    
    /// The index of this surface within 
    /// the Thermal Model surfaces array
    index: usize,

    /// A reference to the original Surface in the 
    /// building model
    surface_index: usize,
        
    /// The interior (i.e. front side) convection coefficient
    rs_i: f64,
    
    /// The exterior (i.e. back side) convection coefficient
    rs_o: f64,

    /// The interior (i.e. front side) resistance before
    /// any layer with mass. It includes the r_si and also
    /// any light-weight material at the front of the construction. 
    /// If the first layer in the construction has mass, then this 
    /// value will be equal to r_si
    full_rsi: f64,

    /// A coefficient with the rho*Cp*dx/dt of the first 
    /// node. This is the right-hand side of the differential
    /// equation we are solving. 
    c_i: f64,

    /// A coefficient with the rho*Cp*dx/dt of the last 
    /// node. This is the right-hand side of the differential
    /// equation we are solving. 
    c_o: f64,

    /// The exterior (i.e. back side) resistance after
    /// the last layer with mass. It includes the r_so and also
    /// any light-weight material at the back of the contruction. 
    /// If the first layer in the construction has mass, then 
    /// this value will be equal to r_si
    full_rso : f64,

    /// The matrix that represents the thermal network
    k_prime : matrix::Matrix,

    /// The temperatures of the nodes
    state : matrix::Matrix,

    /// The number of nodes after discretizing 
    /// the construction
    n_nodes : usize,

    /// Has thermal mass at all?
    massive: bool,

    /// The location of the front boundary zone in the 
    /// Zones array of the Thermal Model
    front_boundary_index: Option<usize>,

    /// The location of the back boundary zone in the 
    /// Zones array of the Thermal Model    
    back_boundary_index: Option<usize>,

    /// The area of the Surface
    area: f64,
}

impl ThermalSurface {

    pub fn new(building: &Building, surface: &Surface, dt: f64, nodes: &Vec<usize>, index: usize ) -> Self {
        
        // Get the surface... or else
        //let surface = building.get_surface(surface_index).unwrap();
        
        // Check if Surface is valid... or else
        ObjectTrait::is_full(surface).unwrap();
        
        // Get construction... or else 
        
        // this (should not fail because the surface is valid)
        let construction_index = surface.get_construction_index().unwrap();
        let construction = building.get_construction(construction_index).unwrap();

        let n_nodes = calc_n_total_nodes(&nodes);
        let (rs_i,rs_o) = calc_convection_coefficients(surface);        
        
        let mut ret = ThermalSurface{
            surface_index: ObjectTrait::index(surface),            
            rs_i: rs_i,
            rs_o: rs_o,
            full_rso: 0.0, // filled when building thermal network
            full_rsi: 0.0, // filled when building thermal network
            c_o: 0.0, // filled when building thermal network
            c_i: 0.0, // filled when building thermal network
            k_prime: Matrix::new(0.0,n_nodes,n_nodes),// filled when building thermal network
            state: Matrix::new(20.0,n_nodes,1),// Set initial conditions... warmup period is required prior to simulation
            n_nodes: n_nodes,            
            massive: true, // filled after building the thermal network
            front_boundary_index: None, // filled when setting boundary
            back_boundary_index: None,
            index: index,
            area : surface.area().unwrap(), // should not fail because surface is full
        };
        
        build_thermal_network(building, &construction, dt, &nodes, rs_i, rs_o, &mut ret.k_prime, &mut ret.full_rsi, &mut ret.full_rso, &mut ret.c_i, &mut ret.c_o).unwrap();
        ret.massive = !(ret.c_o == 0.0 && ret.c_i == 0.);

        // return
        ret
                
    }
    

    pub fn area(&self)->f64{
        self.area
    }

    pub fn rs_i(&self)->f64{
        self.rs_i
    }

    pub fn rs_o(&self)->f64{
        self.rs_o
    }

    fn calc_heat_flow(&self, t_in: f64, t_out: f64)->(f64,f64){
        // Positive is going out of the layer.
        let q_in;
        let q_out;
        let t_si = self.state.get(0,0).unwrap();
        let t_so = self.state.get(self.n_nodes-1,0).unwrap();
        if self.massive {
            q_in = (t_si - t_in)/self.full_rsi;            
            q_out = (t_so - t_out)/self.full_rso;
        }else{            
            q_in = (t_si - t_in)/self.rs_i;            
            q_out = (t_so - t_out)/self.rs_o;            
        }
        return (q_in,q_out)
    }

    pub fn is_massive(&self)->bool{
        self.massive
    }
    pub fn march(&mut self, t_in: f64, t_out: f64)->(f64,f64){

        if self.massive {
            // Update state... T_i+1 = Ti + K_prime*Ti + {t_in/full_rsi/C_i ... 0,0,0... t_out/full_rso/C_o}
            //                              ^^^^^^^^^                    ^^^^^^^^^^^^
            //                        lets call this vector 'a'         These are the F components        
            let mut a = self.k_prime.from_prod_n_diag(&self.state,3).unwrap();// k_prime is tri-diagonal
            // ... 'a' should be a vector

            // Let's add the F components
            let old_value = a.get(0,0).unwrap();
            a.set(0,0,old_value + t_in/self.full_rsi/self.c_i).unwrap();
    
            let old_value = a.get(self.n_nodes-1,0).unwrap();
            a.set(self.n_nodes-1,0,old_value + t_out/self.full_rso/self.c_o).unwrap();
    
            // Let's add a to the state.
            self.state.add_to_this(&a).unwrap();
        }else{
            let q = (t_out - t_in)/self.full_rsi;

            let t_si = t_in  + q*self.rs_i;
            let t_so = t_out - q*self.rs_o;

            self.state.set(0,0,t_si).unwrap();
            self.state.set(self.n_nodes-1,0,t_so).unwrap();
        }

        return self.calc_heat_flow(t_in, t_out);
    }

    pub fn set_front_boundary(&mut self, b: usize) {
        self.front_boundary_index = Some(b);
    }

    pub fn set_back_boundary(&mut self, b: usize) {
        self.back_boundary_index = Some(b);
    }

    pub fn front_boundary(&self) -> Option<usize> {
        self.front_boundary_index
    }

    pub fn back_boundary(&self) -> Option<usize> {
        self.back_boundary_index
    }
}

fn get_first_and_last_massive(n_nodes : &Vec<usize>)->(usize,usize){

    // We need somethig to process!
    if n_nodes.len() == 0{
        panic!("Impossible to check first and last massive layers in empty discretization scheme");
    }

    // Border case
    if n_nodes.len() == 1 {
        if n_nodes[0] == 0 {
            return (0,0);
        }else{
            return (0,1);
        }
    }

    // find the first and last massive layers
    let mut first_massive = n_nodes.len();
    let mut last_massive = 0;
    let mut any_massive = false;
    for (i,nodes) in n_nodes.iter().enumerate(){
        if *nodes > 0 {
            if first_massive > i {
                first_massive = i;
            }
            if last_massive < i {
                last_massive = i;
            }
            any_massive = true;
        }
    }
    if any_massive{
        return (first_massive,last_massive+1);
    }else{
        return (0,0);
    }
}

fn calc_n_total_nodes(n_nodes : &Vec<usize>)->usize{
    

    // We need somethig to process!
    if n_nodes.len() == 0{
        panic!("Wrong discretization scheme!");
    }

    // Border case: Only one element.
    
    if n_nodes.len() == 1 {        
        if n_nodes[0] == 0 {            
            return 2;
        }else{
            return n_nodes[0]+1;
        }
    }
    

    // Otherwise, let's do this.
    let mut n : usize = 1; // the end node.
    let(first_massive,last_massive)=get_first_and_last_massive(n_nodes);

    // No massive layer at all in the construction...
    if first_massive == 0 && last_massive == 0 {
        return 2;
    }

    // Now, process from the first massive to the last.
    // interior and exterior no-mass layers do not contribute
    // nodes.
    let mut in_nomass_layer : bool = false;
    for i in first_massive..last_massive{
        let nodes = n_nodes[i];
        if nodes > 0 {
            // Material with mass.
            in_nomass_layer = false;
            n += nodes;
    
        }else{
            // material with no-mass
            if !in_nomass_layer  {
                n+=1;
                in_nomass_layer = true;
            }
        }
    }    

    return n
}

/// It is assumed to be a sandwich where a potentially massive 
/// wall is between two non-mass layers. These non-mass layers can be 
/// Film convections coefficients or ust a layer of insulation or something
/// with negligible thermal mass.
fn build_thermal_network(building: &Building, c: &Construction, dt: f64, n_nodes: &Vec<usize>, rs_i: f64, rs_o: f64, k_prime: & mut Matrix, full_rsi: &mut f64, full_rso: &mut f64, c_i: & mut f64, c_o: & mut f64 ) ->Result<(),String> {
    
    // This might fail, if there is something wrong with the model
    //let c = building.get_construction(construction_index).unwrap();
    let construction_index = c.index();

    if n_nodes.len() != c.n_layers(){        
        let err = format!("Mismatch between number of layers in construction ({}) and node scheme ({})",c.n_layers(),n_nodes.len());
        return Err(err);
    }
    
    let all_nodes = calc_n_total_nodes(n_nodes);
    let (rows,cols) = k_prime.size();
    if rows != all_nodes || cols != all_nodes {
        let err = format!("Unexpected size of given matrix - found ({},{}) and was expecting ({},{})",rows,cols,all_nodes,all_nodes);
        return Err(err);
    }
        
    // NOW, PROCESS
    ////////////////
    let (first_massive, last_massive) = get_first_and_last_massive(n_nodes);
    
    if first_massive == 0 && last_massive == 0 {
        // no massive layers at all in construction.
        // Simple case... return the following:
        let r = r_value(building, c).unwrap()+rs_i+rs_o;
        *full_rsi = r;
        *full_rso = r;
        k_prime.set(0,0,-1.0/r).unwrap();
        k_prime.set(0,1,1.0/r).unwrap();
        k_prime.set(1,0,1.0/r).unwrap();
        k_prime.set(1,1,-1.0/r).unwrap();
        *c_i=0.0;
        *c_o=0.0;
        return Ok(());
    }else{
        // Calculate inner and outer surface Resistances
        *full_rsi = rs_i;
        for i in 0..first_massive {
            let material_index = c.get_layer_index(i).unwrap();
            let material = building.get_material(material_index).unwrap();            
            let substance_index = material.get_substance_index().unwrap();
            let substance = building.get_substance(substance_index).unwrap();

            

            *full_rsi += material.thickness().unwrap() / substance.thermal_conductivity().unwrap();
        }
        *full_rso = rs_o;
        for i in last_massive..c.n_layers() {        
            let material_index = c.get_layer_index(i).unwrap();
            let material = building.get_material(material_index).unwrap();                        
            let substance_index = material.get_substance_index().unwrap();
            let substance = building.get_substance(substance_index).unwrap();
            *full_rso += material.thickness().unwrap()/substance.thermal_conductivity().unwrap();
        }
    }

    
    // Calculate the rest.
    let mut node : usize = 0;
    let mut n_layer : usize = first_massive;

    while n_layer < last_massive {

        let material_index = c.get_layer_index(n_layer).unwrap();
        let material = building.get_material(material_index).unwrap();                        

        let m = n_nodes[n_layer];                
        if m == 0 { 
            
            // nomass material
            // add up all the R of the no-mass layers that
            // are together            
            let mut r = 0.0; // if the material is no mass, then the first value 
            while n_layer < last_massive && n_nodes[n_layer] == 0 {                
                let material_index = c.get_layer_index(n_layer).unwrap();
                let material = building.get_material(material_index).unwrap();                        
                let dx = material.thickness().unwrap();
                
                let substance_index = material.get_substance_index().unwrap();
                let substance = building.get_substance(substance_index).unwrap();                
                let k =substance.thermal_conductivity().unwrap();

                r += dx / k;
                n_layer += 1;
            }      

            // update values
            let u = 1./r;
            // top left
            let old_value = k_prime.get(node,node).unwrap();
            k_prime.set(node  , node  , old_value - u ).unwrap();
            // top right
            let old_value = k_prime.get(node,node+1).unwrap();
            k_prime.set(node  , node+1, old_value + u ).unwrap();
            // bottom left
            let old_value = k_prime.get(node+1, node).unwrap();
            k_prime.set(node+1, node  , old_value + u ).unwrap();
            // bottom right
            let old_value = k_prime.get(node+1,node+1).unwrap();
            k_prime.set(node+1, node+1, old_value - u ).unwrap();

            // Move one node ahead
            node += 1;
            

        }else{
            // calc U value                                           
            let substance_index = material.get_substance_index().unwrap();
            let substance = building.get_substance(substance_index).unwrap();

            let k = substance.thermal_conductivity().unwrap();                
            let dx = material.thickness().unwrap()/(m as f64);            
            let u : f64 = k/dx;

            for _i in 0..m {                                                                            
                // top left
                let old_value = k_prime.get(node,node).unwrap();
                k_prime.set(node  , node  , old_value - u ).unwrap();
                // top right
                let old_value = k_prime.get(node,node+1).unwrap();
                k_prime.set(node  , node+1, old_value + u ).unwrap();
                // bottom left
                let old_value = k_prime.get(node+1, node).unwrap();
                k_prime.set(node+1, node  , old_value + u ).unwrap();
                // bottom right
                let old_value = k_prime.get(node+1,node+1).unwrap();
                k_prime.set(node+1, node+1, old_value - u ).unwrap();
                
                // advance node.                
                node += 1;
            }
            n_layer += 1;
        }

    }       
    

    // ADD RSI AND RSO 
    // add r_si to first node            
    if *full_rsi > 0.0 {
        let old_value = k_prime.get(0,0).unwrap();        
        k_prime.set(0,0, old_value - 1.0/ (*full_rsi) ).unwrap();        
    }

    if *full_rso > 0.0{        
        let old_value = k_prime.get(all_nodes-1,all_nodes-1).unwrap();
        k_prime.set(all_nodes-1,all_nodes-1, old_value - 1.0/(*full_rso) ).unwrap();        
    }
    
    

    // CALCULATE MASSES
    let mut left_side : Vec<f64> = vec![0.0;all_nodes];
    node = 0;

    for n_layer in 0..c.n_layers() {        

        let layer_index = c.get_layer_index(n_layer).unwrap();
        let material = building.get_material(layer_index).unwrap();                        
        
        let m = n_nodes[n_layer];  

        if m != 0 {
            for _i in 0..m {                   
                // Calc mass                
                let substance_index = material.get_substance_index().unwrap();
                let substance = building.get_substance(substance_index).unwrap();
                
                let rho = substance.density().unwrap();
                let cp = substance.specific_heat_capacity().unwrap();
                let dx = material.thickness().unwrap()/(m as f64);
                let m = rho * cp * dx / dt;                     
                
                left_side[node]+=m/2.0;
                left_side[node+1]+=m/2.0;
               
                // advance node.
                node += 1;
            }
        }else{
            // else, they are zero already...
            // still, we need to move forward one
            // node.            
            if n_layer > 0 && n_nodes[n_layer-1] > 0 { 
                // Only do it if the previous
                // layer was massive
                node+=1;
            }            
        }
    }         
    
    
    // DIVIDE ONE BY THE OTHER 
    for i in 0..all_nodes{
        let mass = left_side[i];
        
        if mass != 0.0 {                        
            // Multiply the whole column by this.
            for j in 0..all_nodes{
                let old_k_value = k_prime.get(i,j).unwrap();
                k_prime.set(i,j, old_k_value / mass).unwrap();                   
            }
        }// Else, masses are already Zero.
                
    }    

    // SET C_I AND C_O
    *c_i = left_side[0];
    *c_o = left_side[all_nodes-1];
    
    return Ok(())
}

/// Given a Maximum thickness (dx) and a minimum timestep (dt), this function
/// will find a good combination of dt and number of nodes in each
/// layer of the construction. N is the iteration number; it should always start
/// at 1.
pub fn find_dt_and_n_nodes(building: &Building, c: &Construction, main_dt: f64, n: usize, max_dx: f64, min_dt: f64) -> (f64,Vec<usize>){

    let dt = main_dt/(n as f64);
    let safety = 1.5;

    // Choose a dx so that dt allows convergence.
    // stability is assured by (alpha * dt / dx^2 <= 1/2 )
    // meaning, we need to satisfy dx >= sqrt(2 * alpha * dt)
    
    // So, for each layer
    let mut n_nodes : Vec<usize> = vec![];

    for n_layer in 0..c.n_layers(){
        let material_index = c.get_layer_index(n_layer).unwrap();
        let material = building.get_material(material_index).unwrap();                        
        let substance_index = material.get_substance_index().unwrap();
        let substance = building.get_substance(substance_index).unwrap();

        // Calculate the optimum_dx
        let thickness = material.thickness().unwrap();
        let alpha = substance.thermal_diffusivity().unwrap();
        let optimum_dx = (safety*2.0*alpha*dt).sqrt();
        let m = (thickness/optimum_dx).floor();            
        let mut dx = thickness/m;
        
        
        // dx cannot be greater than thickness
        if m == 0. {
            dx = thickness;
        }        

        // If dx is larger than the max allowed d_x, try to change timestep        
        if dx > max_dx {
            // check if there is room for reducing dt...
            let next_dt = main_dt/((n+1) as f64);
            if next_dt > min_dt {
                // If there is room for that, do it.
                return find_dt_and_n_nodes(building,c,main_dt,n+1,max_dx,min_dt);                
            }else{
                // otherwise, mark this layer as no-mass

                // Zero nodes indicates no-mass, although
                // a node will be added in the boundary
                // with previous layer, if prev. layer has 
                // mass.
                n_nodes.push(0);                
            }
        } else {
            // "normal" case.            
            n_nodes.push(m as usize)
        }
    }
    
    return (dt,n_nodes)
}



/***********/
/* TESTING */
/***********/



#[cfg(test)]
mod testing{
    use super::*;
    use building_model::material::Material;
    use building_model::substance::Substance;
    use geometry3d::polygon3d::Polygon3D;
    use geometry3d::point3d::Point3D;
    use geometry3d::loop3d::Loop3D;

    fn add_polyurethane(building: &mut Building)->usize{
        let poly_index = building.add_substance("polyurethane".to_string());
        building.set_substance_properties(poly_index, SubstanceProperties{
            density: 17.5, // kg/m3... reverse engineered from paper            
            specific_heat_capacity: 2400., // J/kg.K
            thermal_conductivity: 0.0252, // W/m.K            
        }).unwrap();
        {
            let poly = building.get_substance(poly_index).unwrap();
            assert_eq!(poly.thermal_diffusivity().unwrap(),0.6E-6);
        }      

        poly_index
    }

    fn add_brickwork(building: &mut Building)->usize{
        let brickwork_index = building.add_substance("brickwork".to_string());

        building.set_substance_properties(brickwork_index, SubstanceProperties{
            density: 1700., // kg/m3... reverse engineered from paper            
            specific_heat_capacity: 800., // J/kg.K
            thermal_conductivity: 0.816, // W/m.K            
        }).unwrap();
        {
            let brickwork = building.get_substance(brickwork_index).unwrap();
            assert_eq!(brickwork.thermal_diffusivity().unwrap(),0.6E-6);
        }
        brickwork_index
    }

    fn add_layer(building: &mut Building, substance_index: usize, thickness: f64)->usize{
        let m0_index = building.add_material("123123".to_string());
        building.set_material_substance(m0_index, substance_index).unwrap();
        building.set_material_properties(m0_index, MaterialProperties{
            thickness: thickness
        }).unwrap();

        m0_index
    }

    fn get_wall_1()->Building{
        let mut building = Building::new("The Building".to_string());

        /* SUBSTANCES */

        // Add polyurethane Substance
        let poly_index = add_polyurethane(&mut building);

        /* MATERIAL */
        let m0_thickness = 200.0/1000. as f64;
        let m0_index = add_layer(&mut building,poly_index,m0_thickness);        
        
        /* WALL 1 */
        let c0_index = building.add_construction("Wall 1".to_string());
        building.add_material_to_construction(c0_index, m0_index).unwrap();        
        
        return building;
    }


    fn get_wall_2()->Building{
        let mut building = Building::new("The Building".to_string());

        /* SUBSTANCES */
        // Add polyurethane Substance
        let poly_index = add_polyurethane(&mut building);

        // Add brickwork
        let brickwork_index = add_brickwork(&mut building);

        /* MATERIALS */
        let m0_thickness = 200.0/1000. as f64;
        let m0_index = add_layer(&mut building, poly_index, m0_thickness);    

        let m1_thickness = 110.0/1000. as f64;
        let m1_index = add_layer(&mut building, brickwork_index, m1_thickness);        

        /* WALL 2 */

        let c0_index = building.add_construction("Wall 2".to_string());
        building.add_material_to_construction(c0_index, m0_index).unwrap();        
        building.add_material_to_construction(c0_index, m1_index).unwrap();        

        
        return building;
    }

    fn get_wall_3()->Building{
        let mut building = Building::new("The Building".to_string());

        /* SUBSTANCES */
        // Add polyurethane Substance
        let poly_index = add_polyurethane(&mut building);        

        // Add brickwork
        let brickwork_index = add_brickwork(&mut building);
        

        /* MATERIALS */
        let poly_mat_thickness = 20.0/1000. as f64;
        let poly_mat_index = add_layer(&mut building, poly_index, poly_mat_thickness);
        
        
        let brickwork_mat_thickness = 220.0/1000. as f64;
        let brickwork_mat_index = add_layer(&mut building, brickwork_index, brickwork_mat_thickness);
        
        /* WALL 3 */

        let c0_index = building.add_construction("Wall 3".to_string());
        building.add_material_to_construction(c0_index, poly_mat_index).unwrap();        
        building.add_material_to_construction(c0_index, brickwork_mat_index).unwrap();        
        building.add_material_to_construction(c0_index, poly_mat_index).unwrap();        

        
        return building;
    }

    #[test]
    fn test_new(){

        let mut building = Building::new("A building".to_string());

        /* SUBSTANCES */
        // Add polyurethane Substance
        let poly_index = add_polyurethane(&mut building);
        

        let m0_thickness = 200.0/1000. as f64;
        let m0_index = add_layer(&mut building, poly_index, m0_thickness);
        

        /* CONSTRUCTION */
        let c0_index = building.add_construction("Wall".to_string());
        building.add_material_to_construction(c0_index, m0_index).unwrap();        


        

        /* GEOMETRY */
        let mut the_loop = Loop3D::new();
        let l = 1. as f64;
        the_loop.push( Point3D::new(-l, -l, 0.)).unwrap();
        the_loop.push( Point3D::new(l, -l, 0.)).unwrap();
        the_loop.push( Point3D::new(l, l, 0.)).unwrap();
        the_loop.push( Point3D::new(-l, l, 0.)).unwrap();
        the_loop.close().unwrap();
        let p = Polygon3D::new(the_loop).unwrap();
        

        /* SURFACE */
        let s_index = building.add_surface("Surface".to_string());
        building.set_surface_polygon(s_index, p).unwrap();
        building.set_surface_construction(s_index, c0_index).unwrap();
        
        
        let m0 = building.get_material(m0_index).unwrap();
        let c0 = building.get_construction(c0_index).unwrap();
        let surface = building.get_surface(s_index).unwrap();

        // This should result in a dt of 300, with 8 layers of 0.25m thick.
        let main_dt = 300.0;
        let max_dx = m0.thickness().unwrap()/4.0;
        let min_dt = 1.;
        let (dt,nodes)=find_dt_and_n_nodes(&building,&c0, main_dt, 1, max_dx, min_dt);
        let ts = ThermalSurface::new(&building,&surface,dt,&nodes,0);

        let (rs_i,rs_o)=calc_convection_coefficients(&surface);
        assert!(ts.massive);
        assert_eq!(ts.n_nodes,9);
        assert_eq!(ts.rs_i,rs_i);
        assert_eq!(ts.rs_o,rs_o);
        assert_eq!(ts.area,4.0);



    }

    #[test]
    fn test_get_dt_and_n_nodes_wall_1(){
        let building = get_wall_1();

        let c0_index = 0;
        let c0 = building.get_construction(c0_index).unwrap();
        
        let m0_index = 0;
        let m0 = building.get_material(m0_index).unwrap();
        let m0_thickness = m0.thickness().unwrap();


        // This should result in a dt of 300, with 8 layers of 0.25m thick.
        let main_dt = 300.;
        let (dt,n_nodes) = find_dt_and_n_nodes(&building, &c0, main_dt, 1, m0_thickness/4., 1.);        

        assert_eq!(n_nodes.len(),c0.n_layers());
        assert_eq!(dt,300.);
        assert_eq!(n_nodes[0],8);
        

        // This should result in a dt of 300, with 8 layers of 0.25m thick.
        let main_dt = 300.;
        let (dt,n_nodes) = find_dt_and_n_nodes(&building, &c0, main_dt, 1, m0_thickness/8., 1.);                
        assert_eq!(n_nodes.len(),c0.n_layers());
        assert_eq!(dt,300.);
        assert_eq!(n_nodes[0],8);
        

        // This should result in a dt of 300/2=150, with 12 layers of 0.0166667m thick.
        let main_dt = 300.;
        let (dt,n_nodes) = find_dt_and_n_nodes(&building, &c0, main_dt, 1, m0_thickness/9., 1.);                                                      
        assert_eq!(n_nodes.len(),c0.n_layers());
        assert_eq!(dt,150.);
        assert_eq!(n_nodes[0],12);
        

        // This should result in a dt of 300/4=75, with 17 layers of 0.011...m thick.
        let main_dt = 300.;
        let (dt,n_nodes) = find_dt_and_n_nodes(&building, &c0, main_dt, 1, m0_thickness/15., 1.);                                                              
        assert_eq!(n_nodes.len(),c0.n_layers());
        assert_eq!(dt,75.);
        assert_eq!(n_nodes[0],17);
        

        // We imposed a min_dt of 80, meaning that this should result in a no-mass layer
        let main_dt = 300.;
        let (dt,n_nodes) = find_dt_and_n_nodes(&building, &c0, main_dt, 1, m0_thickness/15., 80.);                                                                                         
        assert_eq!(n_nodes.len(),c0.n_layers());
        assert_eq!(dt,100.);
        assert_eq!(n_nodes[0],0);
        
    }


    #[test]
    fn test_get_dt_and_n_nodes_wall_2(){
        let building = get_wall_2();

        let c0_index = 0;
        let c0 = building.get_construction(c0_index).unwrap();
        
        //let m0_index = 0;
        //let m0 = building.get_material(m0_index).unwrap();
        //let m0_thickness = m0.thickness().unwrap();

        //let m1_index = 0;
        //let m1 = building.get_material(m1_index).unwrap();
        //let m1_thickness = m1.thickness().unwrap();


        // WALL 2

        // This should result in a dt of 300, with 8 layers and 4 layers.
        let main_dt = 300.;
        let (dt,n_nodes) = find_dt_and_n_nodes(&building, &c0, main_dt, 1, 0.03, 1.);                           
        assert_eq!(n_nodes.len(),2);
        assert_eq!(dt,300.);
        assert_eq!(n_nodes[0],8);
        assert_eq!(n_nodes[1],4);
        
        

        // This should result in a dt of 300/2=150, with 12 and 6 layers.
        let main_dt = 300.;
        let (dt,n_nodes) = find_dt_and_n_nodes(&building, &c0, main_dt, 1, 0.02, 1.);
        assert_eq!(n_nodes.len(),c0.n_layers());
        assert_eq!(dt,150.);
        assert_eq!(n_nodes[0],12);
        assert_eq!(n_nodes[1],6);
        

        // This should result in a dt of 300/3=100, with 14 and 8
        let main_dt = 300.;
        let (dt,n_nodes) = find_dt_and_n_nodes(&building, &c0, main_dt, 1, 0.015, 1.);                            
        assert_eq!(n_nodes.len(),c0.n_layers());
        assert_eq!(dt,100.);
        assert_eq!(n_nodes[0],14);
        assert_eq!(n_nodes[1],8);
        

        // We imposed a min_dt of 100, meaning that this should result in a no-mass layer
        //let main_dt = 300.;
        let (dt,n_nodes) = find_dt_and_n_nodes(&building, &c0, main_dt, 1, 0.015, 110.);        
        assert_eq!(n_nodes.len(),c0.n_layers());
        assert_eq!(dt,150.);
        assert_eq!(n_nodes[0],0);
        assert_eq!(n_nodes[1],0);
    }

    #[test]
    fn test_get_dt_and_n_nodes_wall_3(){
        
        let building = get_wall_3();
        
        let c0_index = 0;
        let c0 = building.get_construction(c0_index).unwrap();
        
        //let m0_index = 0;
        //let m0 = building.get_material(m0_index).unwrap();
        //let m0_thickness = m0.thickness().unwrap();

        //let m1_index = 0;
        //let m1 = building.get_material(m1_index).unwrap();
        //let m1_thickness = m1.thickness().unwrap();

        //let m2_index = 0;
        //let m2 = building.get_material(m2_index).unwrap();
        //let m2_thickness = m2.thickness().unwrap();

        // This should result in a dt of 150, with [1,13,1] layers
        let main_dt = 300.;
        let (dt,n_nodes) = find_dt_and_n_nodes(&building, &c0, main_dt, 1, 0.03, 1.);                            
        {
            print!("N-NODES -> [");
            for n in &n_nodes{
                print!("{},",n);
            }
            println!("]");
        }

        assert_eq!(n_nodes.len(),c0.n_layers());
        assert_eq!(dt,300.);
        assert_eq!(n_nodes[0],0);        
        assert_eq!(n_nodes[1],9);
        assert_eq!(n_nodes[2],0);
        
        

        // This should result in a dt of 300/6=50, with [2,23,2] layers
        let main_dt = 300.;
        let (dt,n_nodes) = find_dt_and_n_nodes(&building, &c0, main_dt, 1, 0.015, 1.);                                                        
        assert_eq!(n_nodes.len(),c0.n_layers());
        assert_eq!(dt,50.);
        assert_eq!(n_nodes[0],2);
        assert_eq!(n_nodes[1],23);
        assert_eq!(n_nodes[2],2);
        

        // Limit min_time to 65... This should result in a dt of 300/4=75, with [0, 18, 0] layers
        let main_dt = 300.;
        let (dt,n_nodes) = find_dt_and_n_nodes(&building, &c0, main_dt, 1, 0.015, 65.);                                                                                    
        assert_eq!(n_nodes.len(),c0.n_layers());
        assert_eq!(dt,75.);
        assert_eq!(n_nodes[0],0);
        assert_eq!(n_nodes[1],18);
        assert_eq!(n_nodes[2],0);
                
    }


    #[test]
    fn test_find_n_total_nodes(){
        let n_nodes = vec![2,2];
        let (first,last) = get_first_and_last_massive(&n_nodes);
        assert_eq!(first,0);
        assert_eq!(last,2);
        assert_eq!(5,calc_n_total_nodes(&n_nodes));
        

        let n_nodes = vec![1,0,0,2];
        let (first,last) = get_first_and_last_massive(&n_nodes);
        assert_eq!(first,0);
        assert_eq!(last,4);
        assert_eq!(5,calc_n_total_nodes(&n_nodes));

        let n_nodes = vec![1,0,2];
        let (first,last) = get_first_and_last_massive(&n_nodes);
        assert_eq!(first,0);
        assert_eq!(last,3);
        assert_eq!(5,calc_n_total_nodes(&n_nodes));

        let n_nodes = vec![8];
        let (first,last) = get_first_and_last_massive(&n_nodes);
        assert_eq!(first,0);
        assert_eq!(last,1);
        assert_eq!(9,calc_n_total_nodes(&n_nodes));

        let n_nodes = vec![0];
        let (first,last) = get_first_and_last_massive(&n_nodes);
        assert_eq!(first,0);
        assert_eq!(last,0);
        assert_eq!(2,calc_n_total_nodes(&n_nodes));        

        let n_nodes = vec![1,0];
        let (first,last) = get_first_and_last_massive(&n_nodes);
        assert_eq!(first,0);
        assert_eq!(last,1);
        assert_eq!(2,calc_n_total_nodes(&n_nodes));        


        let n_nodes = vec![0,1];
        let (first,last) = get_first_and_last_massive(&n_nodes);
        assert_eq!(first,1);
        assert_eq!(last,2);
        assert_eq!(2,calc_n_total_nodes(&n_nodes));        

        let n_nodes = vec![0,1,0];
        let (first,last) = get_first_and_last_massive(&n_nodes);
        assert_eq!(first,1);
        assert_eq!(last,2);
        assert_eq!(2,calc_n_total_nodes(&n_nodes));        

        let n_nodes = vec![0,0,0];
        let (first,last) = get_first_and_last_massive(&n_nodes);
        assert_eq!(first,0);
        assert_eq!(last,0);
        assert_eq!(2,calc_n_total_nodes(&n_nodes));        

        let n_nodes = vec![1,0,1];
        let (first,last) = get_first_and_last_massive(&n_nodes);
        assert_eq!(first,0);
        assert_eq!(last,3);
        assert_eq!(4,calc_n_total_nodes(&n_nodes));        
        
    }

    #[test]
    fn test_wall_brickwork(){
        
        let mut building = Building::new("The building".to_string());
        
        // Add brickwork
        let brickwork_index = add_brickwork(&mut building);
     

        /* WALL 1: Single layer, with mass. */
        let brickwork_200_thickness = 200.0/1000. as f64;
        let brickwork_200_index = add_layer(&mut building, brickwork_index, brickwork_200_thickness);
        
        
        let wall_1_index = building.add_construction("Wall 1".to_string());
        building.add_material_to_construction(wall_1_index, brickwork_200_index).unwrap();        
        
        let wall_1 = building.get_construction(wall_1_index).unwrap();
    
        let dt = 156.0;
        let n_nodes : Vec<usize> = vec![4];
        let all_nodes = calc_n_total_nodes(&n_nodes);

        let mut k_prime = Matrix::new(0.0,all_nodes,all_nodes);
        let mut full_rsi=0.0;
        let mut full_rso=0.0;
        let mut c_i=0.0;
        let mut c_o=0.0;    
        build_thermal_network(&building, &wall_1, dt, &n_nodes, 0., 0., &mut k_prime, &mut full_rsi, &mut full_rso, &mut c_i, &mut c_o).unwrap();            
        
            
        let substance = building.get_substance(brickwork_index).unwrap();
        let rho = substance.density().unwrap();
        let cp = substance.specific_heat_capacity().unwrap();
        let k = substance.thermal_conductivity().unwrap();
        
        let n_layer = 0;                        
        let m = n_nodes[n_layer];
        let dx = brickwork_200_thickness/(m as f64);
        let substance = building.get_substance(brickwork_index).unwrap();
        let alpha = substance.thermal_diffusivity().unwrap();
        let mass = rho * cp * dx / dt;

        // check coefficients
        assert_eq!(full_rsi,0.0);// 2.0*dt/(rho*cp*dx));
        assert_eq!(full_rso,0.0);// 2.0*dt/(rho*cp*dx));
        assert_eq!(c_i, rho*cp*dx/(2.0*dt));
        assert_eq!(c_o, rho*cp*dx/(2.0*dt));
        

        
        
        
        
        // Check first node
        let node = 0;
        
        let found = k_prime.get(node,node).unwrap();
        let exp = -2.0*alpha*dt/dx/dx;
        
        //assert!( (found-exp).abs() < 1E-10 );        
        assert_eq!(exp,found);

        let found = k_prime.get(node,node + 1).unwrap();
        let exp = 2.0*k/dx/mass;
        assert!( (found-exp).abs() < 1E-10 );

        // check middle nodes
        for node in 1..all_nodes-1{                        
            let found = k_prime.get(node,node).unwrap();
            let exp = -2.0*alpha*dt/dx/dx;
            //assert!( (found-exp).abs() < 1E-10 );
            assert_eq!(found,exp);

            let found = k_prime.get(node,node + 1).unwrap();
            let exp = k/dx/mass;
            assert!( (found-exp).abs() < 1E-10 );

            let found = k_prime.get(node,node - 1).unwrap();
            let exp = k/dx/mass;
            assert!( (found-exp).abs() < 1E-10 );
        }
        // check end node
        let node = all_nodes-1;                

        let found = k_prime.get(node,node).unwrap();
        let exp = -2.0*alpha*dt/dx/dx;
        assert!( (found-exp).abs() < 1E-10 );        

        let found = k_prime.get(node,node - 1).unwrap();
        let exp = 2.0*k/dx/mass;
        assert!( (found-exp).abs() < 1E-10 );

    }


    #[test]
    fn test_wall_2(){

        let building = get_wall_2();
        let c0_index = 0;
        let c0 = building.get_construction(c0_index).unwrap();
        
        let m0_index = 0;
        let m0 = building.get_material(m0_index).unwrap();
        //let m0_thickness = m0.thickness().unwrap();
        
        let m0_substance_index = m0.get_substance_index().unwrap();
        let m0_substance = building.get_substance(m0_substance_index).unwrap();
        

        let m1_index = 1;
        let m1 = building.get_material(m1_index).unwrap();
        //let m1_thickness = m1.thickness().unwrap();

        let m1_substance_index = m1.get_substance_index().unwrap();
        let m1_substance = building.get_substance(m1_substance_index).unwrap();
        

        let dt = 156.0;
        let n_nodes : Vec<usize> = vec![3,3];
        let all_nodes = calc_n_total_nodes(&n_nodes);

        let cp0=m0_substance.specific_heat_capacity().unwrap();
        let rho0=m0_substance.density().unwrap();
        let dx0=m0.thickness().unwrap()/3.;
        let mass0 = rho0*cp0*dx0/dt;

        let cp1=m1_substance.specific_heat_capacity().unwrap();
        let rho1=m1_substance.density().unwrap();
        let dx1=m1.thickness().unwrap()/3.;
        let mass1 = rho1*cp1*dx1/dt;

        let mut k_prime = Matrix::new(0.0,all_nodes,all_nodes);
        let mut full_rsi=0.0;
        let mut full_rso=0.0;
        let mut c_i=0.0;
        let mut c_o=0.0;        
        build_thermal_network(&building, &c0, dt, &n_nodes, 0.,0., &mut k_prime, &mut full_rsi, &mut full_rso, &mut c_i, &mut c_o).unwrap();

        // check coefficients
        assert_eq!(full_rsi,0.0);// 2.0*dt/(rho*cp*dx));
        assert_eq!(full_rso,0.0);// 2.0*dt/(rho*cp*dx));    
        assert_eq!(c_i, rho0*cp0*dx0/(2.0*dt));
        assert_eq!(c_o, rho1*cp1*dx1/(2.0*dt));


        // FIRST LAYER
        //////////////
        let n_layer = 0;
        //let rho = m1.substance.density;
        //let cp = m1.substance.heat_capacity;
        let m = n_nodes[n_layer];
        let dx = m0.thickness().unwrap()/(m as f64);
        let k = m0_substance.thermal_conductivity().unwrap();
        let alpha = m0_substance.thermal_diffusivity().unwrap();

        // Check first node in layer
        let node = 0;        
        
        let found = k_prime.get(node,node).unwrap();
        let exp = -2.0*alpha*dt/dx/dx;
        assert!( (found-exp).abs() < 1E-10 );        

        let found = k_prime.get(node,node + 1).unwrap();
        let exp = 2.0*k/dx/mass0;        
        assert!( (found-exp).abs() < 1E-10 );

        // check middle nodes
        for node in 1..3{                        

            let found = k_prime.get(node,node).unwrap();
            let exp = -2.0*alpha*dt/dx/dx;
            assert!( (found-exp).abs() < 1E-10 );

            let found = k_prime.get(node,node + 1).unwrap();
            let exp = k/dx/mass0;
            assert!( (found-exp).abs() < 1E-10 );

            let found = k_prime.get(node,node - 1).unwrap();
            let exp = k/dx/mass0;
            assert!( (found-exp).abs() < 1E-10 );
        }
        // check middle node (e.g. node 3), which 
        // is also the beginning of Layer 2
        let node = 3;
        let mass_0 = m0_substance.density().unwrap() * m0_substance.specific_heat_capacity().unwrap()*m0.thickness().unwrap()/(6.0*dt);        
        let mass_1 = m1_substance.density().unwrap() * m1_substance.specific_heat_capacity().unwrap()*m1.thickness().unwrap()/(6.0*dt);        
        
        let found = k_prime.get(node,node).unwrap();
        let u0 = 3.0*m0_substance.thermal_conductivity().unwrap()/m0.thickness().unwrap();
        let u1 = 3.0*m1_substance.thermal_conductivity().unwrap()/m1.thickness().unwrap();
        let exp = -(u0+u1)/(mass_0+mass_1);        
        assert_eq!(found,exp);
        
        let found = k_prime.get(node,node - 1).unwrap();        
        assert_eq!(found, u0/(mass0/2. + mass1/2.0));

        let found = k_prime.get(node,node + 1).unwrap();        
        assert_eq!(found, u1/(mass0/2.0+mass1/2.0));

        // SECOND LAYER
        //////////////
        let n_layer = 1;
        //let rho = m2.substance.density;
        //let cp = m2.substance.heat_capacity;
        let m = n_nodes[n_layer];
        let dx = m1.thickness().unwrap()/(m as f64);
        let k = m1_substance.thermal_conductivity().unwrap();
        let alpha = m1_substance.thermal_diffusivity().unwrap();

        
        // check middle nodes in layer 2
        for node in 4..6{                        

            let found = k_prime.get(node,node).unwrap();
            let exp = -2.0*alpha*dt/dx/dx;
            assert!( (found-exp).abs() < 1E-10 );

            let found = k_prime.get(node,node + 1).unwrap();
            let exp = k/dx/mass1;
            assert!( (found-exp).abs() < 1E-10 );

            let found = k_prime.get(node,node - 1).unwrap();
            let exp = k/dx/mass1;
            assert!( (found-exp).abs() < 1E-10 );
        }
        // check end node 
        let node = 6;        
        
        
        let found = k_prime.get(node,node).unwrap();        
        let u2 = 3.0*m1_substance.thermal_conductivity().unwrap()/m1.thickness().unwrap();
        let exp = -2.0*u2/(mass1);        
        assert_eq!(found,exp);
        
        let found = k_prime.get(node,node - 1).unwrap();        
        assert_eq!(found, 2.0*u1/mass1);
    }

    #[test]
    fn test_wall_1(){

        let building = get_wall_1();
       
        let c0_index = 0;
        let c0 = building.get_construction(c0_index).unwrap();

        let m0_index = 0;
        let m0 = building.get_material(m0_index).unwrap();
        let m0_substance_index = m0.get_substance_index().unwrap();
        let m0_substance = building.get_substance(m0_substance_index).unwrap();
        
        

        let dt = 156.0;
        let n_nodes : Vec<usize> = vec![0];
        let all_nodes = calc_n_total_nodes(&n_nodes);
        let mut k_prime = Matrix::new(0.0,all_nodes,all_nodes);
        
        let mut full_rsi=0.0;
        let mut full_rso=0.0;
        let mut c_i=0.0;
        let mut c_o=0.0;        
        build_thermal_network(&building, &c0, dt, &n_nodes, 0.,0., &mut k_prime, &mut full_rsi, &mut full_rso, &mut c_i, &mut c_o).unwrap();
        
        // check coefficients
        assert_eq!(full_rsi, m0.thickness().unwrap()/m0_substance.thermal_conductivity().unwrap());// 2.0*dt/(rho*cp*dx));
        assert_eq!(full_rso, m0.thickness().unwrap()/m0_substance.thermal_conductivity().unwrap());// 2.0*dt/(rho*cp*dx));
        assert_eq!(c_i, 0.0);
        assert_eq!(c_o, 0.0);     
        
        for row in 0..all_nodes{
            for col in 0..all_nodes{
                let exp = m0_substance.thermal_conductivity().unwrap()/m0.thickness().unwrap();            
                let found = k_prime.get(row,col).unwrap();
                if row == col{
                    assert_eq!(-exp,found);
                }else{
                    assert_eq!(exp,found);
                }
            }            
        }        
    }

    #[test]
    fn test_double_wall_1(){

        let mut building = get_wall_1();

        let c0_index = 0;
        let m0_index = 0;

        // add second layer
        building.add_material_to_construction(c0_index, m0_index).unwrap();

        let c0 = building.get_construction(c0_index).unwrap();
        let m0 = building.get_material(m0_index).unwrap();

        let m0_substance_index = m0.get_substance_index().unwrap();
        let m0_substance = building.get_substance(m0_substance_index).unwrap();
        
        
        let dt = 156.0;
        let n_nodes : Vec<usize> = vec![0,0];
        let all_nodes = calc_n_total_nodes(&n_nodes);
        let r = 2.0*m0.thickness().unwrap()/m0_substance.thermal_conductivity().unwrap();
        let mut k_prime = Matrix::new(0.0,all_nodes,all_nodes);
        let mut full_rsi=0.0;
        let mut full_rso=0.0;
        let mut c_i=0.0;
        let mut c_o=0.0;                
        build_thermal_network(&building, &c0, dt, &n_nodes, 0.,0., &mut k_prime, &mut full_rsi, &mut full_rso, &mut c_i, &mut c_o).unwrap();

        // check coefficients
        assert_eq!(full_rsi, r);
        assert_eq!(full_rso, r);
        assert_eq!(c_i, 0.);
        assert_eq!(c_o, 0.);
                

        for row in 0..all_nodes{
            for col in 0..all_nodes{
                let exp = 1.0/r;
                let found = k_prime.get(row,col).unwrap();
                if row == col{
                    assert_eq!(-exp,found);
                }else{
                    assert_eq!(exp,found);
                }
            }            
        }        
    }


    #[test]
    fn test_wall_5(){

        let mut building = Building::new("A building".to_string());

        /* SUBSTANCES */        
        let poly_index = add_polyurethane(&mut building);        
        let brickwork_index = add_brickwork(&mut building);

        /* MATERIALS */
        let m0_thickness = 200.0/1000. as f64;
        let m0_index = add_layer(&mut building, poly_index, m0_thickness);
                

        let m1_thickness = 200.0/1000. as f64;
        let m1_index = add_layer(&mut building, brickwork_index, m1_thickness);
        
        /* WALL */

        let c0_index = building.add_construction("Wall 2".to_string());
        building.add_material_to_construction(c0_index, m0_index).unwrap();        
        building.add_material_to_construction(c0_index, m1_index).unwrap();        

        
        let m0 = building.get_material(m0_index).unwrap();
        let m0_substance_index = m0.get_substance_index().unwrap();
        let m0_substance = building.get_substance(m0_substance_index).unwrap();

        let m1 = building.get_material(m1_index).unwrap();
        let m1_substance_index = m1.get_substance_index().unwrap();
        let m1_substance = building.get_substance(m1_substance_index).unwrap();
        
        let c0 = building.get_construction(c0_index).unwrap();

        let dt = 156.0;
        let n_nodes : Vec<usize> = vec![1,0];
        let all_nodes = calc_n_total_nodes(&n_nodes);
        assert_eq!(all_nodes,2);

        let cp=m0_substance.specific_heat_capacity().unwrap();
        let rho=m0_substance.density().unwrap();
        let dx=m0.thickness().unwrap();

        let mut k_prime = Matrix::new(0.0,all_nodes,all_nodes);
        let mut full_rsi=0.0;
        let mut full_rso=0.0;
        let mut c_i=0.0;
        let mut c_o=0.0;        
        build_thermal_network(&building, &c0, dt, &n_nodes, 0.,0., &mut k_prime, &mut full_rsi, &mut full_rso, &mut c_i, &mut c_o).unwrap();

        // check coefficients
        assert_eq!(full_rsi,0.0);// 2.0*dt/(rho*cp*dx));
        assert_eq!(full_rso,m1.thickness().unwrap()/m1_substance.thermal_conductivity().unwrap());// 2.0*dt/(rho*cp*dx));
        assert_eq!(c_i, rho*cp*dx/(2.0*dt));
        assert_eq!(c_o, rho*cp*dx/(2.0*dt));

        let half_mass = m0_substance.density().unwrap() * m0_substance.specific_heat_capacity().unwrap() * m0.thickness().unwrap()/(2.0*dt);
        let u = m0_substance.thermal_conductivity().unwrap()/m0.thickness().unwrap();

        // Check first node        
        let found = k_prime.get(0,0).unwrap();
        let exp = -u / half_mass;
        assert_eq!(exp,found);

        // Check last node
        let rso = m1.thickness().unwrap()/m1_substance.thermal_conductivity().unwrap();
        let found = k_prime.get(1,1).unwrap();
        let exp = -(u + 1.0/rso) / half_mass;        
        assert!((exp-found).abs() < 1E-10);
    }

    #[test]
    fn test_wall_6(){

        let mut building = Building::new("A building".to_string());

        /* SUBSTANCES */
    
        let poly_index = add_polyurethane(&mut building);                
        let brickwork_index = add_brickwork(&mut building);

        /* MATERIALS */

        let m0_thickness = 200.0/1000. as f64;
        let m0_index = add_layer(&mut building, brickwork_index, m0_thickness);
        
        let m1_thickness = 200.0/1000. as f64;
        let m1_index = add_layer(&mut building, poly_index, m1_thickness);
        

        /* WALL */

        let c0_index = building.add_construction("Wall".to_string());
        building.add_material_to_construction(c0_index, m1_index).unwrap();        
        building.add_material_to_construction(c0_index, m0_index).unwrap();        

        let c0 = building.get_construction(c0_index).unwrap();

        let m0 = building.get_material(m0_index).unwrap();
        let m0_substance_index = m0.get_substance_index().unwrap();
        let m0_substance = building.get_substance(m0_substance_index).unwrap();

        let m1 = building.get_material(m1_index).unwrap();
        let m1_substance_index = m1.get_substance_index().unwrap();
        let m1_substance = building.get_substance(m1_substance_index).unwrap();
        

        /* TESTS */
        let dt = 156.0;
        let n_nodes : Vec<usize> = vec![0,1];
        let all_nodes = calc_n_total_nodes(&n_nodes);
        assert_eq!(all_nodes,2);

        let cp=m0_substance.specific_heat_capacity().unwrap();
        let rho=m0_substance.density().unwrap();
        let dx=m0.thickness().unwrap();


        let mut k_prime = Matrix::new(0.0,all_nodes,all_nodes);
        let mut full_rsi=0.0;
        let mut full_rso=0.0;
        let mut c_i=0.0;
        let mut c_o=0.0;                
        build_thermal_network(&building, &c0, dt, &n_nodes, 0.,0., &mut k_prime, &mut full_rsi, &mut full_rso, &mut c_i, &mut c_o).unwrap();

        let rsi = m1.thickness().unwrap()/m1_substance.thermal_conductivity().unwrap();

        // check coefficients
        assert_eq!(full_rsi,rsi);// 2.0*dt/(rho*cp*dx));
        assert_eq!(full_rso,0.0);// 2.0*dt/(rho*cp*dx));
        assert_eq!(c_i, rho*cp*dx/(2.0*dt));
        assert_eq!(c_o, rho*cp*dx/(2.0*dt));

        let half_mass = m0_substance.density().unwrap() * m0_substance.specific_heat_capacity().unwrap() * m0.thickness().unwrap()/(2.0*dt);
        let u = m0_substance.thermal_conductivity().unwrap()/m0.thickness().unwrap();

        // Check first node
    
        let found = k_prime.get(0,0).unwrap();
        let exp = -(u + 1.0/rsi) / half_mass;        
        assert!((exp-found).abs() < 1E-10);

        // Check last node                
        let found = k_prime.get(1,1).unwrap();
        let exp = -u / half_mass;        
        assert!((exp-found).abs() < 1E-10);


        
         
                 
    }


    #[test]
    fn test_wall_7(){


        let mut building = Building::new("A building".to_string());

        /* SUBSTANCES */
        
        let poly_index = add_polyurethane(&mut building);            
        let brickwork_index = add_brickwork(&mut building);
        

        /* MATERIALS */

        let m0_thickness = 200.0/1000. as f64;
        let m0_index = add_layer(&mut building, brickwork_index, m0_thickness);
        
        let m1_thickness = 200.0/1000. as f64;
        let m1_index = add_layer(&mut building, poly_index, m1_thickness);

        /* WALL */

        let c0_index = building.add_construction("Wall".to_string());
        building.add_material_to_construction(c0_index, m0_index).unwrap();        
        building.add_material_to_construction(c0_index, m1_index).unwrap();        
        building.add_material_to_construction(c0_index, m1_index).unwrap();        
        building.add_material_to_construction(c0_index, m0_index).unwrap();        

        let c0 = building.get_construction(c0_index).unwrap();

        let m0 = building.get_material(m0_index).unwrap();
        let m0_substance_index = m0.get_substance_index().unwrap();
        let m0_substance = building.get_substance(m0_substance_index).unwrap();

        let m1 = building.get_material(m1_index).unwrap();
        let m1_substance_index = m1.get_substance_index().unwrap();
        let m1_substance = building.get_substance(m1_substance_index).unwrap();
        


        let cp=m0_substance.specific_heat_capacity().unwrap();
        let rho=m0_substance.density().unwrap();
        let dx=m0.thickness().unwrap();

        let dt = 156.0;
        let n_nodes : Vec<usize> = vec![1,0,0,1];
        let all_nodes = calc_n_total_nodes(&n_nodes);
        assert_eq!(all_nodes,4);


        let mut k_prime = Matrix::new(0.0,all_nodes,all_nodes);
        let mut full_rsi=0.0;
        let mut full_rso=0.0;
        let mut c_i=0.0;
        let mut c_o=0.0;        
        build_thermal_network(&building,&c0, dt, &n_nodes, 0.,0., &mut k_prime, &mut full_rsi, &mut full_rso, &mut c_i, &mut c_o).unwrap();

        // check coefficients
        assert_eq!(full_rsi,0.0);// 2.0*dt/(rho*cp*dx));
        assert_eq!(full_rso,0.0);// 2.0*dt/(rho*cp*dx));
        assert_eq!(c_i, rho*cp*dx/(2.0*dt));
        assert_eq!(c_o, rho*cp*dx/(2.0*dt));

        

        let half_mass = m0_substance.density().unwrap() * m0_substance.specific_heat_capacity().unwrap() * m0.thickness().unwrap()/(2.0*dt);
        let u = m0_substance.thermal_conductivity().unwrap()/m0.thickness().unwrap();
        let insulation_r = m1.thickness().unwrap() * 2.0/m1_substance.thermal_conductivity().unwrap();

        // Check first node        
        let found = k_prime.get(0,0).unwrap();
        let exp = -u / half_mass;        
        assert!((exp-found).abs() < 1E-10);


        // Check second node        
        let found = k_prime.get(1,1).unwrap();
        let exp = -(u + 1.0/insulation_r) / half_mass;        
        assert!((exp-found).abs() < 1E-10);
                
        // Check third node        
        let found = k_prime.get(2,2).unwrap();
        let exp = -(u + 1.0/insulation_r) / half_mass;                
        assert!((exp-found).abs()<1E-6);

        // Check last node        
        let found = k_prime.get(3,3).unwrap();
        let exp = -u / half_mass;        
        assert!((exp-found).abs() < 1E-10);
    }
        

    #[test] 
    fn test_calc_heat_flow_with_mass(){


        let mut building = Building::new("A building".to_string());

        /* SUBSTANCES */
        
        let poly_index = add_polyurethane(&mut building);
        
        /* MATERIALS */

        let m0_thickness = 200.0/1000. as f64;
        let m0_index = add_layer(&mut building, poly_index, m0_thickness);
        

        /* CONSTRUCTION */
        let c0_index = building.add_construction("Wall".to_string());
        building.add_material_to_construction(c0_index, m0_index).unwrap();        


        

        /* GEOMETRY */
        let mut the_loop = Loop3D::new();
        let l = 1. as f64;
        the_loop.push( Point3D::new(-l, -l, 0.)).unwrap();
        the_loop.push( Point3D::new(l, -l, 0.)).unwrap();
        the_loop.push( Point3D::new(l, l, 0.)).unwrap();
        the_loop.push( Point3D::new(-l, l, 0.)).unwrap();
        the_loop.close().unwrap();
        let p = Polygon3D::new(the_loop).unwrap();
        

        /* SURFACE */
        let s_index = building.add_surface("Surface".to_string());
        building.set_surface_polygon(s_index, p).unwrap();
        building.set_surface_construction(s_index, c0_index).unwrap();

        
        /* TEST */
        let surface = building.get_surface(s_index).unwrap();
        let m0 = building.get_material(m0_index).unwrap();
        let c0 = building.get_construction(c0_index).unwrap();

        let main_dt = 300.0;
        let max_dx = m0.thickness().unwrap()/4.0;
        let min_dt = 1.0;
        let (dt,nodes)=find_dt_and_n_nodes(&building,&c0, main_dt, 1, max_dx, min_dt);

        let ts = ThermalSurface::new(&building, &surface, dt,&nodes,0);
        assert!(ts.massive);
        assert_eq!(20.0,ts.state.get(0,0).unwrap());
        assert_eq!(20.0,ts.state.get(ts.n_nodes-1,0).unwrap());

        
        
        let t_in = 10.0;
        let q_in_ref = (20.0-t_in)/ts.rs_i;
        
        let t_out = 10.0;
        let q_out_ref = (20.0-t_out)/ts.rs_o;
        let (q_in,q_out)=ts.calc_heat_flow(t_in,t_out);
        
        assert_eq!(q_in,q_in_ref);
        assert_eq!(q_out,q_out_ref);
        
        
        
    }

    #[test]
    fn test_calc_heat_flow_no_mass(){

        let mut building = Building::new("A building".to_string());

        /* SUBSTANCES */
        
        let poly_index = add_polyurethane(&mut building);

        /* MATERIALS */
        let m0_thickness = 200.0/1000. as f64;
        let m0_index = add_layer(&mut building, poly_index, m0_thickness);
    

        /* CONSTRUCTION */
        let c0_index = building.add_construction("Wall".to_string());
        building.add_material_to_construction(c0_index, m0_index).unwrap();        


        

        /* GEOMETRY */
        let mut the_loop = Loop3D::new();
        let l = 1. as f64;
        the_loop.push( Point3D::new(-l, -l, 0.)).unwrap();
        the_loop.push( Point3D::new(l, -l, 0.)).unwrap();
        the_loop.push( Point3D::new(l, l, 0.)).unwrap();
        the_loop.push( Point3D::new(-l, l, 0.)).unwrap();
        the_loop.close().unwrap();
        let p = Polygon3D::new(the_loop).unwrap();
        

        /* SURFACE */
        let s_index = building.add_surface("Surface".to_string());
        building.set_surface_polygon(s_index, p).unwrap();
        building.set_surface_construction(s_index, c0_index).unwrap();

        /* TEST */
        let surface = building.get_surface(s_index).unwrap();
        let m0 = building.get_material(m0_index).unwrap();
        let c0 = building.get_construction(c0_index).unwrap();


        let main_dt = 300.0;
        let max_dx = m0.thickness().unwrap()/15.;
        let min_dt = 80.;
        let (dt,nodes)=find_dt_and_n_nodes(&building,&c0, main_dt, 1, max_dx, min_dt);
        let ts = ThermalSurface::new(&building,&surface,dt,&nodes,0);
        
        assert!(!ts.massive);
        assert_eq!(20.0,ts.state.get(0,0).unwrap());
        assert_eq!(20.0,ts.state.get(ts.n_nodes-1,0).unwrap());
        
        let t_in = 10.0;
        let q_in_ref = (20.0-t_in)/ts.rs_i;
        
        let t_out = 10.0;
        let q_out_ref = (20.0-t_out)/ts.rs_o;
        let (q_in,q_out)=ts.calc_heat_flow(t_in,t_out);
        
        assert_eq!(q_in,q_in_ref);
        assert_eq!(q_out,q_out_ref);
    }

    #[test]
    fn test_calc_heat_flow_mixed_mass(){

        // THIRD TEST -- WITH ONE NO-MASS LAYER IN THE EXTERIOR AND ONE IN THE INTERIOR

        let mut building = Building::new("Some building".to_string());
        
        /* SUBSTANCES */
        let poly_index = add_polyurethane(&mut building);
        let brickwork_index = add_brickwork(&mut building);

        /* MATERIALS */
        let m1_index = add_layer(&mut building, poly_index, 20./1000.);
        let m2_index = add_layer(&mut building, brickwork_index, 220./1000.);

        /* CONSTRUCTION */

        let c_index = building.add_construction("construction".to_string());
        building.add_material_to_construction(c_index, m1_index).unwrap();
        building.add_material_to_construction(c_index, m2_index).unwrap();
        building.add_material_to_construction(c_index, m1_index).unwrap();

        

        /* GEOMETRY */
        
        let mut the_loop = Loop3D::new();
        let l = 1. as f64;
        the_loop.push( Point3D::new(-l, -l, 0.)).unwrap();
        the_loop.push( Point3D::new(l, -l, 0.)).unwrap();
        the_loop.push( Point3D::new(l, l, 0.)).unwrap();
        the_loop.push( Point3D::new(-l, l, 0.)).unwrap();
        the_loop.close().unwrap();
        let p = Polygon3D::new(the_loop).unwrap();

        /* SURFACE */

        let s_index = building.add_surface("Surface 1".to_string());
        building.set_surface_construction(s_index, c_index).unwrap();
        building.set_surface_polygon(s_index, p).unwrap();


        /* TEST */
        let c = building.get_construction(c_index).unwrap();
        let surface = building.get_surface(s_index).unwrap();

        let main_dt = 300.0;
        let max_dx = 0.015;
        let min_dt = 65.;
        let (dt,nodes)=find_dt_and_n_nodes(&building,&c, main_dt, 1, max_dx, min_dt);
        let ts = ThermalSurface::new(&building,&surface,dt,&nodes,0);
        
        assert!(ts.massive);
        assert_eq!(20.0,ts.state.get(0,0).unwrap());
        assert_eq!(20.0,ts.state.get(ts.n_nodes-1,0).unwrap());
        
        let t_in = 10.0;
        let q_in_ref = (20.0-t_in)/ts.full_rsi;
        
        let t_out = 10.0;
        let q_out_ref = (20.0-t_out)/ts.full_rso;
        let (q_in,q_out)=ts.calc_heat_flow(t_in,t_out);
        
        assert_eq!(q_in,q_in_ref);
        assert_eq!(q_out,q_out_ref);
    }

    #[test]
    fn test_march_massive(){
        let mut building = Building::new("Some building".to_string());
        
        /* SUBSTANCES */        
        let brickwork_index = add_brickwork(&mut building);

        /* MATERIALS */
        let m1_index = add_layer(&mut building, brickwork_index, 200./1000.);
        

        /* CONSTRUCTION */

        let c_index = building.add_construction("construction".to_string());
        building.add_material_to_construction(c_index, m1_index).unwrap();
        

        /* GEOMETRY */
        let mut the_loop = Loop3D::new();
        let l = 1. as f64;
        the_loop.push( Point3D::new(-l, -l, 0.)).unwrap();
        the_loop.push( Point3D::new(l, -l, 0.)).unwrap();
        the_loop.push( Point3D::new(l, l, 0.)).unwrap();
        the_loop.push( Point3D::new(-l, l, 0.)).unwrap();
        the_loop.close().unwrap();
        let p = Polygon3D::new(the_loop).unwrap();

        /* SURFACE */
        let s_index = building.add_surface("Surface 1".to_string());
        building.set_surface_polygon(s_index, p).unwrap();
        building.set_surface_construction(s_index, c_index).unwrap();

        /* TEST */
        let m1 = building.get_material(m1_index).unwrap();
        let c = building.get_construction(c_index).unwrap();
        let surface = building.get_surface(s_index).unwrap();

        // FIRST TEST -- 10 degrees on each side
        let main_dt = 300.0;
        let max_dx = m1.thickness().unwrap()/2.0;
        let min_dt = 1.0;
        let (dt,nodes)=find_dt_and_n_nodes(&building,&c, main_dt, 1, max_dx, min_dt);
        let mut ts = ThermalSurface::new(&building, surface, dt, &nodes,0);
        assert!(ts.massive);
       
        // Try marching until q_in and q_out are zero.
        let mut q : f64 = 9000009.0;

        
        let mut counter : usize = 0;
        while q.abs() > 1E-5 {
            let (q_in,q_out)=ts.march(10.0,10.0);                        
            assert!((q_in-q_out).abs() < 1E-5);            
            assert!(q_in  >= 0.);            
            assert!(q_in  < q);            
            q = q_in;            

            counter += 1;
            if counter > 99999 {
                panic!("Exceded number of iterations")
            }
        }

        // all nodes should be at 10.0 now.
        for i in 0..ts.n_nodes{
            let t = ts.state.get(i,0).unwrap();
            assert!((t-10.0).abs()<1E-5);
        }

        // SECOND TEST -- 10 degrees in and 30 out.
        // We expect the final q to be (30-10)/R from
        // outside to inside. I.E. q_in = (30-10)/R, 
        // q_out = -(30-10)/R
        
        // March until q converges
        let mut change : f64 = 99.0;
        let mut counter : usize = 0;
        let mut previous_q : f64 = -125.0;
        let mut final_qin:f64 = -12312.;
        let mut final_qout:f64 = 123123123.;
        while change.abs() > 1E-10 {
            
            let (q_in,q_out)=ts.march(10.0,30.0);                        
            final_qin = q_in;
            final_qout = q_out;
            
            change = (q_in - previous_q).abs();
            previous_q = q_in;

            counter += 1;
            if counter > 99999 {
                panic!("Exceded number of iterations")
            }
        }

        // Expecting 
        assert!(final_qin > 0.0);
        assert!(final_qout < 0.0);
        assert!((final_qin+final_qout).abs()<1E-6);

        let r = r_value(&building, c).unwrap();

        let exp_q = (30.0 - 10.0)/(r + ts.rs_i + ts.rs_o);
        assert!((exp_q-final_qin).abs() < 1E-4);
        assert!((exp_q+final_qout).abs() < 1E-4);
    }


    #[test]
    fn test_march_nomass(){

        let mut building = Building::new("A building".to_string());

        /* SUBSTANCE */
        let brickwork_index = add_brickwork(&mut building);

        /* MATERIAL */
        let m1_index = add_layer(&mut building, brickwork_index, 3./1000.);

        /* CONSTRUCTION */
        let c_index = building.add_construction("Construction".to_string());
        building.add_material_to_construction(c_index, m1_index).unwrap();

        /* GEOMETRY */
        let mut the_loop = Loop3D::new();
        let l = 1. as f64;
        the_loop.push( Point3D::new(-l, -l, 0.)).unwrap();
        the_loop.push( Point3D::new(l, -l, 0.)).unwrap();
        the_loop.push( Point3D::new(l, l, 0.)).unwrap();
        the_loop.push( Point3D::new(-l, l, 0.)).unwrap();
        the_loop.close().unwrap();
        let p = Polygon3D::new(the_loop).unwrap();

        /* SURFACE */
        let s_index = building.add_surface("WALL".to_string());
        building.set_surface_construction(s_index, c_index).unwrap();
        building.set_surface_polygon(s_index, p).unwrap();

        /* TEST */

        let m1 = building.get_material(m1_index).unwrap();
        let c = building.get_construction(c_index).unwrap();
        let surface = building.get_surface(s_index).unwrap();        


        // FIRST TEST -- 10 degrees on each side
        let main_dt = 300.0;
        let max_dx = m1.thickness().unwrap()/2.0;
        let min_dt = 100.0;
        let (dt,nodes)=find_dt_and_n_nodes(&building, c, main_dt, 1, max_dx, min_dt);
        let mut ts = ThermalSurface::new(&building,surface,dt,&nodes,0);
        assert!(!ts.massive);
       
        // Try marching until q_in and q_out are zero.

        let (q_in,q_out)= ts.march(10.0,10.0);                        

        // this should show instantaneous update. So,
        assert_eq!(ts.state.get(0,0).unwrap(),10.0);
        assert_eq!(ts.state.get(1,0).unwrap(),10.0);
        assert_eq!(q_in,0.0);
        assert_eq!(q_out,0.0);
                

        // SECOND TEST -- 10 degrees in and 30 out.
        // We expect the final q to be (30-10)/R from
        // outside to inside. I.E. q_in = (30-10)/R, 
        // q_out = -(30-10)/R
        
        let (q_in,q_out)=ts.march(10.0,30.0);                        

       
        // Expecting 
        assert!(q_in > 0.0);
        assert!(q_out < 0.0);
        assert!((q_in+q_out).abs()<1E-6);

        let exp_q = (30.0 - 10.0)/(r_value(&building,c).unwrap() + ts.rs_i+ts.rs_o);
        assert!((exp_q-q_in).abs() < 1E-4);
        assert!((exp_q+q_out).abs() < 1E-4);
    }
    
}



    
   

