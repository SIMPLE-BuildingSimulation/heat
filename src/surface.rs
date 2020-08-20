use std::rc::Rc;

use matrix::Matrix;
use building::surface::Surface;
use building::construction::Construction;
use convection::*;


use crate::zone::ThermalZone;

/// This is a Surface from the point of 
/// view of our thermal solver.
pub struct ThermalSurface {
    
    /// A reference to the original Surface in the 
    /// building model
    surface : Rc<Surface>,

    /// Area of the surface
    area : f64,
    
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

    /// The position of this surface within 
    /// the Thermal Model surfaces array
    position: usize,
}

impl ThermalSurface {

    pub fn new(surface: Rc<Surface>, dt: f64, nodes: &Vec<usize>, position: usize ) -> Self {
                
        let n_nodes = calc_n_total_nodes(&nodes);
        let (rs_i,rs_o) = calc_convection_coefficients(surface.as_ref());
        let area = surface.area();
        
        let mut ret = ThermalSurface{
            surface: surface,
            area: area,
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
            position: position
        };
        
        build_thermal_network(&ret.surface.construction(), dt, &nodes, rs_i, rs_o, &mut ret.k_prime, &mut ret.full_rsi, &mut ret.full_rso, &mut ret.c_i, &mut ret.c_o).unwrap();
        ret.massive = !(ret.c_o == 0.0 && ret.c_i == 0.);

        // return
        ret
    }
    

    pub fn area(&self)->f64{
        self.surface.area()
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
            //                        lets call this vector 'a'         These arethe F components        
            let mut a = self.k_prime.from_prod_n_diag(&self.state,3).unwrap();// k_prime is tri-diagonal
            // ... a should be a vector

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

fn build_thermal_network(c: &Construction, dt: f64, n_nodes: &Vec<usize>, rs_i: f64, rs_o: f64, k_prime: & mut Matrix, full_rsi: &mut f64, full_rso: &mut f64, c_i: & mut f64, c_o: & mut f64 ) ->Result<(),String> {
    
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
        let r = c.r_value()+rs_i+rs_o;
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
            let layer = c.layer(i).unwrap();
            *full_rsi += layer.thickness()/layer.substance().thermal_conductivity();
        }
        *full_rso = rs_o;
        for i in last_massive..c.n_layers() {        
            let layer = c.layer(i).unwrap();
            *full_rso += layer.thickness()/layer.substance().thermal_conductivity();
        }
    }

    
    // Calculate the rest.
    let mut node : usize = 0;
    let mut n_layer : usize = first_massive;

    while n_layer < last_massive {
        let material = c.layer(n_layer).unwrap();    
        let m = n_nodes[n_layer];                
        if m == 0 { 
            
            // nomass material
            // add up all the R of the no-mass layers that
            // are together            
            let mut r = 0.0; // if the material is no mass, then the first value 
            while n_layer < last_massive && n_nodes[n_layer] == 0 {                
                let layer = c.layer(n_layer).unwrap();
                let dx = layer.thickness();
                let k =layer.substance().thermal_conductivity();
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
            let k = material.substance().thermal_conductivity();                
            let dx = material.thickness()/(m as f64);            
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
        let material = c.layer(n_layer).unwrap();
        let m = n_nodes[n_layer];  

        if m != 0 {
            for _i in 0..m {                   
                // Calc mass
                let rho = material.substance().density();
                let cp = material.substance().heat_capacity();
                let dx = material.thickness()/(m as f64);
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

pub fn find_dt_and_n_nodes(c: &Construction, main_dt: f64, n: usize, max_dx: f64, min_dt: f64) -> (f64,Vec<usize>){

    let dt = main_dt/(n as f64);
    let safety = 1.5;

    // Choose a dx so that dt allows convergence.
    // stability is assured by (alpha * dt / dx^2 <= 1/2 )
    // meaning, we need to satisfy dx >= sqrt(2 * alpha * dt)
    
    // So, for each layer
    let mut n_nodes : Vec<usize> = vec![];

    for n_layer in 0..c.n_layers(){
        let layer = c.layer(n_layer).unwrap();

        // Calculate the optimum_dx
        let alpha = layer.substance().thermal_diffusivity();
        let optimum_dx = (safety*2.0*alpha*dt).sqrt();
        let m = (layer.thickness()/optimum_dx).floor();            
        let mut dx = layer.thickness()/m;
        // dx cannot be smaller than thickness
        if m == 0. {
            dx = layer.thickness();
        }        

        // If dx is larger than the max allowed d_x, try to change timestep        
        if dx > max_dx {
            // check if there is room for reducing dt...
            let next_dt = main_dt/((n+1) as f64);
            if next_dt > min_dt {
                // If there is room for that, do it.
                return find_dt_and_n_nodes(c,main_dt,n+1,max_dx,min_dt)
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
    use building::material::Material;
    use building::substance::Substance;
    use geometry3d::polygon3d::Polygon3D;
    use geometry3d::point3d::Point3D;
    use geometry3d::loop3d::Loop3D;

    #[test]
    fn test_get_dt_and_n_nodes(){
        
        
        let polyurethane = Substance::new(
            "polyurethane".to_string(),
            0.0252, // W/m.K            
            2400., // J/kg.K
            17.5, // kg/m3... reverse engineered from paper            
        );
        assert_eq!(polyurethane.thermal_diffusivity(),0.6E-6);

        let brickwork = Substance::new(
            "brickwork".to_string(),
            0.816, // W/m.K            
            800., // J/kg.K
            1700., // kg/m3... reverse engineered from paper            
        );
        assert_eq!(brickwork.thermal_diffusivity(),0.6E-6);

        /* WALL 1 */
        let m1 = Material::new(Rc::clone(&polyurethane), 200.0/1000.);        
        let c = Construction::new("wall 1".to_string(),vec![Rc::clone(&m1)]);
        

        // This should result in a dt of 300, with 8 layers of 0.25m thick.
        let main_dt = 300.;
        let (dt,n_nodes) = find_dt_and_n_nodes(&c,main_dt,1,m1.thickness()/4.,1.);
        assert_eq!(n_nodes.len(),1);
        assert_eq!(dt,300.);
        assert_eq!(n_nodes[0],8);
        

        // This should result in a dt of 300, with 8 layers of 0.25m thick.
        let main_dt = 300.;
        let (dt,n_nodes) = find_dt_and_n_nodes(&c,main_dt,1,m1.thickness()/8.,1.);
        assert_eq!(n_nodes.len(),1);
        assert_eq!(dt,300.);
        assert_eq!(n_nodes[0],8);
        

        // This should result in a dt of 300/2=150, with 12 layers of 0.0166667m thick.
        let main_dt = 300.;
        let (dt,n_nodes) = find_dt_and_n_nodes(&c,main_dt,1,m1.thickness()/9.,1.);
        assert_eq!(n_nodes.len(),1);
        assert_eq!(dt,150.);
        assert_eq!(n_nodes[0],12);
        

        // This should result in a dt of 300/4=75, with 17 layers of 0.011...m thick.
        let main_dt = 300.;
        let (dt,n_nodes) = find_dt_and_n_nodes(&c,main_dt,1,m1.thickness()/15.,1.);
        assert_eq!(n_nodes.len(),1);
        assert_eq!(dt,75.);
        assert_eq!(n_nodes[0],17);
        

        // We imposed a min_dt of 80, meaning that this should result in a no-mass layer
        let main_dt = 300.;
        let (dt,n_nodes) = find_dt_and_n_nodes(&c,main_dt,1,m1.thickness()/15.,80.);
        assert_eq!(n_nodes.len(),1);
        assert_eq!(dt,100.);
        assert_eq!(n_nodes[0],0);
        

        /* WALL 2 */
        let m1 = Material::new(Rc::clone(&polyurethane),200.0/1000.);
        let m2 = Material::new(Rc::clone(&brickwork),110.0/1000.);
        let c = Construction::new("wall 2".to_string(), vec![Rc::clone(&m1),Rc::clone(&m2)]);
        

        // This should result in a dt of 300, with 8 layers and 4 layers.
        let main_dt = 300.;
        let (dt,n_nodes) = find_dt_and_n_nodes(&c,main_dt,1,0.03,1.);
        assert_eq!(n_nodes.len(),2);
        assert_eq!(dt,300.);
        assert_eq!(n_nodes[0],8);
        assert_eq!(n_nodes[1],4);
        
        

        // This should result in a dt of 300/2=150, with 12 and 6 layers.
        let main_dt = 300.;
        let (dt,n_nodes) = find_dt_and_n_nodes(&c,main_dt, 1 , 0.02,1.);
        assert_eq!(n_nodes.len(),2);
        assert_eq!(dt,150.);
        assert_eq!(n_nodes[0],12);
        assert_eq!(n_nodes[1],6);
        

        // This should result in a dt of 300/3=100, with 14 and 8
        let main_dt = 300.;
        let (dt,n_nodes) = find_dt_and_n_nodes(&c,main_dt,1,0.015,1.);
        assert_eq!(n_nodes.len(),2);
        assert_eq!(dt,100.);
        assert_eq!(n_nodes[0],14);
        assert_eq!(n_nodes[1],8);
        

        // We imposed a min_dt of 100, meaning that this should result in a no-mass layer
        let main_dt = 300.;
        let (dt,n_nodes) = find_dt_and_n_nodes(&c,main_dt,1,0.015,110.);
        assert_eq!(n_nodes.len(),2);
        assert_eq!(dt,150.);
        assert_eq!(n_nodes[0],0);
        assert_eq!(n_nodes[1],0);
        

        /* WALL 3 */
        let m1 = Material::new(Rc::clone(&polyurethane),20.0/1000.);
        let m2 = Material::new(Rc::clone(&brickwork),220.0/1000.);        
        let m3 = Rc::clone(&m1);
        let c = Construction::new("wall 3".to_string(), vec![Rc::clone(&m1),Rc::clone(&m2),Rc::clone(&m3)]);
        

        // This should result in a dt of 150, with [1,13,1] layers
        let main_dt = 300.;
        let (dt,n_nodes) = find_dt_and_n_nodes(&c, main_dt, 1, 0.03, 1.);
        assert_eq!(n_nodes.len(),3);
        assert_eq!(dt,300.);
        assert_eq!(n_nodes[0],0);
        assert_eq!(n_nodes[1],9);
        assert_eq!(n_nodes[2],0);
        
        

        // This should result in a dt of 300/6=50, with [2,23,2] layers
        let main_dt = 300.;
        let (dt,n_nodes) = find_dt_and_n_nodes(&c,main_dt, 1 , 0.015,1.);
        assert_eq!(n_nodes.len(),3);
        assert_eq!(dt,50.);
        assert_eq!(n_nodes[0],2);
        assert_eq!(n_nodes[1],23);
        assert_eq!(n_nodes[2],2);
        

        // Limit min_time to 65... This should result in a dt of 300/4=75, with [0, 18, 0] layers
        let main_dt = 300.;
        let (dt,n_nodes) = find_dt_and_n_nodes(&c,main_dt, 1 , 0.015,65.);
        assert_eq!(n_nodes.len(),3);
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
    fn test_wall_1(){
       
        let brickwork = Substance::new(
            "brickwork".to_string(),
            0.816, // W/m.K            
            800., // J/kg.K
            1700., // kg/m3... reverse engineered from paper            
        ); 
        assert_eq!(brickwork.thermal_diffusivity(),0.6E-6);

        

        /* WALL 1: Single layer, with mass. */
        let m1 = Material::new(Rc::clone(&brickwork),200.0/1000.);
        let c = Construction::new("wall 1".to_string(), vec![Rc::clone(&m1)]);
            
        
    
        let dt = 156.0;
        let n_nodes : Vec<usize> = vec![4];
        let all_nodes = calc_n_total_nodes(&n_nodes);

        let mut k_prime = Matrix::new(0.0,all_nodes,all_nodes);
        let mut full_rsi=0.0;
        let mut full_rso=0.0;
        let mut c_i=0.0;
        let mut c_o=0.0;        
        build_thermal_network(&c, dt, &n_nodes, 0.,0., &mut k_prime, &mut full_rsi, &mut full_rso, &mut c_i, &mut c_o).unwrap();

        let n_layer = 0;
        let rho = m1.substance().density();
        let cp = m1.substance().heat_capacity();
        let m = n_nodes[n_layer];
        let dx = m1.thickness()/(m as f64);
        let k = m1.substance().thermal_conductivity();
        let alpha = m1.substance().thermal_diffusivity();
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
        //assert!( (found-exp).abs() < 1E-10 );        

        let found = k_prime.get(node,node - 1).unwrap();
        let exp = 2.0*k/dx/mass;
        assert!( (found-exp).abs() < 1E-10 );

    }


    #[test]
    fn test_wall_2(){

        let polyurethane = Substance::new(
            "polyurethane".to_string(),
            0.0252, // W/m.K            
            2400., // J/kg.K
            17.5, // kg/m3... reverse engineered from paper            
        ); 
        assert_eq!(polyurethane.thermal_diffusivity(),0.6E-6);

        let brickwork = Substance::new(
            "brickwork".to_string(),
            0.816, // W/m.K            
            800., // J/kg.K
            1700., // kg/m3... reverse engineered from paper            
        ); 
        assert_eq!(brickwork.thermal_diffusivity(),0.6E-6);

        
        /* WALL 2 : Two materials with mass */
        let m1 = Material::new(Rc::clone(&polyurethane),200.0/1000.);
        let m2 = Material::new(Rc::clone(&brickwork),110.0/1000.);
        let c = Construction::new("wall 2".to_string(), vec![Rc::clone(&m1), Rc::clone(&m2)]);            


        let dt = 156.0;
        let n_nodes : Vec<usize> = vec![3,3];
        let all_nodes = calc_n_total_nodes(&n_nodes);

        let cp1=m1.substance().heat_capacity();
        let rho1=m1.substance().density();
        let dx1=m1.thickness()/3.;
        let mass1 = rho1*cp1*dx1/dt;

        let cp2=m2.substance().heat_capacity();
        let rho2=m2.substance().density();
        let dx2=m2.thickness()/3.;
        let mass2 = rho2*cp2*dx2/dt;

        let mut k_prime = Matrix::new(0.0,all_nodes,all_nodes);
        let mut full_rsi=0.0;
        let mut full_rso=0.0;
        let mut c_i=0.0;
        let mut c_o=0.0;        
        build_thermal_network(&c, dt, &n_nodes, 0.,0., &mut k_prime, &mut full_rsi, &mut full_rso, &mut c_i, &mut c_o).unwrap();

        // check coefficients
        assert_eq!(full_rsi,0.0);// 2.0*dt/(rho*cp*dx));
        assert_eq!(full_rso,0.0);// 2.0*dt/(rho*cp*dx));
        assert_eq!(c_i, rho1*cp1*dx1/(2.0*dt));
        assert_eq!(c_o, rho2*cp2*dx2/(2.0*dt));

        // FIRST LAYER
        //////////////
        let n_layer = 0;
        //let rho = m1.substance.density;
        //let cp = m1.substance.heat_capacity;
        let m = n_nodes[n_layer];
        let dx = m1.thickness()/(m as f64);
        let k = m1.substance().thermal_conductivity();
        let alpha = m1.substance().thermal_diffusivity();

        // Check first node in layer
        let node = 0;        
        
        let found = k_prime.get(node,node).unwrap();
        let exp = -2.0*alpha*dt/dx/dx;
        assert!( (found-exp).abs() < 1E-10 );        

        let found = k_prime.get(node,node + 1).unwrap();
        let exp = 2.0*k/dx/mass1;
        assert!( (found-exp).abs() < 1E-10 );

        // check middle nodes
        for node in 1..3{                        

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
        // check middle node (e.g. node 3), which 
        // is also the beginning of Layer 2
        let node = 3;
        let mass_1 = m1.substance().density() * m1.substance().heat_capacity()*m1.thickness()/(6.0*dt);        
        let mass_2 = m2.substance().density() * m2.substance().heat_capacity()*m2.thickness()/(6.0*dt);        
        
        let found = k_prime.get(node,node).unwrap();
        let u1 = 3.0*m1.substance().thermal_conductivity()/m1.thickness();
        let u2 = 3.0*m2.substance().thermal_conductivity()/m2.thickness();
        let exp = -(u1+u2)/(mass_1+mass_2);        
        assert_eq!(found,exp);
        
        let found = k_prime.get(node,node - 1).unwrap();        
        assert_eq!(found, u1/(mass1/2. + mass2/2.0));

        let found = k_prime.get(node,node + 1).unwrap();        
        assert_eq!(found, u2/(mass1/2.0+mass2/2.0));

        // SECOND LAYER
        //////////////
        let n_layer = 1;
        //let rho = m2.substance.density;
        //let cp = m2.substance.heat_capacity;
        let m = n_nodes[n_layer];
        let dx = m2.thickness()/(m as f64);
        let k = m2.substance().thermal_conductivity();
        let alpha = m2.substance().thermal_diffusivity();

        
        // check middle nodes in layer 2
        for node in 4..6{                        

            let found = k_prime.get(node,node).unwrap();
            let exp = -2.0*alpha*dt/dx/dx;
            assert!( (found-exp).abs() < 1E-10 );

            let found = k_prime.get(node,node + 1).unwrap();
            let exp = k/dx/mass2;
            assert!( (found-exp).abs() < 1E-10 );

            let found = k_prime.get(node,node - 1).unwrap();
            let exp = k/dx/mass2;
            assert!( (found-exp).abs() < 1E-10 );
        }
        // check end node 
        let node = 6;        
        
        
        let found = k_prime.get(node,node).unwrap();        
        let u2 = 3.0*m2.substance().thermal_conductivity()/m2.thickness();
        let exp = -2.0*u2/(mass2);        
        assert_eq!(found,exp);
        
        let found = k_prime.get(node,node - 1).unwrap();        
        assert_eq!(found, 2.0*u2/mass2);
    }

    #[test]
    fn test_wall_3(){

         let polyurethane = Substance::new(
            "polyurethane".to_string(),
            0.0252, // W/m.K            
            2400., // J/kg.K
            17.5, // kg/m3... reverse engineered from paper            
         ); 

        assert_eq!(polyurethane.thermal_diffusivity(),0.6E-6);

       
        
        
        /* WALL 3 : Single layer, no mass */
        let m1 = Material::new(Rc::clone(&polyurethane),200.0/1000.);
        let c = Construction::new("wall 3".to_string(), vec![Rc::clone(&m1)]);            
        

        let dt = 156.0;
        let n_nodes : Vec<usize> = vec![0];
        let all_nodes = calc_n_total_nodes(&n_nodes);
        let mut k_prime = Matrix::new(0.0,all_nodes,all_nodes);
        
        let mut full_rsi=0.0;
        let mut full_rso=0.0;
        let mut c_i=0.0;
        let mut c_o=0.0;        
        build_thermal_network(&c, dt, &n_nodes, 0.,0., &mut k_prime, &mut full_rsi, &mut full_rso, &mut c_i, &mut c_o).unwrap();
        
        // check coefficients
        assert_eq!(full_rsi, m1.thickness()/m1.substance().thermal_conductivity());// 2.0*dt/(rho*cp*dx));
        assert_eq!(full_rso, m1.thickness()/m1.substance().thermal_conductivity());// 2.0*dt/(rho*cp*dx));
        assert_eq!(c_i, 0.0);
        assert_eq!(c_o, 0.0);     
        
        for row in 0..all_nodes{
            for col in 0..all_nodes{
                let exp = m1.substance().thermal_conductivity()/m1.thickness();            
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
    fn test_wall_4(){

        let polyurethane = Substance::new(
            "polyurethane".to_string(),
            0.0252, // W/m.K            
            2400., // J/kg.K
            17.5, // kg/m3... reverse engineered from paper            
        ); 
        assert_eq!(polyurethane.thermal_diffusivity(),0.6E-6);

        /* WALL 4 : two equal layers, both with no mass */
        let m1 = Material::new(Rc::clone(&polyurethane), 200.0/1000.);
        let c = Construction::new("wall 4".to_string(),vec![Rc::clone(&m1), Rc::clone(&m1)]);

        
        
        let dt = 156.0;
        let n_nodes : Vec<usize> = vec![0,0];
        let all_nodes = calc_n_total_nodes(&n_nodes);
        let r = 2.0*m1.thickness()/m1.substance().thermal_conductivity();
        let mut k_prime = Matrix::new(0.0,all_nodes,all_nodes);
        let mut full_rsi=0.0;
        let mut full_rso=0.0;
        let mut c_i=0.0;
        let mut c_o=0.0;                
        build_thermal_network(&c, dt, &n_nodes, 0.,0., &mut k_prime, &mut full_rsi, &mut full_rso, &mut c_i, &mut c_o).unwrap();

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

        let polyurethane = Substance::new(
            "polyurethane".to_string(),
            0.0252, // W/m.K            
            2400., // J/kg.K
            17.5, // kg/m3... reverse engineered from paper            
        ); 
        
        assert_eq!(polyurethane.thermal_diffusivity(),0.6E-6);

        let brickwork = Substance::new(
            "brickwork".to_string(),
            0.816, // W/m.K            
            800., // J/kg.K
            1700., // kg/m3... reverse engineered from paper            
        ); 
        assert_eq!(brickwork.thermal_diffusivity(),0.6E-6);

        /* WALL 5 : One layer with mass, one layer without it */
        
        let m1 = Material::new(Rc::clone(&brickwork),200./1000.);
        let m2 = Material::new(Rc::clone(&polyurethane),200./1000.);        
        let c = Construction::new("wall 5".to_string(),vec![Rc::clone(&m1), Rc::clone(&m2)]);

        let dt = 156.0;
        let n_nodes : Vec<usize> = vec![1,0];
        let all_nodes = calc_n_total_nodes(&n_nodes);
        assert_eq!(all_nodes,2);

        let cp=m1.substance().heat_capacity();
        let rho=m1.substance().density();
        let dx=m1.thickness();

        let mut k_prime = Matrix::new(0.0,all_nodes,all_nodes);
        let mut full_rsi=0.0;
        let mut full_rso=0.0;
        let mut c_i=0.0;
        let mut c_o=0.0;        
        build_thermal_network(&c, dt, &n_nodes, 0.,0., &mut k_prime, &mut full_rsi, &mut full_rso, &mut c_i, &mut c_o).unwrap();

        // check coefficients
        assert_eq!(full_rsi,0.0);// 2.0*dt/(rho*cp*dx));
        assert_eq!(full_rso,m2.thickness()/m2.substance().thermal_conductivity());// 2.0*dt/(rho*cp*dx));
        assert_eq!(c_i, rho*cp*dx/(2.0*dt));
        assert_eq!(c_o, rho*cp*dx/(2.0*dt));

        let half_mass = m1.substance().density() * m1.substance().heat_capacity() * m1.thickness()/(2.0*dt);
        let u = m1.substance().thermal_conductivity()/m1.thickness();

        // Check first node        
        let found = k_prime.get(0,0).unwrap();
        let exp = -u / half_mass;
        assert_eq!(exp,found);

        // Check last node
        let rso = m2.thickness()/m2.substance().thermal_conductivity();
        let found = k_prime.get(1,1).unwrap();
        let exp = -(u + 1.0/rso) / half_mass;        
        assert!((exp-found).abs() < 1E-10);
    }

    #[test]
    fn test_wall_6(){

        let polyurethane = Substance::new(
            "polyurethane".to_string(),
            0.0252, // W/m.K            
            2400., // J/kg.K
            17.5, // kg/m3... reverse engineered from paper            
        ); 
        assert_eq!(polyurethane.thermal_diffusivity(),0.6E-6);

        let brickwork = Substance::new(
            "brickwork".to_string(),
            0.816, // W/m.K            
            800., // J/kg.K
            1700., // kg/m3... reverse engineered from paper            
        );
        assert_eq!(brickwork.thermal_diffusivity(),0.6E-6);

                

        /* WALL 6 : One layer without mass, one layer with it */
        
        let m1 = Material::new(Rc::clone(&brickwork),200./1000.);
        let m2 = Material::new(Rc::clone(&polyurethane),200./1000.);
        

        let c = Construction::new("wall 6".to_string(),vec![Rc::clone(&m2), Rc::clone(&m1)]);

        let dt = 156.0;
        let n_nodes : Vec<usize> = vec![0,1];
        let all_nodes = calc_n_total_nodes(&n_nodes);
        assert_eq!(all_nodes,2);

        let cp=m1.substance().heat_capacity();
        let rho=m1.substance().density();
        let dx=m1.thickness();


        let mut k_prime = Matrix::new(0.0,all_nodes,all_nodes);
        let mut full_rsi=0.0;
        let mut full_rso=0.0;
        let mut c_i=0.0;
        let mut c_o=0.0;                
        build_thermal_network(&c, dt, &n_nodes, 0.,0., &mut k_prime, &mut full_rsi, &mut full_rso, &mut c_i, &mut c_o).unwrap();

        let rsi = m2.thickness()/m2.substance().thermal_conductivity();

        // check coefficients
        assert_eq!(full_rsi,rsi);// 2.0*dt/(rho*cp*dx));
        assert_eq!(full_rso,0.0);// 2.0*dt/(rho*cp*dx));
        assert_eq!(c_i, rho*cp*dx/(2.0*dt));
        assert_eq!(c_o, rho*cp*dx/(2.0*dt));

        let half_mass = m1.substance().density() * m1.substance().heat_capacity() * m1.thickness()/(2.0*dt);
        let u = m1.substance().thermal_conductivity()/m1.thickness();

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

        let polyurethane = Substance::new(
            "polyurethane".to_string(),
            0.0252, // W/m.K            
            2400., // J/kg.K
            17.5, // kg/m3... reverse engineered from paper            
        ); 
        assert_eq!(polyurethane.thermal_diffusivity(),0.6E-6);

        let brickwork = Substance::new(
            "brickwork".to_string(),
            0.816, // W/m.K            
            800., // J/kg.K
            1700., // kg/m3... reverse engineered from paper            
        ); 
        assert_eq!(brickwork.thermal_diffusivity(),0.6E-6);

        /* WALL 7 : One layer with mass, two layers without it, and one with it again */
        
        let m1 = Material::new(Rc::clone(&brickwork),200./1000.);
        let m2 = Material::new(Rc::clone(&polyurethane),200./1000.);
        

        let c = Construction::new("wall 7".to_string(),vec![Rc::clone(&m1), Rc::clone(&m2), Rc::clone(&m2), Rc::clone(&m1)]);

        let cp=m1.substance().heat_capacity();
        let rho=m1.substance().density();
        let dx=m1.thickness();

        let dt = 156.0;
        let n_nodes : Vec<usize> = vec![1,0,0,1];
        let all_nodes = calc_n_total_nodes(&n_nodes);
        assert_eq!(all_nodes,4);


        let mut k_prime = Matrix::new(0.0,all_nodes,all_nodes);
        let mut full_rsi=0.0;
        let mut full_rso=0.0;
        let mut c_i=0.0;
        let mut c_o=0.0;        
        build_thermal_network(&c, dt, &n_nodes, 0.,0., &mut k_prime, &mut full_rsi, &mut full_rso, &mut c_i, &mut c_o).unwrap();

        // check coefficients
        assert_eq!(full_rsi,0.0);// 2.0*dt/(rho*cp*dx));
        assert_eq!(full_rso,0.0);// 2.0*dt/(rho*cp*dx));
        assert_eq!(c_i, rho*cp*dx/(2.0*dt));
        assert_eq!(c_o, rho*cp*dx/(2.0*dt));

        

        let half_mass = m1.substance().density() * m1.substance().heat_capacity() * m1.thickness()/(2.0*dt);
        let u = m1.substance().thermal_conductivity()/m1.thickness();
        let insulation_r = m2.thickness() * 2.0/m2.substance().thermal_conductivity();

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
    fn test_new(){

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

        let m1 = Material::new(Rc::clone(&polyurethane),200./1000.);

        let c = Construction::new("wall 1".to_string(),vec![Rc::clone(&m1)]);

        // Physical surface
        let surface = Surface::new(p,Rc::clone(&c));        

        // This should result in a dt of 300, with 8 layers of 0.25m thick.
        let main_dt = 300.0;
        let max_dx = m1.thickness()/4.0;
        let min_dt = 1.;
        let (dt,nodes)=find_dt_and_n_nodes(&c, main_dt, 1, max_dx, min_dt);
        let ts = ThermalSurface::new(Rc::clone(&surface),dt,&nodes,0);

        let (rs_i,rs_o)=calc_convection_coefficients(&surface);
        assert!(ts.massive);
        assert_eq!(ts.n_nodes,9);
        assert_eq!(ts.rs_i,rs_i);
        assert_eq!(ts.rs_o,rs_o);
        assert_eq!(ts.area,4.0);



    }
    #[test] 
    fn test_calc_heat_flow_with_mass(){

        
        // Materials

        let polyurethane = Substance::new(
            "polyurethane".to_string(),
            0.0252, // W/m.K            
            2400., // J/kg.K
            17.5, // kg/m3... reverse engineered from paper            
        ); 

        let m1 = Material::new(Rc::clone(&polyurethane),200./1000.);

        let c = Construction::new("wall 1".to_string(), vec![Rc::clone(&m1)]);

        


        // FIRST TEST -- WITH MASS

        // Geometry
        let mut the_loop = Loop3D::new();
        let l = 1. as f64;
        the_loop.push( Point3D::new(-l, -l, 0.)).unwrap();
        the_loop.push( Point3D::new(l, -l, 0.)).unwrap();
        the_loop.push( Point3D::new(l, l, 0.)).unwrap();
        the_loop.push( Point3D::new(-l, l, 0.)).unwrap();
        the_loop.close().unwrap();
        let p = Polygon3D::new(the_loop).unwrap();
        
        // Physical surface
        let surface = Surface::new(p,Rc::clone(&c));        

        let main_dt = 300.0;
        let max_dx = m1.thickness()/4.0;
        let min_dt = 1.0;
        let (dt,nodes)=find_dt_and_n_nodes(&c, main_dt, 1, max_dx, min_dt);
        let ts = ThermalSurface::new(Rc::clone(&surface),dt,&nodes,0);
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
        // Materials

        let polyurethane = Substance::new(
            "polyurethane".to_string(),
            0.0252, // W/m.K            
            2400., // J/kg.K
            17.5, // kg/m3... reverse engineered from paper            
        );

        let m1 = Material::new(Rc::clone(&polyurethane),200./1000.);

        let c = Construction::new("wall 1".to_string(),vec![Rc::clone(&m1)]);

        // SECOND TEST -- NO MASS

        // Geometry 
        let mut the_loop = Loop3D::new();
        let l = 1. as f64;
        the_loop.push( Point3D::new(-l, -l, 0.)).unwrap();
        the_loop.push( Point3D::new(l, -l, 0.)).unwrap();
        the_loop.push( Point3D::new(l, l, 0.)).unwrap();
        the_loop.push( Point3D::new(-l, l, 0.)).unwrap();
        the_loop.close().unwrap();
        let p = Polygon3D::new(the_loop).unwrap();
        
        // Physical surface
        let surface = Surface::new(p,Rc::clone(&c));        

        let main_dt = 300.0;
        let max_dx = m1.thickness()/15.;
        let min_dt = 80.;
        let (dt,nodes)=find_dt_and_n_nodes(&c, main_dt, 1, max_dx, min_dt);
        let ts = ThermalSurface::new(Rc::clone(&surface),dt,&nodes,0);
        
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
        
        // Materials

        let polyurethane = Substance::new(
            "polyurethane".to_string(),
            0.0252, // W/m.K            
            2400., // J/kg.K
            17.5, // kg/m3... reverse engineered from paper            
        );

        let brickwork = Substance::new(
            "brickwork".to_string(),
            0.816, // W/m.K            
            800., // J/kg.K
            1700., // kg/m3... reverse engineered from paper            
        );

        // Geometry
        let mut the_loop = Loop3D::new();
        let l = 1. as f64;
        the_loop.push( Point3D::new(-l, -l, 0.)).unwrap();
        the_loop.push( Point3D::new(l, -l, 0.)).unwrap();
        the_loop.push( Point3D::new(l, l, 0.)).unwrap();
        the_loop.push( Point3D::new(-l, l, 0.)).unwrap();
        the_loop.close().unwrap();
        let p = Polygon3D::new(the_loop).unwrap();


        
        
        let m1 = Material::new(Rc::clone(&polyurethane),20./1000.);
        let m2 = Material::new(Rc::clone(&brickwork),220./1000.);
        
        let m3 = Rc::clone(&m1);
        
        let c = Construction::new("wall 3".to_string(),vec![Rc::clone(&m1),Rc::clone(&m2),Rc::clone(&m3)]);
        

        // Physical surface
        let surface = Surface::new(p,Rc::clone(&c));        
        
        let main_dt = 300.0;
        let max_dx = 0.015;
        let min_dt = 65.;
        let (dt,nodes)=find_dt_and_n_nodes(&c, main_dt, 1, max_dx, min_dt);
        let ts = ThermalSurface::new(Rc::clone(&surface),dt,&nodes,0);
        
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

        let brickwork = Substance::new(
            "brickwork".to_string(),
            0.816, // W/m.K            
            800., // J/kg.K
            1700., // kg/m3... reverse engineered from paper            
        );
        


        let m1 = Material::new(Rc::clone(&brickwork),200./1000.);
        let c = Construction::new("wall 1".to_string(),vec![Rc::clone(&m1)]);
        
        // Physical surface
        let surface = Surface::new(p,Rc::clone(&c));        


        // FIRST TEST -- 10 degrees on each side
        let main_dt = 300.0;
        let max_dx = m1.thickness()/2.0;
        let min_dt = 1.0;
        let (dt,nodes)=find_dt_and_n_nodes(&c, main_dt, 1, max_dx, min_dt);
        let mut ts = ThermalSurface::new(Rc::clone(&surface),dt,&nodes,0);
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

        let exp_q = (30.0 - 10.0)/(c.r_value()+ts.rs_i+ts.rs_o);
        assert!((exp_q-final_qin).abs() < 1E-4);
        assert!((exp_q+final_qout).abs() < 1E-4);
    }


    #[test]
    fn test_march_nomass(){
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

        let brickwork = Substance::new(
            "brickwork".to_string(),
            0.816, // W/m.K            
            800., // J/kg.K
            1700., // kg/m3... reverse engineered from paper            
        );


        let m1 = Material::new(Rc::clone(&brickwork), 3./1000.);
        let c = Construction::new("wall 1".to_string(),vec![Rc::clone(&m1)]);

        // Physical surface
        let surface = Surface::new(p,Rc::clone(&c));        


        // FIRST TEST -- 10 degrees on each side
        let main_dt = 300.0;
        let max_dx = m1.thickness()/2.0;
        let min_dt = 100.0;
        let (dt,nodes)=find_dt_and_n_nodes(&c, main_dt, 1, max_dx, min_dt);
        let mut ts = ThermalSurface::new(Rc::clone(&surface),dt,&nodes,0);
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

        let exp_q = (30.0 - 10.0)/(c.r_value()+ts.rs_i+ts.rs_o);
        assert!((exp_q-q_in).abs() < 1E-4);
        assert!((exp_q+q_out).abs() < 1E-4);
    }
}



    
   

