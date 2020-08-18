use building::construction::Construction;
use matrix::Matrix;
use building::material::Material;
use building::substance::Substance;

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

fn build_thermal_network(c: &Construction, dt: f64, n_nodes: &Vec<usize>, r_si: f64, r_so: f64, k_prime: & mut Matrix, k_f : & mut Matrix ) ->Result<(),String> {
    
    if n_nodes.len() != c.layers.len(){        
        let err = format!("Mismatch between number of layers in construction ({}) and node scheme ({})",c.layers.len(),n_nodes.len());
        return Err(err);
    }
    
    let all_nodes = calc_n_total_nodes(n_nodes);
    let (rows,cols) = k_prime.size();
    if rows != all_nodes || cols != all_nodes {
        let err = format!("Unexpected size of given matrix - found ({},{}) and was expecting ({},{})",rows,cols,all_nodes,all_nodes);
        return Err(err);
    }
    let (rows,cols) = k_f.size();
    if rows != all_nodes || cols != 1 {
        let err = format!("Unexpected size of given matrix - found ({},{}) and was expecting ({},{})",rows,cols,all_nodes,1);
        return Err(err);
    }
    
    // set matrices to Zero.
    for row in 0..rows{
        k_f.set(row,0,0.0).unwrap();
        for col in 0..cols{
            k_prime.set(row,col,0.0).unwrap();
        }
    }
    
    // NOW, PROCESS
    ////////////////
    let (first_massive, mut last_massive) = get_first_and_last_massive(n_nodes);
    
    if first_massive == 0 && last_massive == 0 {
        // no massive layers at all in construction.
        // Handle this case through the algorithm below.
        last_massive = c.layers.len();
    }

    // Calculate inner and outer surface Resistances
    let mut full_rsi = r_si;
    for i in 0..first_massive {
        let layer = c.layers[i];
        full_rsi += layer.thickness/layer.substance.thermal_conductivity;
    }
    let mut full_rso = r_so;
    for i in last_massive..c.layers.len() {        
        let layer = c.layers[i];
        full_rso += layer.thickness/layer.substance.thermal_conductivity;
    }
    
    // Calculate the rest.
    let mut node : usize = 0;
    let mut n_layer : usize = first_massive;

    while n_layer < last_massive {
        let material = c.layers[n_layer];    
        let m = n_nodes[n_layer];                
        if m == 0 { 
            
            // nomass material
            // add up all the R of the no-mass layers that
            // are together            
            let mut r = 0.0; // if the material is no mass, then the first value 
            while n_layer < last_massive && n_nodes[n_layer] == 0 {                
                let layer = c.layers[n_layer];
                let dx = layer.thickness;
                let k =layer.substance.thermal_conductivity;
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
            let k = material.substance.thermal_conductivity;                
            let dx = material.thickness/(m as f64);            
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
    if full_rsi > 0.0 {
        let old_value = k_prime.get(0,0).unwrap();        
        k_prime.set(0,0, old_value - 1.0/full_rsi ).unwrap();

        k_f.set(0,0,1.0/full_rsi).unwrap();
    }

    if full_rso > 0.0{        
        let old_value = k_prime.get(all_nodes-1,all_nodes-1).unwrap();
        k_prime.set(all_nodes-1,all_nodes-1, old_value - 1.0/full_rso ).unwrap();

        k_f.set(all_nodes-1, 0 ,1.0/full_rso).unwrap();        
    }
    
    

    // CALCULATE MASSES
    let mut left_side : Vec<f64> = vec![0.0;all_nodes];
    node = 0;

    for (n_layer,material) in c.layers.iter().enumerate(){        
        let m = n_nodes[n_layer];  

        if m != 0 {
            for _i in 0..m {                   
                // Calc mass
                let rho = material.substance.density;
                let cp = material.substance.heat_capacity;
                let dx = material.thickness/(m as f64);
                let m = rho * cp * dx / dt;                     
                
                left_side[node]+=m/2.0;
                left_side[node+1]+=m/2.0;
                //c.set(node,  0, c.get(node  ,0).unwrap()+m/2.).unwrap();
                //c.set(node+1,0, c.get(node+1,0).unwrap()+m/2.).unwrap();            
               
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
            let mass_value = 1./mass;        

            let old_value = k_f.get(i,0).unwrap();
            k_f.set(i, 0, old_value * mass_value).unwrap(); 

            let old_k_value = k_prime.get(i,i).unwrap();
            k_prime.set(i,i, old_k_value * mass_value).unwrap();                   
        }// Else, masses are already Zero.
        
        // This algorth implies that, if a construction
        // has no mass at all, a K_prime becomes a 2x2 matrix         
        // where k_prime[0,0]=-U, k_prime[0,1]=U, etc.
    }    
    
    return Ok(())
}

fn find_dt_and_n_nodes(c: &Construction, main_dt: f64, n: usize, max_dx: f64, min_dt: f64) -> (f64,Vec<usize>){

    let dt = main_dt/(n as f64);
    let safety = 1.5;

    // Choose a dx so that dt allows convergence.
    // stability is assured by (alpha * dt / dx^2 <= 1/2 )
    // meaning, we need to satisfy dx >= sqrt(2 * alpha * dt)
    
    // So, for each layer
    let mut n_nodes : Vec<usize> = vec![];

    for layer in c.layers.iter(){

        // Calculate the optimum_dx
        let alpha = layer.substance.thermal_diffusivity();
        let optimum_dx = (safety*2.0*alpha*dt).sqrt();
        let m = (layer.thickness/optimum_dx).floor();            
        let mut dx = layer.thickness/m;
        // dx cannot be smaller than thickness
        if m == 0. {
            dx = layer.thickness;
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

    #[test]
    fn test_get_dt_and_n_nodes(){
        
        
        let polyurethane = Substance {
            name: "polyurethane".to_string(),
            thermal_conductivity: 0.0252, // W/m.K            
            heat_capacity: 2400., // J/kg.K
            density: 17.5, // kg/m3... reverse engineered from paper            
        };
        assert_eq!(polyurethane.thermal_diffusivity(),0.6E-6);

        let brickwork = Substance {
            name: "brickwork".to_string(),
            thermal_conductivity: 0.816, // W/m.K            
            heat_capacity: 800., // J/kg.K
            density: 1700., // kg/m3... reverse engineered from paper            
        };
        assert_eq!(brickwork.thermal_diffusivity(),0.6E-6);

        /* WALL 1 */
        let m1 = Material {
            thickness: 200.0/1000.,
            substance: &polyurethane,
        };

        let c = Construction{
            name: "wall 1".to_string(),
            layers: vec![&m1],
        };

        // This should result in a dt of 300, with 8 layers of 0.25m thick.
        let main_dt = 300.;
        let (dt,n_nodes) = find_dt_and_n_nodes(&c,main_dt,1,m1.thickness/4.,1.);
        assert_eq!(n_nodes.len(),1);
        assert_eq!(dt,300.);
        assert_eq!(n_nodes[0],8);
        

        // This should result in a dt of 300, with 8 layers of 0.25m thick.
        let main_dt = 300.;
        let (dt,n_nodes) = find_dt_and_n_nodes(&c,main_dt,1,m1.thickness/8.,1.);
        assert_eq!(n_nodes.len(),1);
        assert_eq!(dt,300.);
        assert_eq!(n_nodes[0],8);
        

        // This should result in a dt of 300/2=150, with 12 layers of 0.0166667m thick.
        let main_dt = 300.;
        let (dt,n_nodes) = find_dt_and_n_nodes(&c,main_dt,1,m1.thickness/9.,1.);
        assert_eq!(n_nodes.len(),1);
        assert_eq!(dt,150.);
        assert_eq!(n_nodes[0],12);
        

        // This should result in a dt of 300/4=75, with 17 layers of 0.011...m thick.
        let main_dt = 300.;
        let (dt,n_nodes) = find_dt_and_n_nodes(&c,main_dt,1,m1.thickness/15.,1.);
        assert_eq!(n_nodes.len(),1);
        assert_eq!(dt,75.);
        assert_eq!(n_nodes[0],17);
        

        // We imposed a min_dt of 80, meaning that this should result in a no-mass layer
        let main_dt = 300.;
        let (dt,n_nodes) = find_dt_and_n_nodes(&c,main_dt,1,m1.thickness/15.,80.);
        assert_eq!(n_nodes.len(),1);
        assert_eq!(dt,100.);
        assert_eq!(n_nodes[0],0);
        

        /* WALL 2 */
        let m1 = Material {
            thickness: 200.0/1000.,
            substance: &polyurethane,
        };

        let m2 = Material{
            thickness: 110.0/1000.,
            substance: &brickwork,
        };

        let c = Construction{
            name: "wall 2".to_string(),
            layers: vec![&m1,&m2],
        };

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
        let m1 = Material {
            thickness: 20.0/1000.,
            substance: &polyurethane,
        };

        let m2 = Material{
            thickness: 220.0/1000.,
            substance: &brickwork,
        };

        let m3 = &m1;

        let c = Construction{
            name: "wall 3".to_string(),
            layers: vec![&m1,&m2,&m3],
        };

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
    fn test_build_thermal_network(){

        let polyurethane = Substance {
            name: "polyurethane".to_string(),
            thermal_conductivity: 0.0252, // W/m.K            
            heat_capacity: 2400., // J/kg.K
            density: 17.5, // kg/m3... reverse engineered from paper            
        };
        assert_eq!(polyurethane.thermal_diffusivity(),0.6E-6);

        let brickwork = Substance {
            name: "brickwork".to_string(),
            thermal_conductivity: 0.816, // W/m.K            
            heat_capacity: 800., // J/kg.K
            density: 1700., // kg/m3... reverse engineered from paper            
        };
        assert_eq!(brickwork.thermal_diffusivity(),0.6E-6);


        /* WALL 1: Single layer, with mass. */
        let m1 = Material {
            thickness: 200.0/1000.,
            substance: &polyurethane,
        };

        let c = Construction{
            name: "wall 1".to_string(),
            layers: vec![&m1],
        };

    
        let dt = 156.0;
        let n_nodes : Vec<usize> = vec![4];
        let all_nodes = calc_n_total_nodes(&n_nodes);

        let mut k_f = Matrix::new(0.0,all_nodes,1);
        let mut k_prime = Matrix::new(0.0,all_nodes,all_nodes);
        build_thermal_network(&c, dt, &n_nodes, 0.,0., &mut k_prime, &mut k_f).unwrap();

        let n_layer = 0;
        //let rho = m1.substance.density;
        //let cp = m1.substance.heat_capacity;
        let m = n_nodes[n_layer];
        let dx = m1.thickness/(m as f64);
        let k = m1.substance.thermal_conductivity;
        let alpha = m1.substance.thermal_diffusivity();

        // Check first node
        let node = 0;
        assert_eq!(k_f.get(node,0).unwrap(),0.0);// 2.0*dt/(rho*cp*dx));
        
        let found = k_prime.get(node,node).unwrap();
        let exp = -2.0*alpha*dt/dx/dx;
        assert!( (found-exp).abs() < 1E-10 );        

        let found = k_prime.get(node,node + 1).unwrap();
        let exp = k/dx;
        assert!( (found-exp).abs() < 1E-10 );

        // check middle nodes
        for node in 1..all_nodes-1{            
            assert_eq!(k_f.get(node,0).unwrap(), 0.0);//dt/(rho*cp*dx));

            let found = k_prime.get(node,node).unwrap();
            let exp = -2.0*alpha*dt/dx/dx;
            assert!( (found-exp).abs() < 1E-10 );

            let found = k_prime.get(node,node + 1).unwrap();
            let exp = k/dx;
            assert!( (found-exp).abs() < 1E-10 );

            let found = k_prime.get(node,node - 1).unwrap();
            let exp = k/dx;
            assert!( (found-exp).abs() < 1E-10 );
        }
        // check end node
        let node = all_nodes-1;        
        assert_eq!(k_f.get(node,0).unwrap(), 0.0);//2.0*dt/(rho*cp*dx));

        let found = k_prime.get(node,node).unwrap();
        let exp = -2.0*alpha*dt/dx/dx;
        assert!( (found-exp).abs() < 1E-10 );        

        let found = k_prime.get(node,node - 1).unwrap();
        let exp = k/dx;
        assert!( (found-exp).abs() < 1E-10 );


        /* WALL 2 : Two materials with mass */
        let m1 = Material {
            thickness: 200.0/1000.,
            substance: &polyurethane,
        };

        let m2 = Material{
            thickness: 110.0/1000.,
            substance: &brickwork,
        };

        let c = Construction{
            name: "wall 2".to_string(),
            layers: vec![&m1,&m2],
        };

        let dt = 156.0;
        let n_nodes : Vec<usize> = vec![3,3];
        let all_nodes = calc_n_total_nodes(&n_nodes);

        let mut k_f = Matrix::new(0.0,all_nodes,1);
        let mut k_prime = Matrix::new(0.0,all_nodes,all_nodes);
        build_thermal_network(&c, dt, &n_nodes, 0.,0., &mut k_prime, &mut k_f).unwrap();

        // FIRST LAYER
        //////////////
        let n_layer = 0;
        //let rho = m1.substance.density;
        //let cp = m1.substance.heat_capacity;
        let m = n_nodes[n_layer];
        let dx = m1.thickness/(m as f64);
        let k = m1.substance.thermal_conductivity;
        let alpha = m1.substance.thermal_diffusivity();

        // Check first node in layer
        let node = 0;
        assert_eq!(k_f.get(node,0).unwrap(), 0.0);//2.0*dt/(rho*cp*dx));
        
        let found = k_prime.get(node,node).unwrap();
        let exp = -2.0*alpha*dt/dx/dx;
        assert!( (found-exp).abs() < 1E-10 );        

        let found = k_prime.get(node,node + 1).unwrap();
        let exp = k/dx;
        assert!( (found-exp).abs() < 1E-10 );

        // check middle nodes
        for node in 1..3{            
            assert_eq!(k_f.get(node,0).unwrap(),0.0);// dt/(rho*cp*dx));

            let found = k_prime.get(node,node).unwrap();
            let exp = -2.0*alpha*dt/dx/dx;
            assert!( (found-exp).abs() < 1E-10 );

            let found = k_prime.get(node,node + 1).unwrap();
            let exp = k/dx;
            assert!( (found-exp).abs() < 1E-10 );

            let found = k_prime.get(node,node - 1).unwrap();
            let exp = k/dx;
            assert!( (found-exp).abs() < 1E-10 );
        }
        // check middle node (e.g. node 3), which 
        // is also the beginning of Layer 2
        let node = 3;
        let mass_1 = m1.substance.density * m1.substance.heat_capacity*m1.thickness/(6.0*dt);        
        let mass_2 = m2.substance.density * m2.substance.heat_capacity*m2.thickness/(6.0*dt);        

        assert_eq!(k_f.get(node,0).unwrap(), 0.0);//1.0/(mass_1+mass_2));
        let found = k_prime.get(node,node).unwrap();
        let u1 = 3.0*m1.substance.thermal_conductivity/m1.thickness;
        let u2 = 3.0*m2.substance.thermal_conductivity/m2.thickness;
        let exp = -(u1+u2)/(mass_1+mass_2);        
        assert_eq!(found,exp);
        
        let found = k_prime.get(node,node - 1).unwrap();        
        assert_eq!(found, u1);

        let found = k_prime.get(node,node + 1).unwrap();        
        assert_eq!(found, u2);

        // SECOND LAYER
        //////////////
        let n_layer = 1;
        //let rho = m2.substance.density;
        //let cp = m2.substance.heat_capacity;
        let m = n_nodes[n_layer];
        let dx = m2.thickness/(m as f64);
        let k = m2.substance.thermal_conductivity;
        let alpha = m2.substance.thermal_diffusivity();

        
        // check middle nodes in layer 2
        for node in 4..6{            
            assert_eq!(k_f.get(node,0).unwrap(),0.0);// dt/(rho*cp*dx));

            let found = k_prime.get(node,node).unwrap();
            let exp = -2.0*alpha*dt/dx/dx;
            assert!( (found-exp).abs() < 1E-10 );

            let found = k_prime.get(node,node + 1).unwrap();
            let exp = k/dx;
            assert!( (found-exp).abs() < 1E-10 );

            let found = k_prime.get(node,node - 1).unwrap();
            let exp = k/dx;
            assert!( (found-exp).abs() < 1E-10 );
        }
        // check end node (e.g. node 3)
        let node = 6;        
        let mass_2 = m2.substance.density * m2.substance.heat_capacity*m2.thickness/(6.0*dt);        

        assert_eq!(k_f.get(node,0).unwrap(),0.0);// 1.0/mass_2);
        let found = k_prime.get(node,node).unwrap();        
        let u2 = 3.0*m2.substance.thermal_conductivity/m2.thickness;
        let exp = -u2/(mass_2);        
        assert_eq!(found,exp);
        
        let found = k_prime.get(node,node - 1).unwrap();        
        assert_eq!(found, u2);

        /* WALL 3 : Single layer, no mass */
        let m1 = Material {
            thickness: 200.0/1000.,
            substance: &polyurethane,
        };

        let c = Construction{
            name: "wall 3".to_string(),
            layers: vec![&m1],
        };

        let dt = 156.0;
        let n_nodes : Vec<usize> = vec![0];
        let all_nodes = calc_n_total_nodes(&n_nodes);
        let mut k_f = Matrix::new(0.0,all_nodes,1);
        let mut k_prime = Matrix::new(0.0,all_nodes,all_nodes);
        build_thermal_network(&c, dt, &n_nodes, 0.,0., &mut k_prime, &mut k_f).unwrap();        
        
        
        for i in 0..all_nodes{
            let found = k_f.get(i,0).unwrap();
            let exp = 0.;
            assert_eq!(exp,found);
        }        
        
        for row in 0..all_nodes{
            for col in 0..all_nodes{
                let exp = m1.substance.thermal_conductivity/m1.thickness;            
                let found = k_prime.get(row,col).unwrap();
                if row == col{
                    assert_eq!(-exp,found);
                }else{
                    assert_eq!(exp,found);
                }
            }            
        }        

        /* WALL 4 : two equal layers, both with no mass */
        let m1 = Material {
            thickness: 200.0/1000.,
            substance: &polyurethane,
        };

        let c = Construction{
            name: "wall 4".to_string(),
            layers: vec![&m1, &m1],
        };

        let dt = 156.0;
        let n_nodes : Vec<usize> = vec![0,0];
        let all_nodes = calc_n_total_nodes(&n_nodes);
        let mut k_f = Matrix::new(0.0,all_nodes,1);
        let mut k_prime = Matrix::new(0.0,all_nodes,all_nodes);
        build_thermal_network(&c, dt, &n_nodes, 0.,0., &mut k_prime, &mut k_f).unwrap();        


        for i in 0..all_nodes{
            let found = k_f.get(i,0).unwrap();
            let exp = 0.;
            assert_eq!(exp,found);
        }        

        for row in 0..all_nodes{
            for col in 0..all_nodes{
                let exp = m1.substance.thermal_conductivity/(2.0*m1.thickness);            
                let found = k_prime.get(row,col).unwrap();
                if row == col{
                    assert_eq!(-exp,found);
                }else{
                    assert_eq!(exp,found);
                }
            }            
        }        

        /* WALL 5 : One layer with mass, one layer without it */
        
        let m1 = Material{
            thickness: 200.0/1000.,
            substance: &brickwork,
        };
        let m2 = Material {
            thickness: 200.0/1000.,
            substance: &polyurethane,
        };

        let c = Construction{
            name: "wall 5".to_string(),
            layers: vec![&m1, &m2],
        };

        let dt = 156.0;
        let n_nodes : Vec<usize> = vec![1,0];
        let all_nodes = calc_n_total_nodes(&n_nodes);
        assert_eq!(all_nodes,2);
        let mut k_f = Matrix::new(0.0,all_nodes,1);
        let mut k_prime = Matrix::new(0.0,all_nodes,all_nodes);
        build_thermal_network(&c, dt, &n_nodes, 0.,0., &mut k_prime, &mut k_f).unwrap();        


        let half_mass = m1.substance.density * m1.substance.heat_capacity * m1.thickness/(2.0*dt);
        let u = m1.substance.thermal_conductivity/m1.thickness;

        // Check first node
        let found = k_f.get(0,0).unwrap();
        let exp = 0.;
        assert_eq!(exp,found);

        let found = k_prime.get(0,0).unwrap();
        let exp = -u / half_mass;
        assert_eq!(exp,found);

        // Check last node
        let found = k_f.get(1,0).unwrap();
        let rso = m2.thickness/m2.substance.thermal_conductivity;
        let exp =  1.0/ half_mass / rso;
        assert_eq!(exp,found);
        
        let found = k_prime.get(1,1).unwrap();
        let exp = -(u + 1.0/rso) / half_mass;        
        assert!((exp-found).abs() < 1E-10);
         
        

        /* WALL 6 : One layer without mass, one layer with it */
        
        let m1 = Material{
            thickness: 200.0/1000.,
            substance: &brickwork,
        };
        let m2 = Material {
            thickness: 200.0/1000.,
            substance: &polyurethane,
        };

        let c = Construction{
            name: "wall 6".to_string(),
            layers: vec![&m2, &m1],
        };

        let dt = 156.0;
        let n_nodes : Vec<usize> = vec![0,1];
        let all_nodes = calc_n_total_nodes(&n_nodes);
        assert_eq!(all_nodes,2);
        let mut k_f = Matrix::new(0.0,all_nodes,1);
        let mut k_prime = Matrix::new(0.0,all_nodes,all_nodes);
        build_thermal_network(&c, dt, &n_nodes, 0.,0., &mut k_prime, &mut k_f).unwrap();        


        let half_mass = m1.substance.density * m1.substance.heat_capacity * m1.thickness/(2.0*dt);
        let u = m1.substance.thermal_conductivity/m1.thickness;

        // Check first node
        let found = k_f.get(0,0).unwrap();
        let rsi = m2.thickness/m2.substance.thermal_conductivity;
        let exp =  1.0/ half_mass / rsi;        
        assert_eq!(exp,found);

        let found = k_prime.get(0,0).unwrap();
        let exp = -(u + 1.0/rsi) / half_mass;        
        assert!((exp-found).abs() < 1E-10);

        // Check last node
        let found = k_f.get(1,0).unwrap();        
        let exp =  0.0;
        assert_eq!(exp,found);
        
        let found = k_prime.get(1,1).unwrap();
        let exp = -u / half_mass;        
        assert!((exp-found).abs() < 1E-10);


        /* WALL 7 : One layer with mass, two layers without it, and one with it again */
        
        let m1 = Material{
            thickness: 200.0/1000.,
            substance: &brickwork,
        };
        let m2 = Material {
            thickness: 200.0/1000.,
            substance: &polyurethane,
        };

        let c = Construction{
            name: "wall 7".to_string(),
            layers: vec![&m1, &m2, &m2, &m1],
        };

        let dt = 156.0;
        let n_nodes : Vec<usize> = vec![1,0,0,1];
        let all_nodes = calc_n_total_nodes(&n_nodes);
        assert_eq!(all_nodes,4);
        let mut k_f = Matrix::new(0.0,all_nodes,1);
        let mut k_prime = Matrix::new(0.0,all_nodes,all_nodes);
        build_thermal_network(&c, dt, &n_nodes, 0.,0., &mut k_prime, &mut k_f).unwrap();        

        

        let half_mass = m1.substance.density * m1.substance.heat_capacity * m1.thickness/(2.0*dt);
        let u = m1.substance.thermal_conductivity/m1.thickness;
        let insulation_r = m2.thickness * 2.0/m2.substance.thermal_conductivity;

        // Check first node
        let found = k_f.get(0,0).unwrap();        
        let exp =  0.0;
        assert_eq!(exp,found);

        let found = k_prime.get(0,0).unwrap();
        let exp = -u / half_mass;        
        assert!((exp-found).abs() < 1E-10);


        // Check second node
        let found = k_f.get(1,0).unwrap();        
        let exp =  0.0;
        assert_eq!(exp,found);

        let found = k_prime.get(1,1).unwrap();
        let exp = -(u + 1.0/insulation_r) / half_mass;        
        assert!((exp-found).abs() < 1E-10);
                
        // Check third node
        let found = k_f.get(2,0).unwrap();        
        let exp =  0.0;
        assert_eq!(exp,found);

        let found = k_prime.get(2,2).unwrap();
        let exp = -(u + 1.0/insulation_r) / half_mass;                
        assert_eq!(exp,found);

        // Check last node
        let found = k_f.get(3,0).unwrap();        
        let exp =  0.0;
        assert_eq!(exp,found);

        let found = k_prime.get(3,3).unwrap();
        let exp = -u / half_mass;        
        assert!((exp-found).abs() < 1E-10);
         
        

        
        
        
    
        
    }
        
}



    
   

