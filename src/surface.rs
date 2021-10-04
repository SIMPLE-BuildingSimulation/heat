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

use std::rc::Rc;
use matrix::Matrix;

use crate::Float;
use simple_model::construction::Construction;
use simple_model::boundary::Boundary;
use simple_model::fenestration::Fenestration;
use simple_model::surface::Surface;
use simple_model::simulation_state::{SimulationStateHeader, SimulationState};
use simple_model::simulation_state_element::SimulationStateElement;
use convection::*;

use crate::construction::*;

pub enum ThermalSurface{
    Fenestration(Rc<Fenestration>,ThermalSurfaceData),
    Surface(Rc<Surface>,ThermalSurfaceData),
    // Floor(ThermalSurfaceData)
}

impl ThermalSurface{
    pub fn new_surface(
        // model: &SimpleModel,
        state: &mut SimulationStateHeader,
        surface: &Rc<Surface>,
        dt: Float,
        n_elements: &[usize],
        // index: usize,
    ) -> Result<Self, String> {
        
        let ref_surface_index = *surface.index().unwrap();
        let construction = &surface.construction;
        let (rs_front, rs_back) = calc_convection_coefficients(&**surface);
        let area = surface.area();

        // surface,
        let (first_node,last_node,data)= ThermalSurfaceData::new(            
            state,
            construction,
            dt,
            area,
            rs_front,
            rs_back,
            n_elements,
            ref_surface_index, 
            false // not fenestration           
        )?;
        surface.set_first_node_temperature_index(first_node);
        surface.set_last_node_temperature_index(last_node);
        Ok(Self::Surface(Rc::clone(surface),data))
    }

    pub fn new_fenestration(
        // model: &SimpleModel,
        state: &mut SimulationStateHeader,
        fenestration: &Rc<Fenestration>,
        dt: Float,
        n_elements: &[usize],
        // index: usize,
    ) -> Result<Self, String> {
        
        let ref_surface_index = *fenestration.index().unwrap();        
        let construction = &fenestration.construction;
        let area = fenestration.area(); // should not fail because surface is full
        let (rs_front, rs_back) = calc_convection_coefficients_for_fenestration(fenestration);

        let (first_node,last_node,data) = ThermalSurfaceData::new(
            // model,
            state,
            construction,
            dt,
            area,
            rs_front,
            rs_back,
            n_elements,
            ref_surface_index,            
            true, // it is a fenestration
        )?;
        fenestration.set_first_node_temperature_index(first_node);
        fenestration.set_last_node_temperature_index(last_node);
        Ok(Self::Fenestration(Rc::clone(fenestration), data))
    }

    /// Borrows the [`ThermalSurfaceData`] of this surface
    fn data(&self)->&ThermalSurfaceData{
        match &self{
            Self::Surface(_,d)=>d,
            Self::Fenestration(_,d)=>d,
        }
    }

    /// Borrows a mutable version of [`ThermalSurfaceData`] of this surface
    fn mut_data(&mut self)->&mut ThermalSurfaceData{
        match self{
            Self::Surface(_,d)=>d,
            Self::Fenestration(_,d)=>d,
        }
    }

    /// Gets the `area` of the surface
    pub fn area(&self) -> Float {
        self.data().area
    }

    /// Gets the `rs_front`
    pub fn rs_front(&self) -> Float {
        self.data().rs_front
    }

    /// Gets the `rs_back`
    pub fn rs_back(&self) -> Float {
        self.data().rs_back
    }

    pub fn set_front_boundary(&mut self, b: Boundary) {
        self.mut_data().front_boundary = Some(b);
    }

    pub fn set_back_boundary(&mut self, b: Boundary) {
        self.mut_data().back_boundary = Some(b);
    }

    pub fn front_boundary(&self) -> &Option<Boundary> {
        &self.data().front_boundary
    }

    pub fn back_boundary(&self) -> &Option<Boundary> {
        &self.data().back_boundary
    }

    /// Checks whether a wall has thermal mass
    pub fn is_massive(&self) -> bool {
        self.data().massive
    }

    pub fn front_temperature(&self, state: &SimulationState) -> Float {
        match &self{
            Self::Surface(s,_)=>{
                s.first_node_temperature(state).expect("Trying to get front temperature from Surface without first_node state index")
            },
            Self::Fenestration(s,_) =>{
                s.first_node_temperature(state).expect("Trying to get front Fenestration from Surface without first_node state index")
            }
        }
    }

    pub fn back_temperature(&self, state: &SimulationState) -> Float {
        match &self{
            Self::Surface(s,_)=>{
                s.last_node_temperature(state).expect("Trying to get back temperature from Surface without last_node state index")
            },
            Self::Fenestration(s,_) =>{
                s.last_node_temperature(state).expect("Trying to get back temperature from Fenestration without last_node state index")
            }
        }
    }

    pub fn set_node_temperatures(
        &self,        
        state: &mut SimulationState,
        matrix: &Matrix,
    ) {
        let (ini, fin) = match &self{
            Self::Surface(s,d)=>{
                let ini = s.first_node_temperature_index().expect("Surface to SET node temperatures to does not have a First Node state index");
                let fin = ini + d.n_nodes;
                (ini,fin)
            },
            Self::Fenestration(s,d) =>{
                let ini = s.first_node_temperature_index().expect("Fenestration to SET node temperatures to does not have a First Node state index");
                let fin = ini + d.n_nodes;
                (ini,fin)
            }
        };
        for (node_index, i) in (ini..fin).enumerate(){
            let new_t = matrix.get(node_index, 0).unwrap();
            state[i] = new_t;

        }
    }
    
    /// Retrieves the state of the Surface as a Matrix
    /// object.
    pub fn get_node_temperatures(
        &self,
        state: &SimulationState,        
    ) -> Matrix {
        let (ini, fin) = match &self{
            Self::Surface(s,d)=>{
                let ini = s.first_node_temperature_index().expect("Surface to GET node temperatures from does not have a First Node state index");
                let fin = ini + d.n_nodes;
                (ini,fin)
            },
            Self::Fenestration(s,d) =>{
                let ini = s.first_node_temperature_index().expect("Fenestration to GET node temperatures from does not have a First Node state index");
                let fin = ini + d.n_nodes;
                (ini,fin)
            }
        };
        let n_nodes = fin - ini;
        let mut ret = Matrix::new(0.0, n_nodes, 1);
        for (node_index,i) in (ini..fin).enumerate(){
            let temp = state[i];
            ret.set(node_index, 0, temp).unwrap();
        }
        ret
    }

    /// Calculates the heat flow out of the layer, based
    /// on the inside and outside temperatures
    fn calc_heat_flow(
        &self,        
        state: &SimulationState,
        t_front: Float,
        t_back: Float,
    ) -> (Float, Float) {
        // Positive is going out of the layer.
        let data = self.data();
        let q_in;
        let q_out;
        let t_si = self.front_temperature(state);
        let t_so = self.back_temperature( state);
        if data.massive {
            q_in = (t_si - t_front) / data.full_rs_front;
            q_out = (t_so - t_back) / data.full_rs_back;
        } else {
            q_in = (t_si - t_front) / data.rs_front;
            q_out = (t_so - t_back) / data.rs_back;
        }
        // return
        (q_in, q_out)
    }
    

    /// Marches one timestep. Returns front and back heat flow    
    pub fn march(
        &self,        
        state: &mut SimulationState,
        t_front: Float,
        t_back: Float,
        _dt: Float,
    ) -> (Float, Float) {
        let mut temperatures = self.get_node_temperatures( state);
        let data = self.data();

        // println!("T front = {} | T back = {} ", t_front, t_back);
        // println!(" OLD TEMPS = {}", temperatures);

        if data.massive {

            if let Some(func) = &data.kt4_func {

                // First
                let mut k1 = func(&temperatures, t_front, t_back);

                
                // returning "temperatures + k1" is Euler... continuing is 
                // Runge–Kutta 4th order
                
                // Second                
                let mut aux = k1.from_scale(0.5).unwrap(); //  aux = k1 /2
                aux.add_to_this(&temperatures).unwrap(); // aux = T + k1/2
                let mut k2 = func(&aux, t_front, t_back);

                // Third... put the result into `aux`
                k2.scale(0.5, &mut aux).unwrap(); //  aux = k2 /2                
                aux.add_to_this(&temperatures).unwrap(); // T + k2/2
                let mut k3 = func(&aux, t_front, t_back);

                // Fourth... put the result into `aux`
                k3.scale(1., &mut aux).unwrap(); //  aux = k3 
                aux.add_to_this(&temperatures).unwrap(); // aux = T + k3
                let mut k4 = func(&aux, t_front, t_back);
                
                // Scale them and add them all up 
                k1.scale_this(1./6.);
                k2.scale_this(1./3.);
                k3.scale_this(1./3.);
                k4.scale_this(1./6.);

                k1.add_to_this(&k2).unwrap();
                k1.add_to_this(&k3).unwrap();
                k1.add_to_this(&k4).unwrap();

                
                // Let's add it to the temperatures.
                temperatures.add_to_this(&k1).unwrap();
            }else{
                unreachable!()
            }            
            
        } else {
            // full_rs_front is, indeed, the whole R
            let q = (t_back - t_front) / data.total_r;

            let ts_front = t_front + q * data.rs_front;
            let ts_back = t_back - q * data.rs_back;

            temperatures.set(0, 0, ts_front).unwrap();
            temperatures.set(data.n_nodes - 1, 0, ts_back).unwrap();
        }
        // println!(" Temperatures_after = {}", temperatures);
        // Set state
        self.set_node_temperatures(state, &temperatures);

        // return
        self.calc_heat_flow( state, t_front, t_back)
    }
}

/// This is a Surface from the point of view of our thermal solver.
/// Since this module only calculate heat transfer (and not short-wave solar
/// radiation, e.g., light), both simple_model::Fenestration and simple_model::Surface
/// are treated in the same way.
pub struct ThermalSurfaceData {

    
    /// The front side convection coefficient
    rs_front: Float,

    /// The back side convection coefficient
    rs_back: Float,

    total_r : Float,

    kt4_func: Option<Box<dyn Fn(&Matrix, Float, Float)->Matrix>>,


    /// The interior (i.e. front side) resistance before
    /// any layer with mass. It includes the r_si and also
    /// any light-weight material at the front of the construction.
    /// If the first layer in the construction has mass, then this
    /// value will be equal to r_si
    full_rs_front: Float,

    /// The exterior (i.e. back side) resistance after
    /// the last layer with mass. It includes the r_so and also
    /// any light-weight material at the back of the contruction.
    /// If the first layer in the construction has mass, then
    /// this value will be equal to r_si
    full_rs_back: Float,

    // The location of the first temperature node
    // in the SimulationState
    // front_side_node_index : usize,
    /// The number of nodes after discretizing
    /// the construction
    n_nodes: usize,

    /// Has thermal mass at all?
    massive: bool,

    /// The location of the front boundary zone in the
    /// Zones array of the Thermal Model
    front_boundary: Option<Boundary>,

    /// The location of the back boundary zone in the
    /// Zones array of the Thermal Model    
    back_boundary: Option<Boundary>,

    /// The area of the Surface
    area: Float,
}

impl ThermalSurfaceData {

    /// Constructs a new ThermalSurface object.        
    fn new(        
        state: &mut SimulationStateHeader,        
        construction: &Rc<Construction>,
        dt: Float,
        area: Float,
        rs_front: Float,
        rs_back: Float,
        n_elements: &[usize],
        ref_surface_index: usize,
        is_fenestration: bool,
    ) -> Result<(usize, usize, Self), String> {
        
        // Calculate number of nodes in that construction.
        let n_nodes = calc_n_total_nodes(&n_elements)?;

        // Push elements to the SimulationStateHeader
        let mut first_node : Option<usize> = None;
        let mut last_node : Option<usize> = None;        
        for i in 0..n_nodes {
            let node = if is_fenestration {                
                state.push(SimulationStateElement::FenestrationNodeTemperature(ref_surface_index,i),22.0)
            } else {
                state.push(SimulationStateElement::SurfaceNodeTemperature(ref_surface_index,i),22.0,)
            };
            if i == 0 {
                first_node = Some(node);
            }else if i == n_nodes - 1 {
                last_node = Some(node);
            }
        }

        // Build ThermalSurface 
        let (first_massive, last_massive) = get_first_and_last_massive_elements(n_elements)?;
        let mut ret = ThermalSurfaceData {                        
            rs_front,
            rs_back,
            total_r : construction.r_value().unwrap() + rs_front + rs_back, 
            full_rs_front : calc_full_rs_front(construction, rs_front, first_massive),
            full_rs_back : calc_full_rs_back(construction, rs_back, last_massive),            
            n_nodes,
            massive: true,        
            front_boundary: None, // filled when setting boundary
            back_boundary: None,        
            area,   
            kt4_func: None,         
        };

        
        if first_massive == 0 && last_massive == 0 {
            // No-mass surface
            ret.massive = false;

        }else{
            // massive surface

            let func = build_thermal_network(                
                construction,
                first_massive,
                last_massive,
                dt,
                n_nodes,
                &n_elements,
                rs_front,
                rs_back,
                ret.full_rs_front,
                ret.full_rs_back,
            )
            .unwrap();
            ret.kt4_func = Some(Box::new(func));

        }
        // Build the thermal network for this surface
        
        
        // return
        debug_assert!(first_node.is_some());
        debug_assert!(last_node.is_some());
        Ok((first_node.unwrap(), last_node.unwrap(),ret))
    }

    

    
} // END OF SURFACE IMPL

/***********/
/* TESTING */
/***********/

#[cfg(test)]
mod testing {
    use super::*;
    use geometry3d::loop3d::Loop3D;
    use geometry3d::point3d::Point3D;
    use geometry3d::polygon3d::Polygon3D;
    use simple_model::substance::Substance;
    use simple_model::material::Material;
    use simple_model::construction::Construction;
    use simple_model::model::SimpleModel;

    fn add_polyurethane(model: &mut SimpleModel) -> Rc<Substance> {

        let mut poly = Substance::new("polyurethane".to_string());
        poly   .set_density(17.5)// kg/m3... reverse engineered from paper
                    .set_specific_heat_capacity(2400.)// J/kg.K
                    .set_thermal_conductivity(0.0252);// W/m.K

        let ret = model.add_substance(poly);                        
        assert_eq!(ret.thermal_diffusivity().unwrap(), 0.6E-6);        
        ret
    }

    fn add_brickwork(model: &mut SimpleModel) -> Rc<Substance> {
        let mut brickwork = Substance::new("brickwork".to_string());

        brickwork   .set_density(1700.)// kg/m3... reverse engineered from paper
                    .set_specific_heat_capacity(800.)// J/kg.K
                    .set_thermal_conductivity(0.816);// W/m.K

        
        let ret = model.add_substance(brickwork);                        
            
        assert!((ret.thermal_diffusivity().unwrap() - 0.6E-6).abs() < 0.00000001);
        
        ret
    }

    fn add_material(model: &mut SimpleModel, substance: &Rc<Substance>, thickness: Float) -> Rc<Material> {
        let mat = Material::new("123123".to_string(), Rc::clone(substance), thickness);
        
        model.add_material(mat)
    }

    

    #[test]
    fn test_new() {
        let mut model = SimpleModel::new("A model".to_string());

        /* SUBSTANCES */
        // Add polyurethane Substance
        let poly = add_polyurethane(&mut model);

        let m0_thickness = 200.0 / 1000. as Float;
        let m0 = add_material(&mut model, &poly, m0_thickness);

        /* CONSTRUCTION */
        let mut c0 = Construction::new("Wall".to_string());
        c0.layers.push(Rc::clone(&m0));
        let c0 = model.add_construction(c0);
        

        /* GEOMETRY */
        let mut the_loop = Loop3D::new();
        let l = 1. as Float;
        the_loop.push(Point3D::new(-l, -l, 0.)).unwrap();
        the_loop.push(Point3D::new(l, -l, 0.)).unwrap();
        the_loop.push(Point3D::new(l, l, 0.)).unwrap();
        the_loop.push(Point3D::new(-l, l, 0.)).unwrap();
        the_loop.close().unwrap();
        let p = Polygon3D::new(the_loop).unwrap();

        /* SURFACE */
        let s = Surface::new("Surface".to_string(), p, Rc::clone(&c0));        
        let surface = model.add_surface(s);

        // This should result in a dt of 300, with 8 layers of 0.25m thick.
        let main_dt = 300.0;
        let max_dx = m0.thickness / 4.0;
        let min_dt = 1.;
        let (n_subdivisions, nodes) =
            discretize_construction( &c0, main_dt, max_dx, min_dt);
        let dt = main_dt / n_subdivisions as Float;
        let mut state_header = SimulationStateHeader::new();
        let ts =
            ThermalSurface::new_surface( &mut state_header, &surface, dt, &nodes).unwrap();

        let (rs_front, rs_back) = calc_convection_coefficients(&surface);
        assert!(ts.data().massive);
        // assert_eq!(ts.data().n_nodes, 9);
        assert_eq!(ts.data().rs_front, rs_front);
        assert_eq!(ts.data().rs_back, rs_back);
        assert_eq!(ts.data().area, 4.0);
    }

    
    #[test]
    fn test_find_n_total_nodes() {
        let n_nodes = vec![2, 2];
        let (first, last) = get_first_and_last_massive_elements(&n_nodes).unwrap();
        assert_eq!(first, 0);
        assert_eq!(last, 2);
        assert_eq!(5, calc_n_total_nodes(&n_nodes).unwrap());

        let n_nodes = vec![1, 0, 0, 2];
        let (first, last) = get_first_and_last_massive_elements(&n_nodes).unwrap();
        assert_eq!(first, 0);
        assert_eq!(last, 4);
        assert_eq!(5, calc_n_total_nodes(&n_nodes).unwrap());

        let n_nodes = vec![1, 0, 2];
        let (first, last) = get_first_and_last_massive_elements(&n_nodes).unwrap();
        assert_eq!(first, 0);
        assert_eq!(last, 3);
        assert_eq!(5, calc_n_total_nodes(&n_nodes).unwrap());

        let n_nodes = vec![8];
        let (first, last) = get_first_and_last_massive_elements(&n_nodes).unwrap();
        assert_eq!(first, 0);
        assert_eq!(last, 1);
        assert_eq!(9, calc_n_total_nodes(&n_nodes).unwrap());

        let n_nodes = vec![0];
        let (first, last) = get_first_and_last_massive_elements(&n_nodes).unwrap();
        assert_eq!(first, 0);
        assert_eq!(last, 0);
        assert_eq!(2, calc_n_total_nodes(&n_nodes).unwrap());

        let n_nodes = vec![1, 0];
        let (first, last) = get_first_and_last_massive_elements(&n_nodes).unwrap();
        assert_eq!(first, 0);
        assert_eq!(last, 1);
        assert_eq!(2, calc_n_total_nodes(&n_nodes).unwrap());

        let n_nodes = vec![0, 1];
        let (first, last) = get_first_and_last_massive_elements(&n_nodes).unwrap();
        assert_eq!(first, 1);
        assert_eq!(last, 2);
        assert_eq!(2, calc_n_total_nodes(&n_nodes).unwrap());

        let n_nodes = vec![0, 1, 0];
        let (first, last) = get_first_and_last_massive_elements(&n_nodes).unwrap();
        assert_eq!(first, 1);
        assert_eq!(last, 2);
        assert_eq!(2, calc_n_total_nodes(&n_nodes).unwrap());

        let n_nodes = vec![0, 0, 0];
        let (first, last) = get_first_and_last_massive_elements(&n_nodes).unwrap();
        assert_eq!(first, 0);
        assert_eq!(last, 0);
        assert_eq!(2, calc_n_total_nodes(&n_nodes).unwrap());

        let n_nodes = vec![1, 0, 1];
        let (first, last) = get_first_and_last_massive_elements(&n_nodes).unwrap();
        assert_eq!(first, 0);
        assert_eq!(last, 3);
        assert_eq!(4, calc_n_total_nodes(&n_nodes).unwrap());
    }

    #[test]
    fn test_calc_heat_flow_with_mass() {
        let mut model = SimpleModel::new("A model".to_string());

        /* SUBSTANCES */

        let poly = add_polyurethane(&mut model);

        /* MATERIALS */

        let m0_thickness = 200.0 / 1000. as Float;
        let m0 = add_material(&mut model, &poly, m0_thickness);

        /* CONSTRUCTION */
        let mut c0 = Construction::new("Wall".to_string());
        c0.layers.push(Rc::clone(&m0));
        let c0 = model.add_construction(c0);

        /* GEOMETRY */
        let mut the_loop = Loop3D::new();
        let l = 1. as Float;
        the_loop.push(Point3D::new(-l, -l, 0.)).unwrap();
        the_loop.push(Point3D::new(l, -l, 0.)).unwrap();
        the_loop.push(Point3D::new(l, l, 0.)).unwrap();
        the_loop.push(Point3D::new(-l, l, 0.)).unwrap();
        the_loop.close().unwrap();
        let p = Polygon3D::new(the_loop).unwrap();

        /* SURFACE */
        let surface = Surface::new("Surface".to_string(), p, Rc::clone(&c0));
        let surface = model.add_surface(surface);
        

        /* TEST */        

        let main_dt = 300.0;
        let max_dx = m0.thickness / 4.0;
        let min_dt = 1.0;
        let (n_subdivisions, nodes) =
            discretize_construction( &c0, main_dt, max_dx, min_dt);

        let dt = main_dt / n_subdivisions as Float;
        let mut state_header = SimulationStateHeader::new();
        let ts =
            ThermalSurface::new_surface( &mut state_header, &surface, dt, &nodes).unwrap();
        assert!(ts.data().massive);

        let state = state_header.take_values().unwrap();
        

        // TEST
        let temperatures = ts.get_node_temperatures(&state);

        assert_eq!(22.0, temperatures.get(0, 0).unwrap());
        assert_eq!(22.0, temperatures.get(ts.data().n_nodes - 1, 0).unwrap());

        let t_front = 10.0;
        let q_in_ref = (22.0 - t_front) / ts.data().rs_front;

        let t_back = 10.0;
        let q_out_ref = (22.0 - t_back) / ts.data().rs_back;
        let (q_in, q_out) = ts.calc_heat_flow( &state, t_front, t_back);

        assert_eq!(q_in, q_in_ref);
        assert_eq!(q_out, q_out_ref);
    }

    #[test]
    fn test_calc_heat_flow_no_mass() {
        let mut model = SimpleModel::new("A model".to_string());

        /* SUBSTANCES */

        let poly = add_polyurethane(&mut model);

        /* MATERIALS */
        let m0_thickness = 200.0 / 1000.;
        let m0 = add_material(&mut model, &poly, m0_thickness);

        /* CONSTRUCTION */
        let mut c0 = Construction::new("Wall".to_string());
        c0.layers.push(Rc::clone(&m0));
        let c0 = model.add_construction(c0);

        /* GEOMETRY */
        let mut the_loop = Loop3D::new();
        let l = 1. as Float;
        the_loop.push(Point3D::new(-l, -l, 0.)).unwrap();
        the_loop.push(Point3D::new(l, -l, 0.)).unwrap();
        the_loop.push(Point3D::new(l, l, 0.)).unwrap();
        the_loop.push(Point3D::new(-l, l, 0.)).unwrap();
        the_loop.close().unwrap();
        let p = Polygon3D::new(the_loop).unwrap();

        /* SURFACE */
        let surface = Surface::new("Surface".to_string(), p, Rc::clone(&c0));
        let surface = model.add_surface(surface);

        
        let main_dt = 300.0;
        let max_dx = m0.thickness / 15.;
        let min_dt = 80.;
        let (n_subdivisions, nodes) = discretize_construction( &c0, main_dt, max_dx, min_dt);

        let dt = main_dt / n_subdivisions as Float;
        let mut state_header = SimulationStateHeader::new();
        let ts = ThermalSurface::new_surface( &mut state_header, &surface, dt, &nodes).unwrap();

        let state = state_header.take_values().unwrap();
        

        // TEST

        assert!(!ts.data().massive);
        let temperatures = ts.get_node_temperatures(&state);
        assert_eq!(22.0, temperatures.get(0, 0).unwrap());
        assert_eq!(22.0, temperatures.get(ts.data().n_nodes - 1, 0).unwrap());

        let t_front = 10.0;
        let q_in_ref = (22.0 - t_front) / ts.data().rs_front;

        let t_back = 10.0;
        let q_out_ref = (22.0 - t_back) / ts.data().rs_back;
        let (q_in, q_out) = ts.calc_heat_flow( &state, t_front, t_back);

        assert_eq!(q_in, q_in_ref);
        assert_eq!(q_out, q_out_ref);
    }

    #[test]
    fn test_calc_heat_flow_mixed_mass() {
        // THIRD TEST -- WITH ONE NO-MASS LAYER IN THE EXTERIOR AND ONE IN THE INTERIOR

        let mut model = SimpleModel::new("Some model".to_string());

        /* SUBSTANCES */
        let poly = add_polyurethane(&mut model);
        let brickwork = add_brickwork(&mut model);

        /* MATERIALS */
        let m1 = add_material(&mut model, &poly, 20. / 1000.);
        let m2 = add_material(&mut model, &brickwork, 220. / 1000.);

        /* CONSTRUCTION */
        let mut c = Construction::new("construction".to_string());
        c.layers.push(Rc::clone(&m1));
        c.layers.push(Rc::clone(&m2));
        c.layers.push(Rc::clone(&m1));
        let c = model.add_construction(c);

        /* GEOMETRY */

        let mut the_loop = Loop3D::new();
        let l = 1. as Float;
        the_loop.push(Point3D::new(-l, -l, 0.)).unwrap();
        the_loop.push(Point3D::new(l, -l, 0.)).unwrap();
        the_loop.push(Point3D::new(l, l, 0.)).unwrap();
        the_loop.push(Point3D::new(-l, l, 0.)).unwrap();
        the_loop.close().unwrap();
        let p = Polygon3D::new(the_loop).unwrap();

        /* SURFACE */
        let s = Surface::new("Surface 1".to_string(), p, Rc::clone(&c));
        let surface = model.add_surface(s);
        

        /* TEST */        
        let main_dt = 300.0;
        let max_dx = 0.015;
        let min_dt = 65.;
        let (n_subdivisions, nodes) = discretize_construction( &c, main_dt, max_dx, min_dt);

        let dt = main_dt / n_subdivisions as Float;
        let mut state_header = SimulationStateHeader::new();
        let ts = ThermalSurface::new_surface( &mut state_header, &surface, dt, &nodes).unwrap();

        let  state = state_header.take_values().unwrap();
        

        // TEST

        assert!(ts.data().massive);
        let temperatures = ts.get_node_temperatures(&state);
        assert_eq!(22.0, temperatures.get(0, 0).unwrap());
        assert_eq!(22.0, temperatures.get(ts.data().n_nodes - 1, 0).unwrap());

        let t_front = 10.0;
        let q_in_ref = (22.0 - t_front) / ts.data().full_rs_front;

        let t_back = 10.0;
        let q_out_ref = (22.0 - t_back) / ts.data().full_rs_back;
        let (q_in, q_out) = ts.calc_heat_flow( &state, t_front, t_back);

        assert_eq!(q_in, q_in_ref);
        assert_eq!(q_out, q_out_ref);
    }

    #[test]
    fn test_march_massive() {
        let mut model = SimpleModel::new("Some model".to_string());

        /* SUBSTANCES */
        let brickwork = add_brickwork(&mut model);

        /* MATERIALS */
        let m1 = add_material(&mut model, &brickwork, 20. / 1000.);

        /* CONSTRUCTION */
        let mut c = Construction::new("construction".to_string());
        c.layers.push(Rc::clone(&m1));
        let c = model.add_construction(c);
        
        /* GEOMETRY */
        let mut the_loop = Loop3D::new();
        let l = 1. as Float;
        the_loop.push(Point3D::new(-l, -l, 0.)).unwrap();
        the_loop.push(Point3D::new(l, -l, 0.)).unwrap();
        the_loop.push(Point3D::new(l, l, 0.)).unwrap();
        the_loop.push(Point3D::new(-l, l, 0.)).unwrap();
        the_loop.close().unwrap();
        let p = Polygon3D::new(the_loop).unwrap();

        /* SURFACE */
        let s = Surface::new("Surface 1".to_string(), p, Rc::clone(&c));
        let surface = model.add_surface(s);

        // FIRST TEST -- 10 degrees on each side
        let main_dt = 300.0;
        let max_dx = m1.thickness / 2.0;
        let min_dt = 1.0;
        let (n_subdivisions, nodes) = discretize_construction( &c, main_dt, max_dx, min_dt);
        let dt = main_dt / n_subdivisions as Float;

        let mut state_header = SimulationStateHeader::new();
        let ts = ThermalSurface::new_surface( &mut state_header, &surface, dt, &nodes).unwrap();
        assert!(ts.data().massive);

        let mut state = state_header.take_values().unwrap();

        

        // TEST

        // Try marching until q_in and q_out are zero.
        let mut q: Float = 9999000009.0;        
        let mut counter: usize = 0;
        while q.abs() > 0.00015 {
            let (q_in, q_out) = ts.march( &mut state, 10.0, 10.0, dt);
            // the same amount of heat needs to leave in each direction
            // println!("q_in = {}, q_out = {} | diff = {}", q_in, q_out, (q_in - q_out).abs());
            // assert!((q_in - q_out).abs() < 1E-5);
            

            // q_front is positive
            assert!(q_in >= 0.);
            assert!(q_out >= 0.);

            // q_in needs to be getting smaller
            // assert!(q_in < q);            
            q = q_in;

            counter += 1;
            if counter > 9999999 {
                panic!("Exceded number of iterations... q.abs() = {}", q.abs())
            }
        }

        // all nodes should be at 10.0 now.
        let temperatures = ts.get_node_temperatures(&state);
        for i in 0..ts.data().n_nodes {
            let t = temperatures.get(i, 0).unwrap();
            assert!((t - 10.0).abs() < 0.00015);
        }

        // SECOND TEST -- 10 degrees in and 30 out.
        // We expect the final q to be (30-10)/R from
        // outside to inside. I.E. q_in = (30-10)/R,
        // q_out = -(30-10)/R

        // March until q converges
        let mut change: Float = 99.0;
        let mut counter: usize = 0;
        let mut previous_q: Float = -125.0;
        let mut final_qin: Float = -12312.;
        let mut final_qout: Float = 123123123.;
        while change.abs() > 1E-10 {
            let (q_in, q_out) = ts.march( &mut state, 10.0, 30.0, dt);
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
        assert!((final_qin + final_qout).abs() < 0.00033);

        let r = c.r_value().unwrap();

        let exp_q = (30.0 - 10.0) / (r + ts.data().rs_front + ts.data().rs_back);
        assert!((exp_q - final_qin).abs() < 0.00033);
        assert!((exp_q + final_qout).abs() < 0.00033);
    }

    #[test]
    fn test_march_nomass() {
        let mut model = SimpleModel::new("A model".to_string());

        /* SUBSTANCE */
        let polyurethane = add_polyurethane(&mut model);

        /* MATERIAL */
        let m1 = add_material(&mut model, &polyurethane, 3. / 1000.);

        /* CONSTRUCTION */
        let mut c = Construction::new("Construction".to_string());
        c.layers.push(Rc::clone(&m1));
        let c = model.add_construction(c);
        
        
        /* GEOMETRY */
        let mut the_loop = Loop3D::new();
        let l = 1. as Float;
        the_loop.push(Point3D::new(-l, -l, 0.)).unwrap();
        the_loop.push(Point3D::new(l, -l, 0.)).unwrap();
        the_loop.push(Point3D::new(l, l, 0.)).unwrap();
        the_loop.push(Point3D::new(-l, l, 0.)).unwrap();
        the_loop.close().unwrap();
        let p = Polygon3D::new(the_loop).unwrap();

        /* SURFACE */
        let s = Surface::new("WALL".to_string(), p, Rc::clone(&c));        
        let surface = model.add_surface(s);
        let mut state_header = SimulationStateHeader::new();

        /* TEST */
        

        // FIRST TEST -- 10 degrees on each side
        let main_dt = 300.0;
        let max_dx = m1.thickness / 2.0;
        let min_dt = 100.0;
        let (n_subdivisions, nodes) = discretize_construction( &c, main_dt, max_dx, min_dt);

        let dt = main_dt / n_subdivisions as Float;

        let ts = ThermalSurface::new_surface( &mut state_header, &surface, dt, &nodes).unwrap();
        assert!(!ts.data().massive);

        let mut state = state_header.take_values().unwrap();
        

        // TEST
        
        // Try marching until q_in and q_out are zero.

        let (q_in, q_out) = ts.march( &mut state, 10.0, 10.0, dt);

        // this should show instantaneous update. So,
        let temperatures = ts.get_node_temperatures( &state);
        assert_eq!(temperatures.get(0, 0).unwrap(), 10.0);
        assert_eq!(temperatures.get(1, 0).unwrap(), 10.0);
        assert_eq!(q_in, 0.0);
        assert_eq!(q_out, 0.0);

        // SECOND TEST -- 10 degrees in and 30 out.
        // We expect the final q to be (30-10)/R from
        // outside to inside. I.E. q_in = (30-10)/R,
        // q_out = -(30-10)/R

        let (q_in, q_out) = ts.march(&mut state, 10.0, 30.0, dt);

        // Expecting
        assert!(q_in > 0.0);
        assert!(q_out < 0.0);
        assert!((q_in + q_out).abs() < 1E-6);

        let exp_q = (30.0 - 10.0) / (c.r_value().unwrap() + ts.data().rs_front + ts.data().rs_back);
        assert!((exp_q - q_in).abs() < 1E-4);
        assert!((exp_q + q_out).abs() < 1E-4);
    }
}
