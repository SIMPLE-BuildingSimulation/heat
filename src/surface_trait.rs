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

use crate::construction::Discretization;
use crate::gas::Gas;
use crate::Float;
use matrix::Matrix;
use simple_model::{Boundary, Construction, SimulationStateHeader, Surface};
use simple_model::{SimulationState, SimulationStateElement};
use std::rc::Rc;

pub trait SurfaceTrait {
    fn add_front_convection_state(
        &self,
        state: &mut SimulationStateHeader,
        ref_surface_index: usize,
    );
    fn add_back_convection_state(
        &self,
        state: &mut SimulationStateHeader,
        ref_surface_index: usize,
    );

    fn add_front_convective_heatflow_state(
        &self,
        state: &mut SimulationStateHeader,
        ref_surface_index: usize,
    );
    fn add_back_convective_heatflow_state(
        &self,
        state: &mut SimulationStateHeader,
        ref_surface_index: usize,
    );

    fn add_front_solar_irradiance_state(
        &self,
        state: &mut SimulationStateHeader,
        ref_surface_index: usize,
    );
    fn add_back_solar_irradiance_state(
        &self,
        state: &mut SimulationStateHeader,
        ref_surface_index: usize,
    );

    fn add_front_ir_irradiance_state(
        &self,
        state: &mut SimulationStateHeader,
        ref_surface_index: usize,
    );
    fn add_back_ir_irradiance_state(
        &self,
        state: &mut SimulationStateHeader,
        ref_surface_index: usize,
    );

    fn add_node_temperature_states(
        &self,
        state: &mut SimulationStateHeader,
        ref_surface_index: usize,
        n_nodes: usize,
    );

    fn first_node_temperature_index(&self)->usize;
    fn last_node_temperature_index(&self)->usize;

    fn get_node_temperatures(&self, state: &SimulationState)->Matrix{
        let ini = self.first_node_temperature_index();
        let fin = self.last_node_temperature_index();
        let n_nodes = fin - ini;
        let mut ret = Matrix::new(0.0, n_nodes, 1);
        for (node_index, i) in (ini..fin).enumerate() {
            let temp = state[i];
            ret.set(node_index, 0, temp).unwrap();
        }
        ret

    }

    fn set_node_temperatures(&self, state: &mut SimulationState, matrix: &Matrix) {
        let ini = self.first_node_temperature_index();
        let fin = self.last_node_temperature_index();
        
        
        for (node_index, i) in (ini..fin).enumerate() {
            let new_t = matrix.get(node_index, 0).unwrap();
            state[i] = new_t;
        }
    }
}

impl SurfaceTrait for Surface {

    fn first_node_temperature_index(&self)->usize{
        self.first_node_temperature_index().expect("Could not get first node index in surface")
    }
    fn last_node_temperature_index(&self)->usize{
        self.last_node_temperature_index().expect("Could not get last node index in surface")
    }

    fn add_front_convection_state(
        &self,
        state: &mut SimulationStateHeader,
        ref_surface_index: usize,
    ) {
        let i = state.push(
            SimulationStateElement::SurfaceFrontConvectionCoefficient(ref_surface_index),
            0.1,
        );
        self.set_front_convection_coefficient_index(i);
    }

    fn add_back_convection_state(
        &self,
        state: &mut SimulationStateHeader,
        ref_surface_index: usize,
    ) {
        let i = state.push(
            SimulationStateElement::SurfaceBackConvectionCoefficient(ref_surface_index),
            0.1,
        );
        self.set_back_convection_coefficient_index(i);
    }

    fn add_front_convective_heatflow_state(
        &self,
        state: &mut SimulationStateHeader,
        ref_surface_index: usize,
    ) {
        let i = state.push(
            SimulationStateElement::SurfaceFrontConvectiveHeatFlow(ref_surface_index),
            0.1,
        );
        self.set_front_convective_heat_flow_index(i);
    }
    fn add_back_convective_heatflow_state(
        &self,
        state: &mut SimulationStateHeader,
        ref_surface_index: usize,
    ) {
        let i = state.push(
            SimulationStateElement::SurfaceBackConvectiveHeatFlow(ref_surface_index),
            0.1,
        );
        self.set_back_convective_heat_flow_index(i);
    }

    fn add_front_solar_irradiance_state(
        &self,
        state: &mut SimulationStateHeader,
        ref_surface_index: usize,
    ) {
        let i = state.push(
            SimulationStateElement::SurfaceFrontSolarIrradiance(ref_surface_index),
            0.0,
        );
        self.set_front_incident_solar_irradiance_index(i);
    }
    fn add_back_solar_irradiance_state(
        &self,
        state: &mut SimulationStateHeader,
        ref_surface_index: usize,
    ) {
        let i = state.push(
            SimulationStateElement::SurfaceBackSolarIrradiance(ref_surface_index),
            0.0,
        );
        self.set_back_incident_solar_irradiance_index(i);
    }

    fn add_front_ir_irradiance_state(
        &self,
        state: &mut SimulationStateHeader,
        ref_surface_index: usize,
    ) {
        let i = state.push(
            SimulationStateElement::SurfaceFrontIRIrradiance(ref_surface_index),
            0.0,
        );
        self.set_front_ir_irradiance_index(i);
    }
    fn add_back_ir_irradiance_state(
        &self,
        state: &mut SimulationStateHeader,
        ref_surface_index: usize,
    ) {
        let i = state.push(
            SimulationStateElement::SurfaceBackIRIrradiance(ref_surface_index),
            0.0,
        );
        self.set_back_ir_irradiance_index(i);
    }

    fn add_node_temperature_states(
        &self,
        state: &mut SimulationStateHeader,
        ref_surface_index: usize,
        n_nodes: usize,
    ) {
        // let n_nodes = d.segments.len();
        let first_node = state.len();
        for node_index in 0..n_nodes {
            state.push(
                SimulationStateElement::SurfaceNodeTemperature(ref_surface_index, node_index),
                22.0,
            );
        }
        let last_node = state.len();
        self.set_first_node_temperature_index(first_node);
        self.set_last_node_temperature_index(last_node);
    }
}



/// This is a Surface from the point of view of our thermal solver.
/// Since this module only calculate heat transfer (and not short-wave solar
/// radiation, e.g., light), both simple_model::Fenestration and simple_model::Surface
/// are treated in the same way.
pub struct ThermalSurfaceData<T: SurfaceTrait> 
{
    parent: Rc<T>,

    // / The front side convection coefficient
    // rs_front: Float,

    // // / The back side convection coefficient
    // // rs_back: Float,

    // /// Total resistance of the construction, WITHOUT Rsi and Rso
    // total_r: Float,

    /// The [`Discretization`] that represents this `ThermalSurfaceData`
    discretization: Rc<Discretization>,

    // /// The interior (i.e. front side) resistance before
    // /// any layer with mass. It includes any light-weight material at the front of the construction.
    // /// If the first layer in the construction has mass, then this
    // /// value will be equal to 0. This value DOES NOT include the r_si
    // pub r_front: Float,

    // /// The exterior (i.e. back side) resistance after
    // /// the last layer with mass. It includes
    // /// any light-weight material at the back of the contruction.
    // /// If the last layer in the construction has mass, then
    // /// this value will be equal to 0. This value DOES NOT include the r_so.
    // pub r_back: Float,

    // The location of the first temperature node
    // in the SimulationState
    // front_side_node_index : usize,
    /// The number of nodes after discretizing
    /// the construction
    n_nodes: usize,

    // /// Has thermal mass at all?
    // massive: bool,
    /// The location of the front boundary zone in the
    /// Zones array of the Thermal Model
    front_boundary: Option<Boundary>,

    /// The location of the back boundary zone in the
    /// Zones array of the Thermal Model    
    back_boundary: Option<Boundary>,

    /// The area of the Surface
    area: Float,
}

impl<T: SurfaceTrait> ThermalSurfaceData<T> 
{
    fn new(
        state: &mut SimulationStateHeader,
        ref_surface_index: usize,
        parent: &Rc<T>,
        area: Float,
        construction: &Rc<Construction>,
        model_dt: Float,
        // max_dx: Float,
        // min_dt: Float,
        // n_elements: &[usize],
        gases: &[Gas],
        discretization: &Rc<Discretization>,
    ) -> Result<ThermalSurfaceData<T>, String> {
        // Set Front and Back state
        parent.add_front_convection_state(state, ref_surface_index);
        parent.add_back_convection_state(state, ref_surface_index);

        parent.add_front_convective_heatflow_state(state, ref_surface_index);
        parent.add_back_convective_heatflow_state(state, ref_surface_index);

        parent.add_front_solar_irradiance_state(state, ref_surface_index);
        parent.add_back_solar_irradiance_state(state, ref_surface_index);

        parent.add_front_ir_irradiance_state(state, ref_surface_index);
        parent.add_back_ir_irradiance_state(state, ref_surface_index);

        // let d = Discretization::new(construction, gases, model_dt, max_dx, min_dt);

        // // Add node data.
        let n_nodes = discretization.segments.len();
        parent.add_node_temperature_states(state, ref_surface_index, n_nodes);

        // Build resulting
        Ok(ThermalSurfaceData {
            parent: parent.clone(),
            area,
            n_nodes,
            discretization: discretization.clone(),
            front_boundary: None,
            back_boundary: None,
        })
    }

    /// Marches one timestep. Returns front and back heat flow    
    pub fn march(
        &self,
        state: &mut SimulationState,
        t_front: Float,
        t_back: Float,
        _dt: Float,
    ) -> (Float, Float) {
        let mut temperatures = self.parent.get_node_temperatures(state);
        let (rows, _cols) = temperatures.size();

        // Calculate and set Front and Back convection coefficients
        let (rs_front, rs_back) = (0.1, 0.1);
        // match &self {
        //     Self::Fenestration(fen, _data) => {
        //         let (front, back) = calc_convection_coefficients_for_fenestration(&**fen);
        //         fen.set_front_convection_coefficient(state, front);
        //         fen.set_back_convection_coefficient(state, back);
        //         (front, back)
        //     }
        //     Self::Surface(sur, _data) => {
        //         let (front, back) = calc_convection_coefficients(&**sur);
        //         sur.set_front_convection_coefficient(state, front);
        //         sur.set_back_convection_coefficient(state, back);
        //         (front, back)
        //     }
        // };

        // Calculate and set Front and Back Solar Irradiance
        let (solar_front, solar_back) = (0., 0.);
        // match &self {
        //     Self::Fenestration(fen, _data) => {
        //         let front = fen
        //             .front_incident_solar_irradiance(state)
        //             .expect("Could not get front solar irradiance");
        //         let back = fen
        //             .back_incident_solar_irradiance(state)
        //             .expect("Could not get front solar irradiance");
        //         (front, back)
        //     }
        //     Self::Surface(sur, _data) => {
        //         let front = sur
        //             .front_incident_solar_irradiance(state)
        //             .expect("Could not get front solar irradiance");
        //         let back = sur
        //             .back_incident_solar_irradiance(state)
        //             .expect("Could not get back solar irradiance");
        //         (front, back)
        //     }
        // };

        // Calculate and set Front and Back IR Irradiance
        let (ir_front, ir_back) = (0., 0.);
        // match &self {
        //     Self::Fenestration(fen, _data) => {
        //         let front = fen
        //             .front_ir_irradiance(state)
        //             .expect("Could not get front IR irradiance");
        //         let back = fen
        //             .back_ir_irradiance(state)
        //             .expect("Could not get back IR irradiance");
        //         (front, back)
        //     }
        //     Self::Surface(sur, _data) => {
        //         let front = sur
        //             .front_ir_irradiance(state)
        //             .expect("Could not get front IR irradiance");
        //         let back = sur
        //             .back_ir_irradiance(state)
        //             .expect("Could not get back IR irradiance");
        //         (front, back)
        //     }
        // };

        // if data.massive {
        if let Some(func) = &self.discretization.march {
            func(&mut temperatures,
                t_front,
                t_back,
                rs_front,
                rs_back,
                solar_front,
                solar_back,
                ir_front,
                ir_back);
        } else {
            unreachable!()
        }
        // } else {
        //     // full_rs_front is, indeed, the whole R
        //     let q = (t_back - t_front) / (data.total_r + rs_front + rs_back);

        //     let ts_front = t_front + q * rs_front;
        //     let ts_back = t_back - q * rs_back;

        //     temperatures.set(0, 0, ts_front).unwrap();
        //     temperatures.set(data.n_nodes - 1, 0, ts_back).unwrap();
        // }

        dbg!("temperatures =");
        println!("{}", temperatures);
        // Set state
        self.parent.set_node_temperatures(state, &temperatures);

        // Calc heat flow
        let ts_front = temperatures.get(0, 0).unwrap();
        let ts_back = temperatures.get(rows-1, 0).unwrap();

        let flow_front = (ts_front - t_front) / rs_front;
        let flow_back = (ts_back - t_back) / rs_back;

        // Set state
        // self.parent.set_front_convective_heat_flow();


        (flow_front, flow_back)
    }
}

pub type ThermalSurface = ThermalSurfaceData<simple_model::Surface>;

/***********/
/* TESTING */
/***********/

#[cfg(test)]
mod testing {

    use super::*;
    use geometry3d::{Loop3D, Point3D, Polygon3D};

    use simple_model::{
        substance::Normal as NormalSubstance, Construction, Material, SimpleModel, Substance,
        Surface,
    };

    fn add_polyurethane(model: &mut SimpleModel) -> Substance {
        let mut poly = NormalSubstance::new("polyurethane".to_string());
        poly.set_density(17.5) // kg/m3... reverse engineered from paper
            .set_specific_heat_capacity(2400.) // J/kg.K
            .set_thermal_conductivity(0.0252); // W/m.K

        assert_eq!(poly.thermal_diffusivity().unwrap(), 0.6E-6);
        let ret = model.add_substance(poly.wrap());
        ret
    }

    fn add_brickwork(model: &mut SimpleModel) -> Substance {
        let mut brickwork = NormalSubstance::new("brickwork".to_string());

        brickwork
            .set_density(1700.) // kg/m3... reverse engineered from paper
            .set_specific_heat_capacity(800.) // J/kg.K
            .set_thermal_conductivity(0.816); // W/m.K

        assert!((brickwork.thermal_diffusivity().unwrap() - 0.6E-6).abs() < 0.00000001);
        let ret = model.add_substance(brickwork.wrap());

        ret
    }

    fn add_material(
        model: &mut SimpleModel,
        substance: Substance,
        thickness: Float,
    ) -> Rc<Material> {
        let mat = Material::new("123123".to_string(), substance.clone(), thickness);

        model.add_material(mat)
    }

    #[test]
    fn test_new_discretization() {
        let mut model = SimpleModel::new("A model".to_string());

        /* SUBSTANCES */
        // Add polyurethane Substance
        let poly = add_polyurethane(&mut model);

        let m0_thickness = 200.0 / 1000. as Float;
        let m0 = add_material(&mut model, poly, m0_thickness);

        /* CONSTRUCTION */
        let mut c0 = Construction::new("Wall".to_string());
        c0.materials.push(Rc::clone(&m0));
        let c0 = model.add_construction(c0);

        

        // This should result in a dt of 300, with 8 layers of 0.25m thick.
        let main_dt = 300.0;
        let max_dx = m0.thickness / 4.0;
        let min_dt = 1.;
        let gases = vec![];
        let  d = Discretization::new(&c0, &gases, main_dt, max_dx, min_dt); 
        assert!(d.is_massive);        
        assert!(d.is_static);
    }

    

   

    #[test]
    fn test_march_massive() {
        let mut model = SimpleModel::new("Some model".to_string());

        /* SUBSTANCES */
        let brickwork = add_brickwork(&mut model);

        /* MATERIALS */
        let m1 = add_material(&mut model, brickwork, 20. / 1000.);

        /* CONSTRUCTION */
        let mut c = Construction::new("construction".to_string());
        c.materials.push(Rc::clone(&m1));
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
        let gases = vec![];
        let mut d = Discretization::new(&c, &gases, main_dt, max_dx, min_dt); 
        let dt = main_dt / d.tstep_subdivision as Float;
        d.build_thermal_network(&gases, &c, dt, 0.0, 0.0);
        let d = Rc::new(d);
        let dt = main_dt / d.tstep_subdivision as Float;

        let mut state_header = SimulationStateHeader::new();
        let ts = ThermalSurface::new(&mut state_header, 0, &surface, surface.area(), &c, main_dt, &vec![], &d).unwrap();
        assert!(d.is_massive);

        let mut state = state_header.take_values().unwrap();

        // TEST

        // Try marching until q_in and q_out are zero.
        let mut q: Float = 9999000009.0;
        let mut counter: usize = 0;
        while q.abs() > 0.00015 {
            let (q_in, q_out) = ts.march(&mut state, 10.0, 10.0, dt);
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
        let temperatures = ts.parent.get_node_temperatures(&state);
        let (n_nodes, ..) = temperatures.size();
        for i in 0..n_nodes{
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
        let mut final_qfront: Float = -12312.;
        let mut final_qback: Float = 123123123.;
        while change.abs() > 1E-10 {
            let (q_front, q_back) = ts.march(&mut state, 10.0, 30.0, dt);
            final_qfront = q_front;
            final_qback = q_back;

            change = (q_front - previous_q).abs();
            previous_q = q_front;

            counter += 1;
            if counter > 99999 {
                panic!("Exceded number of iterations")
            }
        }

        // Expecting
        assert!(final_qfront > 0.0);
        assert!(final_qback < 0.0);
        assert!((final_qfront + final_qback).abs() < 0.00033);

        let n_nodes = d.segments.len();
        let r = d.r_between(0, n_nodes, &Vec::with_capacity(n_nodes), &vec![]);

        let rs_front = surface.front_convection_coefficient(&state).unwrap();
        let rs_back = surface.back_convection_coefficient(&state).unwrap();
        let exp_q = (30.0 - 10.0) / (r + rs_front + rs_back);
        assert!((exp_q - final_qfront).abs() < 0.00033, "In: exp: {} | final: {}", exp_q, final_qfront);
        assert!((exp_q + final_qback).abs() < 0.00033, "Out: exp: {} | final: {}", exp_q, final_qback);
    }

    #[test]
    fn test_march_nomass() {
        // https://github.com/LBNL-ETA/Windows-CalcEngine/blob/a80854f067e72d6b960330c163b5dbf5a56e66e1/src/Common/src/LinearSolver.cpp
        let mut model = SimpleModel::new("A model".to_string());

        /* SUBSTANCE */
        let polyurethane = add_polyurethane(&mut model);

        /* MATERIAL */
        let m1 = add_material(&mut model, polyurethane, 3. / 1000.);

        /* CONSTRUCTION */
        let mut c = Construction::new("Construction".to_string());
        c.materials.push(Rc::clone(&m1));
        c.materials.push(Rc::clone(&m1));
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

        let main_dt = 3.0;
        let max_dx = m1.thickness / 7.0;
        let min_dt = 10.0;
        let gases = vec![];
        let mut d = Discretization::new(&c, &gases, main_dt, max_dx, min_dt); 
        let dt = main_dt / d.tstep_subdivision as Float;
        d.build_thermal_network(&gases, &c, dt, 0.0, 0.0);
        let d = Rc::new(d);
        
        
        let ts = ThermalSurface::new(&mut state_header, 0, &surface, surface.area(), &c, main_dt, &vec![], &d).unwrap();
        // assert!(!d.is_massive);
        
        let mut state = state_header.take_values().unwrap();
        
        
        // FIRST TEST -- 10 degrees on each side

        // // Try marching until q_in and q_out are zero.

        // let (q_in, q_out) = ts.march(&mut state, 10.0, 10.0, dt);

        // // this should show instantaneous update. So,
        // let temperatures = ts.parent.get_node_temperatures(&state);
        // assert!( (temperatures.get(0, 0).unwrap()- 10.0).abs()<2E-3);
        // assert!( (temperatures.get(1, 0).unwrap()- 10.0).abs()<2E-3);
        // assert!( q_in.abs()<2E-2, "q_in is {}", q_in);
        // assert!( q_out.abs()<2E-2);

        // SECOND TEST -- 10 degrees in and 30 out.
        // We expect the final q to be (30-10)/R from
        // outside to inside. I.E. q_in = (30-10)/R,
        // q_out = -(30-10)/R

        let (q_in, q_out) = ts.march(&mut state, 10.0, 30.0, dt);

        // Expecting
        assert!(q_in > 0.0);
        assert!(q_out < 0.0);
        assert!((q_in + q_out).abs() < 0.05, "q_in = {} | q_out = {} | delta = {}", q_in, q_out, (q_in+q_out).abs());

        let rs_front = surface.front_convection_coefficient(&state).unwrap();
        let rs_back = surface.back_convection_coefficient(&state).unwrap();
        // let exp_q = (30.0 - 10.0) / (r_value(&c).unwrap() + rs_front + rs_back);
        // assert!((exp_q - q_in).abs() < 1E-4);
        // assert!((exp_q + q_out).abs() < 1E-4);
    }
}
