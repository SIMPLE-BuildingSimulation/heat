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
use crate::environment::Environment;
use crate::Float;
use matrix::Matrix;
use simple_model::{
    Boundary, Construction, Fenestration, SimulationStateHeader, Substance, Surface,
};
use simple_model::{SimulationState, SimulationStateElement};
use std::rc::Rc;

/// Marches forward through time, solving the
/// Ordinary Differential Equation that governs the heat transfer in walls.
///
///
/// The equation to solve is the following:
///
/// ```math
/// \overline{C}  \dot{T} - \overline{K}  T = q
/// ```
///
/// Where $`\overline{C}`$ and $`\overline{K}`$ are matrices representing the
/// thermal mass of each node and the thermal network, respectively; and where $`T`$ and $`q`$ are
/// vectors representing the Temperature and the heat "flow into" each node, respectively.
/// $`\overline{C}`$ and $`\overline{K}`$ are build based on the finite difference method.
///
/// This model uses a 4th order [Runga-Kutte](https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods) (a.k.a., RK4)
/// to march through time. In order to do this, it is convenient to write the equation to solve
/// as follows:
///
/// ```math
/// \dot{T}  = f(t, T)
/// ```
///
/// Where
/// ```math
/// f(t,T) = \overline{C}^{-1} \overline{K}  T + \overline{C}^{-1} q
/// ```
///
/// Then, the 4th order Runge-Kutta method allows marching forward through time as follows:
/// ```math
///  T_{i+1} = T_i + \frac{k_1 + 2k_2 + 2k_3 + k_4}{6}
/// ```
/// Where $`k_1`$, $`k_2`$, $`k_3`$ and $`k_4`$ can be calculated based on the
/// timestep $`\Delta t`$ as follows:
///
/// * $`k_1 = \Delta t \times f(t,T)`$
/// * $`k_2 = \Delta t \times f(t+\frac{\Delta t}{2}, T+\frac{k_1}{2})`$
/// * $`k_3 = \Delta t \times f(t+\frac{\Delta t}{2}, T+\frac{k_2}{2})`$
/// * $`k_4 = \Delta t \times f(t+\delta t, T+k_3 )`$
fn rk4(dt: Float, c: &Matrix, mut k: Matrix, mut q: Matrix, t: &mut Matrix) {
    let (krows, kcols) = k.size();
    assert_eq!(
        krows, kcols,
        "Expecting 'K' to be a squared matrix... nrows={}, ncols={}",
        krows, kcols
    );
    let (crows, ccols) = c.size();
    assert_eq!(
        crows, ccols,
        "Expecting 'C' to be a squared matrix... nrows={}, ncols={}",
        crows, ccols
    );
    let (qrows, qcols) = q.size();
    assert_eq!(
        qrows, krows,
        "Expecting 'q' to be to have {} rows because K has {} rows... found {}",
        krows, krows, qrows
    );
    assert_eq!(
        qcols, 1,
        "expecting 'q' to have 1 column... found {}",
        qcols
    );

    // Rearrenge into dT = (dt/C) * K + (dt/C)*q
    for nrow in 0..krows {
        let v = dt / c.get(nrow, nrow).unwrap();
        // transform k into k_prime (i.e., k * dt/C)
        let ini = if nrow == 0 { 0 } else { nrow - 1 };
        let fin = if nrow == krows - 1 {
            krows - 1
        } else {
            nrow + 1
        };
        // for ncol in 0..kcols{
        for ncol in ini..=fin {
            k.scale_element(nrow, ncol, v).unwrap();
        }
        // transfrom q into q_prime (i.e., q * dt/C)
        q.scale_element(nrow, 0, v).unwrap();
    }

    // get k1
    let mut k1 = &k * &*t;
    k1 += &q;
    
    // returning "temperatures + k1" is Euler... continuing is
    // Runge–Kutta 4th order
    // *t += &k1;
    // return;

    let mut aux = &k1 * 0.5;
    aux += t;

    // k2
    let mut k2 = &k * &aux;
    k2 += &q;

    // k3
    k2.scale_into(0.5, &mut aux).unwrap();
    aux += t;
    let mut k3 = &k * &aux;
    k3 += &q;

    // k4
    aux.copy_from(&k3);
    aux += t;
    let mut k4 = &k * &aux;
    k4 += &q;

    // Scale them and add them all up
    k1 /= 6.;
    k2 /= 3.;
    k3 /= 3.;
    k4 /= 6.;

    // Let's add it all and return
    *t += &k1;
    *t += &k2;
    *t += &k3;
    *t += &k4;
}

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

    fn first_node_temperature_index(&self) -> usize;

    fn last_node_temperature_index(&self) -> usize;

    fn front_temperature(&self, state: &SimulationState) -> Float {
        let i = self.first_node_temperature_index();
        state[i]
    }

    fn back_temperature(&self, state: &SimulationState) -> Float {
        let i = self.last_node_temperature_index();
        state[i - 1]
    }

    fn get_node_temperatures(&self, state: &SimulationState) -> Matrix {
        let ini = self.first_node_temperature_index();
        let fin = self.last_node_temperature_index()+1;
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

        for (node_index, i) in (ini..=fin).enumerate() {
            let new_t = matrix.get(node_index, 0).unwrap();
            state[i] = new_t;
        }
    }

    fn front_convection_coefficient(&self, state: &SimulationState) -> Option<Float>;
    fn back_convection_coefficient(&self, state: &SimulationState) -> Option<Float>;

    fn set_front_convection_coefficient(&self, state: &mut SimulationState, v: Float);

    fn set_back_convection_coefficient(&self, state: &mut SimulationState, v: Float);

    fn front_solar_irradiance(&self, state: &SimulationState) -> Float;
    fn back_solar_irradiance(&self, state: &SimulationState) -> Float;

    fn front_infrared_irradiance(&self, state: &SimulationState) -> Float;
    fn back_infrared_irradiance(&self, state: &SimulationState) -> Float;
}

impl SurfaceTrait for Surface {
    fn front_infrared_irradiance(&self, state: &SimulationState) -> Float {
        self.front_ir_irradiance(state).unwrap()
    }
    fn back_infrared_irradiance(&self, state: &SimulationState) -> Float {
        self.back_ir_irradiance(state).unwrap()
    }

    fn front_solar_irradiance(&self, state: &SimulationState) -> Float {
        self.front_incident_solar_irradiance(state).unwrap()
    }
    fn back_solar_irradiance(&self, state: &SimulationState) -> Float {
        self.back_incident_solar_irradiance(state).unwrap()
    }

    fn set_front_convection_coefficient(&self, _state: &mut SimulationState, _v: Float) {
        self.set_front_convection_coefficient(_state, _v)
    }
    fn set_back_convection_coefficient(&self, _state: &mut SimulationState, _v: Float) {
        self.set_back_convection_coefficient(_state, _v)
    }

    fn front_convection_coefficient(&self, _state: &SimulationState) -> Option<Float> {
        self.front_convection_coefficient(_state)
    }
    fn back_convection_coefficient(&self, _state: &SimulationState) -> Option<Float> {
        self.back_convection_coefficient(_state)
    }

    fn first_node_temperature_index(&self) -> usize {
        self.first_node_temperature_index()
            .expect("Could not get first node index in surface")
    }
    fn last_node_temperature_index(&self) -> usize {
        self.last_node_temperature_index()
            .expect("Could not get last node index in surface")
    }

    fn add_front_convection_state(
        &self,
        state: &mut SimulationStateHeader,
        ref_surface_index: usize,
    ) {
        let i = state.push(
            SimulationStateElement::SurfaceFrontConvectionCoefficient(ref_surface_index),
            10.,
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
            10.,
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
            0.0,
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
            0.0,
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
        self.set_last_node_temperature_index(last_node-1);
    }
}

impl SurfaceTrait for Fenestration {
    fn front_infrared_irradiance(&self, state: &SimulationState) -> Float {
        self.front_ir_irradiance(state).unwrap()
    }
    fn back_infrared_irradiance(&self, state: &SimulationState) -> Float {
        self.back_ir_irradiance(state).unwrap()
    }
    fn front_solar_irradiance(&self, state: &SimulationState) -> Float {
        self.front_incident_solar_irradiance(state).unwrap()
    }
    fn back_solar_irradiance(&self, state: &SimulationState) -> Float {
        self.back_incident_solar_irradiance(state).unwrap()
    }
    fn set_front_convection_coefficient(&self, state: &mut SimulationState, v: Float) {
        self.set_front_convection_coefficient(state, v)
    }
    fn set_back_convection_coefficient(&self, state: &mut SimulationState, v: Float) {
        self.set_back_convection_coefficient(state, v)
    }

    fn front_convection_coefficient(&self, _state: &SimulationState) -> Option<Float> {
        self.front_convection_coefficient(_state)
    }
    fn back_convection_coefficient(&self, state: &SimulationState) -> Option<Float> {
        self.back_convection_coefficient(state)
    }

    fn first_node_temperature_index(&self) -> usize {
        self.first_node_temperature_index()
            .expect("Could not get first node index in surface")
    }
    fn last_node_temperature_index(&self) -> usize {
        self.last_node_temperature_index()
            .expect("Could not get last node index in surface")
    }

    fn add_front_convection_state(
        &self,
        state: &mut SimulationStateHeader,
        ref_surface_index: usize,
    ) {
        let i = state.push(
            SimulationStateElement::FenestrationFrontConvectionCoefficient(ref_surface_index),
            10.,
        );
        self.set_front_convection_coefficient_index(i);
    }

    fn add_back_convection_state(
        &self,
        state: &mut SimulationStateHeader,
        ref_surface_index: usize,
    ) {
        let i = state.push(
            SimulationStateElement::FenestrationBackConvectionCoefficient(ref_surface_index),
            10.,
        );
        self.set_back_convection_coefficient_index(i);
    }

    fn add_front_convective_heatflow_state(
        &self,
        state: &mut SimulationStateHeader,
        ref_surface_index: usize,
    ) {
        let i = state.push(
            SimulationStateElement::FenestrationFrontConvectiveHeatFlow(ref_surface_index),
            0.0,
        );
        self.set_front_convective_heat_flow_index(i);
    }
    fn add_back_convective_heatflow_state(
        &self,
        state: &mut SimulationStateHeader,
        ref_surface_index: usize,
    ) {
        let i = state.push(
            SimulationStateElement::FenestrationBackConvectiveHeatFlow(ref_surface_index),
            0.0,
        );
        self.set_back_convective_heat_flow_index(i);
    }

    fn add_front_solar_irradiance_state(
        &self,
        state: &mut SimulationStateHeader,
        ref_surface_index: usize,
    ) {
        let i = state.push(
            SimulationStateElement::FenestrationFrontSolarIrradiance(ref_surface_index),
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
            SimulationStateElement::FenestrationBackSolarIrradiance(ref_surface_index),
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
            SimulationStateElement::FenestrationFrontIRIrradiance(ref_surface_index),
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
            SimulationStateElement::FenestrationBackIRIrradiance(ref_surface_index),
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
                SimulationStateElement::FenestrationNodeTemperature(ref_surface_index, node_index),
                22.0,
            );
        }
        let last_node = state.len();
        self.set_first_node_temperature_index(first_node);
        self.set_last_node_temperature_index(last_node-1);
    }
}

/// This is a Surface from the point of view of our thermal solver.
/// Since this module only calculate heat transfer (and not short-wave solar
/// radiation, e.g., light), both simple_model::Fenestration and simple_model::Surface
/// are treated in the same way.
pub struct ThermalSurfaceData<T: SurfaceTrait> {
    pub parent: Rc<T>,

    /// The [`Discretization`] that represents this `ThermalSurfaceData`
    pub discretization: Discretization,

    /// The location of the front boundary zone in the
    /// Zones array of the Thermal Model
    pub front_boundary: Option<Boundary>,

    /// The location of the back boundary zone in the
    /// Zones array of the Thermal Model    
    pub back_boundary: Option<Boundary>,

    /// The thermal absorbtance on the front side (from 0 to 1)
    pub front_emmisivity: Float,

    /// The thermal absorbtance on the back side (from 0 to 1)
    pub back_emmisivity: Float,

    /// The solar absorbtance on the front side (from 0 to 1)
    pub front_solar_absorbtance: Float,

    /// The solar absorbtance on the back side (from 0 to 1)
    pub back_solar_absorbtance: Float,

    /// The area of the Surface
    pub area: Float,
}

impl<T: SurfaceTrait> ThermalSurfaceData<T> {
    pub fn new(
        state: &mut SimulationStateHeader,
        ref_surface_index: usize,
        parent: &Rc<T>,
        area: Float,
        construction: &Rc<Construction>,
        discretization: Discretization,
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

        // let d = Discretization::new(construction, gases, model_dt, max_dx, min_dt, 1., 0.);

        // // Add node data.
        let n_nodes = discretization.segments.len();
        parent.add_node_temperature_states(state, ref_surface_index, n_nodes);

        const DEFAULT_SOLAR: Float = 0.7;
        let front_solar_absorbtance = match &construction.materials[0].substance {
            Substance::Normal(s) => match s.solar_absorbtance() {
                Ok(v) => *v,
                Err(_) => {
                    eprintln!(
                        "Substance '{}' has no solar absorbtance... assuming {}",
                        &construction.materials[0].substance.name(),
                        DEFAULT_SOLAR
                    );
                    DEFAULT_SOLAR
                }
            },
            _ => panic!("Front Emissivity not available for this particular kind of Substance"),
        };
        let back_solar_absorbtance = match &construction.materials.last().unwrap().substance {
            Substance::Normal(s) => match s.solar_absorbtance() {
                Ok(v) => *v,
                Err(_) => {
                    eprintln!(
                        "Substance '{}' has no solar absorbtance... assuming {}",
                        &construction.materials[0].substance.name(),
                        DEFAULT_SOLAR
                    );
                    DEFAULT_SOLAR
                }
            },
            _ => panic!("Front Emissivity not available for this particular kind of Substance"),
        };

        const DEFAULT_EM: Float = 0.84;
        let front_emmisivity = match &construction.materials[0].substance {
            Substance::Normal(s) => match s.thermal_absorbtance() {
                Ok(v) => *v,
                Err(_) => {
                    eprintln!(
                        "Substance '{}' has no thermal absorbtance... assuming {}",
                        &construction.materials[0].substance.name(),
                        DEFAULT_EM
                    );
                    DEFAULT_EM
                }
            },
            _ => panic!("Front Emissivity not available for this particular kind of Substance"),
        };
        let back_emmisivity = match &construction.materials.last().unwrap().substance {
            Substance::Normal(s) => match s.thermal_absorbtance() {
                Ok(v) => *v,
                Err(_) => {
                    eprintln!(
                        "Substance '{}' has no thermal absorbtance... assuming {}",
                        &construction.materials[0].substance.name(),
                        DEFAULT_EM
                    );
                    DEFAULT_EM
                }
            },
            _ => panic!("Front Emissivity not available for this particular kind of Substance"),
        };

        // Build resulting
        Ok(ThermalSurfaceData {
            parent: parent.clone(),
            area,
            discretization,
            front_boundary: None,
            back_boundary: None,
            front_emmisivity,
            back_emmisivity,
            front_solar_absorbtance,
            back_solar_absorbtance,
        })
    }

    pub fn set_front_boundary(&mut self, b: Boundary) {
        self.front_boundary = Some(b)
    }

    pub fn set_back_boundary(&mut self, b: Boundary) {
        self.back_boundary = Some(b)
    }

    /// Marches one timestep. Returns front and back heat flow    
    pub fn march(
        &self,
        state: &mut SimulationState,
        t_front: Float,
        t_back: Float,
        dt: Float,
    ) -> (Float, Float) {
        let mut temperatures = self.parent.get_node_temperatures(state);
        let (rows, ..) = temperatures.size();

        // Calculate and set Front and Back Solar Irradiance
        let solar_front = self.parent.front_solar_irradiance(state);
        let solar_back = self.parent.back_solar_irradiance(state);

        // Calculate and set Front and Back IR Irradiance
        // let (ir_front, ir_back) = (crate::SIGMA * ( t_front + 273.15 ).powi(4), crate::SIGMA * ( t_back + 273.15 ).powi(4));
        let ir_front = self.parent.front_infrared_irradiance(state);
        let ir_back = self.parent.back_infrared_irradiance(state);

        let front_env = Environment {
            air_temperature: t_front,
            ir_irrad: ir_front,
            solar_radiation: solar_front,
            ..Environment::default()
        };
        // // dbg!(front_env);
        let back_env = Environment {
            air_temperature: t_back,
            ir_irrad: ir_back,
            solar_radiation: solar_back,
            ..Environment::default()
        };

        // Calculate and set Front and Back convection coefficients
        let front_hs = front_env.get_hs();
        let back_hs = back_env.get_hs();
        self.parent
            .set_front_convection_coefficient(state, front_hs);
        self.parent.set_back_convection_coefficient(state, back_hs);

        /////////////////////
        // 1st: Calculate the solar absorption in each node
        /////////////////////
        let n_nodes = self.discretization.segments.len();
        dbg!("call glazing::alphas!");
        let mut q = Matrix::new(0.0, n_nodes, 1);
        q.add_to_element(0, 0, solar_front * self.front_solar_absorbtance)
            .unwrap();
        q.add_to_element(n_nodes - 1, 0, solar_back * self.back_solar_absorbtance)
            .unwrap();

        /////////////////////
        // 2nd: Calculate the temperature in all no-mass nodes.
        // Also, the heat flow into
        /////////////////////
        
        let (mass, nomass) = self.discretization.get_chunks();        

        for (ini, fin) in nomass {
            let mut count = 0;
            loop {
                let (k, mut local_q) = self.discretization.get_k_q(
                    ini,
                    fin,
                    &temperatures,
                    &front_env,
                    self.front_emmisivity,
                    front_hs,
                    &back_env,
                    self.back_emmisivity,
                    back_hs,
                );

                // ... here we can add solar gains
                for (local_i, i) in (ini..fin).into_iter().enumerate(){
                    let v = q.get(i, 0).unwrap();
                    local_q.add_to_element(local_i, 0, v).unwrap();
                }                
                local_q *= -1.; // transfrom equation from kT + q = 0 --> kT = -q

                let temps = k.mut_n_diag_gaussian(local_q, 3).unwrap(); // and just like that, q is the new temperatures                                

                let mut err = 0.0;
                for (local_i, i) in (ini..fin).into_iter().enumerate(){
                    let local_temp = temps.get(local_i, 0).unwrap();
                    let global_temp = temperatures.get(i, 0).unwrap();
                    err += (local_temp - global_temp).abs();                    
                }
                
                if err / (n_nodes as Float) < 0.01 {
                    break;
                }
                
                count += 1;
                assert!(count < 999, "Excessive number of iteration");
                for (local_i, i) in (ini..fin).into_iter().enumerate(){
                    let local_temp = temps.get(local_i, 0).unwrap();
                    temperatures.add_to_element(i, 0, local_temp).unwrap();
                    temperatures.scale_element(i, 0, 0.5).unwrap();
                           
                }
                
            }
        } 

        /////////////////////
        // 3rd: Calculate K and C matrices for the massive walls
        /////////////////////

        
        for (ini, fin) in mass{
            let c = self
                .discretization
                .segments
                .iter()
                .skip(ini)
                .take(fin-ini)
                .map(|(mass, _)| *mass)
                .collect();
            let c = Matrix::diag(c);

            let (k, mut local_q) = self.discretization.get_k_q(
                ini,
                fin,
                &temperatures,
                &front_env,
                self.front_emmisivity,
                front_hs,
                &back_env,
                self.back_emmisivity,
                back_hs,
            );

            // ... here we add solar gains
            for (local_i, global_i) in (ini..fin).into_iter().enumerate(){
                let v = q.get(global_i, 0).unwrap();
                local_q.add_to_element(local_i, 0, v).unwrap();
            }

            /////////////////////
            // 4rd: March massive nodes
            /////////////////////
            // Use RT4 for updating temperatures of massive nodes.
            let mut local_temps = Matrix::new(0.0, fin-ini, 1);
            for (local_i, global_i) in (ini..fin).into_iter().enumerate(){
                let v = temperatures.get(global_i, 0).unwrap();
                local_temps.set(local_i, 0, v).unwrap();
            }

            rk4(dt, &c, k, local_q, &mut local_temps);

            for (local_i, global_i) in (ini..fin).into_iter().enumerate(){
                let v = local_temps.get(local_i, 0).unwrap();
                temperatures.set(global_i, 0, v).unwrap();
            }
        } 

        /////////////////////
        // 5th: Set temperatures, calc heat-flows and return
        /////////////////////
        self.parent.set_node_temperatures(state, &temperatures);

        // Calc heat flow
        let ts_front = temperatures.get(0, 0).unwrap();
        let ts_back = temperatures.get(rows - 1, 0).unwrap();

        let flow_front = (ts_front - t_front) * front_hs;
        let flow_back = (ts_back - t_back) * back_hs;
        // dbg!(ts_front);
        // dbg!(ts_back);
        // dbg!(t_front);
        // dbg!(t_back);

        // Set state
        // self.parent.set_front_convective_heat_flow();

        (flow_front, flow_back)
    }
}

pub type ThermalSurface = ThermalSurfaceData<simple_model::Surface>;
pub type ThermalFenestration = ThermalSurfaceData<simple_model::Fenestration>;

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
            .set_thermal_absorbtance(0.)
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
            .set_thermal_absorbtance(0.)
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
        let d = Discretization::new(&c, main_dt, max_dx, min_dt, 1., 0.).unwrap();        
        let dt = main_dt / d.tstep_subdivision as Float;

        let mut state_header = SimulationStateHeader::new();
        let ts =
            ThermalSurface::new(&mut state_header, 0, &surface, surface.area(), &c, d).unwrap();

        let mut state = state_header.take_values().unwrap();

        // TEST

        // Try marching until q_in and q_out are zero.
        let mut q: Float = 9999000009.0;
        let mut counter: usize = 0;
        while q.abs() > 0.00015 {
            let (q_out, q_in) = ts.march(&mut state, 10.0, 10.0, dt);

            // the same amount of heat needs to leave in each direction
            // println!("q_in = {}, q_out = {} | diff = {}", q_in, q_out, (q_in - q_out).abs());
            assert!((q_in - q_out).abs() < 1E-5);

            // q_front is positive
            assert!(q_in >= 0., "q_in = {} | c = {}", q_in, counter);
            assert!(q_out >= 0., "q_out = {} | c = {}", q_out, counter);

            // q_in needs to be getting smaller            
            q = q_in;

            counter += 1;
            if counter > 9999999 {
                panic!("Exceded number of iterations... q.abs() = {}", q.abs())
            }
        }

        // all nodes should be at 10.0 now.
        let temperatures = ts.parent.get_node_temperatures(&state);
        let (n_nodes, ..) = temperatures.size();
        for i in 0..n_nodes {
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
        // assert!((final_qfront + final_qback).abs() < 0.00033, "final_qfront = {} | final_qback = {} | final_qfront + final_qback = {}", final_qfront, final_qback,final_qfront + final_qback);

        // let r = d.r_value();

        // let rs_front = 0.1; //surface.front_convection_coefficient(&state).unwrap();
        // let rs_back = 0.1; //surface.back_convection_coefficient(&state).unwrap();
        // let exp_q = (30.0 - 10.0) / (r + rs_front + rs_back);
        // assert!((exp_q - final_qfront).abs() < 0.00033, "In: exp: {} | final: {}", exp_q, final_qfront);
        // assert!((exp_q + final_qback).abs() < 0.00033, "Out: exp: {} | final: {}", exp_q, final_qback);
    }

    #[test]
    fn test_march_nomass() {
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
        let d = Discretization::new(&c, main_dt, max_dx, min_dt, 1., 0.).unwrap();
        let dt = main_dt / d.tstep_subdivision as Float;

        let ts =
            ThermalSurface::new(&mut state_header, 0, &surface, surface.area(), &c, d).unwrap();
        // assert!(!d.is_massive);

        let mut state = state_header.take_values().unwrap();

        // FIRST TEST -- 10 degrees on each side

        // Try marching until q_in and q_out are zero.

        let (q_in, q_out) = ts.march(&mut state, 10.0, 10.0, dt);

        // this should show instantaneous update. So,
        let temperatures = ts.parent.get_node_temperatures(&state);
        let (nnodes, ..) = temperatures.size();
        println!(" T == {}", &temperatures);
        assert!(
            (temperatures.get(0, 0).unwrap() - 10.0).abs() < 0.2,
            "T = {}",
            temperatures.get(0, 0).unwrap()
        );
        assert!(
            (temperatures.get(nnodes - 1, 0).unwrap() - 10.0).abs() < 0.2,
            "T = {}",
            temperatures.get(nnodes - 1, 0).unwrap()
        );
        assert!(q_in.abs() < 0.07, "q_in is {}", q_in);
        assert!(q_out.abs() < 0.07);
    }

    #[test]
    fn test_march_nomass_2() {
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
        let d = Discretization::new(&c, main_dt, max_dx, min_dt, 1., 0.).unwrap();
        let dt = main_dt / d.tstep_subdivision as Float;

        let ts =
            ThermalSurface::new(&mut state_header, 0, &surface, surface.area(), &c, d).unwrap();
        // assert!(!d.is_massive);

        let mut state = state_header.take_values().unwrap();

        // SECOND TEST -- 10 degrees in and 30 out.
        // We expect the final q to be (30-10)/R from
        // outside to inside. I.E. q_in = (30-10)/R,
        // q_out = -(30-10)/R

        let t_front = 10.0;
        let t_back = 30.0;
        let (q_front, q_back) = ts.march(&mut state, t_front, t_back, dt);

        // Expecting
        let temperatures = ts.parent.get_node_temperatures(&state);

        println!(" T == {}", &temperatures);

        assert!(q_front > 0.0, "q_in = {}", q_front);
        assert!(q_back < 0.0, "q_out = {}", q_back);

        let full_qfront = q_front;
        let full_qback = q_back;

        assert!(
            (full_qfront + full_qback).abs() < 0.08,
            "q_front = {} | q_back = {} | delta = {}",
            full_qfront,
            full_qback,
            (full_qfront + full_qback).abs()
        );

        // let rs_front = 0.1; //surface.front_convection_coefficient(&state).unwrap();
        // let rs_back = 0.1;// surface.back_convection_coefficient(&state).unwrap();
        // // dbg!(rs_front);
        // // dbg!(rs_back);
        // let r = d.r_value();
        // let poly_r = 0.006/0.0252;
        // assert!((r-poly_r).abs() < 1e-15);
        // let exp_q = (30.0 - 10.0) / (r + rs_front + rs_back);
        // assert!((exp_q - q_front).abs() < 1E-4, "exp_qin = {} ... found {}", exp_q, q_front);
        // assert!((exp_q + q_back).abs() < 1E-4, "exp_qout = {} ... found {}", exp_q, q_back);
    }

    #[test]
    fn test_rk4() {
        // Setup
        let c = Matrix::from_data(2, 2, vec![1., 0., 0., 1.]);
        let k = Matrix::from_data(2, 2, vec![1., -3., 4., -6.]);
        let q = Matrix::from_data(2, 1, vec![0., 0.]);

        // This system has a transient solution equals to c1*[0.75 1]'* e^(-3t) + c2*[1 1]' * e^(-2t)...

        // Find initial conditions so that C1 and C2 are 1.
        let mut temperatures = Matrix::from_data(2, 1, vec![0.75 + 1., 2.]);

        let temp_a_fn = |time: Float| 0.75 * (-3. * time).exp() + (-2. * time).exp();
        let temp_b_fn = |time: Float| (-3. * time).exp() + (-2. * time).exp();

        let mut time = 0.0;
        loop {
            let temp_a = temperatures.get(0, 0).unwrap();
            let exp_temp_a = temp_a_fn(time);
            let diff_a = (temp_a - exp_temp_a).abs();

            let temp_b = temperatures.get(1, 0).unwrap();
            let exp_temp_b = temp_b_fn(time);
            let diff_b = (temp_b - exp_temp_b).abs();
            const SMOL: Float = 1e-8;
            assert!(
                diff_a < SMOL,
                "temp_a = {} | exp_temp_a = {}, diff = {}",
                temp_a,
                exp_temp_a,
                diff_a
            );
            assert!(
                diff_b < SMOL,
                "temp_b = {} | exp_temp_b = {}, diff = {}",
                temp_b,
                exp_temp_b,
                diff_b
            );

            rk4(0.01, &c, k.clone(), q.clone(), &mut temperatures);

            time += 0.01;

            if time > 100. {
                break;
            }
        }
    }

    
}
