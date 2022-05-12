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

use crate::{Float, gas};
use matrix::Matrix;
use simple_model::{Construction, Substance};
use std::rc::Rc;
use crate::gas::Gas;


/// Represents a thermal resistance through a Wall.
#[derive(Debug, Clone, Copy)]
pub enum Resistance {
    /// A nornal (i.e., $`t/\lambda`$)resistance
    Solid(Float),

    /// A cavity, comprised of a gas
    Cavity(usize),

    /// The resistance is a surface coefficient.
    Back,

    /// Undefined yet
    None,
}

impl std::default::Default for Resistance {
    fn default() -> Self {
        Resistance::None
    }
}

/// This is the signature of a function that allows updating the temperatures of a 
/// `Discretization` by marching forward in time. The inputs are the following:
/// 
/// * `nodes_temps: &mut Matrix`,
/// * `air_temp_front: Float`,
/// * `air_temp_back: Float`,
/// * `rs_front: Float`,
/// * `rs_back: Float`,
/// * `solar_irradiance_front: Float`,
/// * `solar_irradiance_back: Float`,
/// * `ir_irradiance_front: Float`,
/// * `ir_irradiance_back: Float`
type ConstructionMarchFunction = dyn Fn(&mut Matrix, Float, Float, Float, Float, Float, Float, Float, Float);

/// Represents the discretization of a [`Construction`] for heat transfer
/// calculation purposes
pub struct Discretization {
    /// Contains the node mass and the `Resistance` of each segment
    pub segments: Vec<(Float, Resistance)>,

    /// Contains the minimum number of timesteps per model timestep that 
    /// this discretization requires to ensure numerical stability and accuracy.
    /// 
    /// While the caller model—e.g., a Multiphysics Simulation model—may attempt
    /// to simulate at a certain timestep, some of the constructions in it may 
    /// require smaller timesteps to ensure numerical stability and accuracy.
    /// This means that—on each timestep in the caller model—the thermal model needs
    /// to perform `time_subdivision` sub-timesteps. 
    pub tstep_subdivision: usize,

    /// Indicates whether the thermal network for this `Discretization` is static, 
    /// meaning that the Thermal Network will not change when marching through time.
    /// 
    /// The thermal network can also be dynamic; for example, when there are non-linearities
    /// in the system. An example of this are cavities—which contain Radiation heat transfer,
    /// meaning that the R-value of a cavity depends on temperature—and also
    /// materials that change their thermal properties with temperature (e.g., Phase Change
    /// Materials or materials that get wet).
    pub is_static: bool,

    /// Indicates if any of the nodes have any thermal mass
    pub is_massive: bool,

    /// The `ConstructionMarchFunction` that allows marching forward
    /// in time
    pub march: Option<Box<ConstructionMarchFunction>>,
}

impl Discretization {

    /// Creates a new `Discretization`.
    /// 
    /// It first calculates the `tstep_subdivision` and the number of elements 
    /// on each layer of Construction by calling `discretize_construction()`; and then builds the 
    /// `Discretization` by calling `build()`.
    pub fn new(
        construction: &Rc<Construction>,
        gases: &[Gas],
        model_dt: Float,
        max_dx: Float,
        min_dt: Float,
    )->Self{
        let (tstep_subdivision, n_elements) = Self::discretize_construction(construction, model_dt, max_dx, min_dt);        
        Self::build(construction, gases, tstep_subdivision, &n_elements)
    }

    /// Creates the `segments` of the `Discretization`.
    fn build(
        construction: &Rc<Construction>,
        gases: &[Gas],
        tstep_subdivision: usize,
        n_elements: &[usize],        
    )->Self{
        debug_assert_eq!(n_elements.len(), construction.materials.len());

        // Let's start with an empty set of segments
        let mut n_nodes : usize = n_elements.iter().sum();        
        n_nodes = n_nodes.max(construction.materials.len()); // At least one per layer... but Zero means  "no_mass wall"

        let mut segments : Vec<(Float, Resistance)> = vec![ (0.0, Resistance::default()); n_nodes+1];

        let mut n_segment = 0;
        for (n_layer,n) in n_elements.iter().enumerate(){
            let mut n = *n;
            let material = &construction.materials[n_layer];

            // get the mass of each segment.
            let mass = if n == 0 {
                0.0
            }else{
                match &material.substance {
                    Substance::Normal(s)=> {
                        let dx = material.thickness / n as Float; 
                        let rho = s.density().expect(
                            "Trying to calculate C_Matrix with a substance without 'density'",
                        );
                        let cp = s.specific_heat_capacity().expect("Trying to calculate C_Matrix with a substance without 'specific heat capacity'");
                        let m = rho * cp * dx; // dt;
                        m    
                    },
                    Substance::Gas(_s)=>{
                        0.
                    }
                }
            };

            if n == 0 {
                n = 1;
            }
            // Now iterate all segments... 
            for _ in 0..n {
                match &material.substance {
                    Substance::Normal(s)=> {
                        
                        // Add mass to this and next nodes (if it is NoMass, it is Zero)
                        segments[n_segment].0 += mass/2.;                        
                        segments[n_segment+1].0 += mass/2.;
                        
                        
                        // Add resistance
                        let dx = material.thickness; 
                        let k = s.thermal_conductivity().expect(&format!("Substance '{}' in material '{}' in Construction '{}' has no thermal conductivity, but we need it", s.name(), material.name(), construction.name()));

                        segments[n_segment].1 = Resistance::Solid(dx/k);

    
                    },
                    Substance::Gas(s)=>{

                        
                        // Search by name
                        let i = {
                            let mut index : Option<usize> = None;
                            for (i,g) in gases.iter().enumerate() {
                                if &g.name == s.name() {
                                    index = Some(i);
                                }
                            }
                            index.unwrap()                                       
                        };
                        segments[n_segment].1 = Resistance::Cavity(i);
                    }
                }
                n_segment+=1;
            }
            // Process last one
            segments[n_nodes].1 = Resistance::Back;
        }

        let is_static = segments.iter().enumerate().all(|(index, (_mass, resistance))| {
            // Last one does not count
            if index == n_nodes {
                true
            }else if let Resistance::Solid(_) = resistance {
                true
            }else{
                false
            }
        });
        let is_massive = segments.iter().any(|(mass, _resistance)| *mass > 0.0);

        Self {
            is_static,
            is_massive,
            segments,
            tstep_subdivision,
            march: None,// Added when build_thermal_network()
        }
    }

    
    /// Given a Maximum element thickness ($`\Delta x_{max}`$) and a minimum timestep ($`\Delta t_{min}`$), this function
    /// will find an arguibly good (i.e., stable and accurate) combination of $`\Delta t`$ and number of elements in each
    /// layer of the construction. 
    ///
    /// This function recursively increases the model's timestep subdivisions (`n`) in order to reduce $`\Delta t`$ to numbers
    /// that respect the restrictions of (1) stability, (2) $`\Delta x_{max}`$, and (3) $`\Delta t_{min}`$. In other words,
    /// it searches (by testing $`\Delta t_{model}/1`$, $`\Delta t_{model}/2`$, $`\Delta t_{model}/3`$, ... $`\Delta t_{model}/n`$)
    /// for the minimum `n` that respects this restrictions
    ///
    /// # The math behind it
    ///
    /// The first thing to know is that the walls in this module march
    /// through time using a 4th order [Runga-Kutte](https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods)
    /// (a.k.a., RK4). The second thing to know is that the RK4 method is
    /// more stable than the [Euler method](https://en.wikipedia.org/wiki/Euler_method),
    /// and thus the restrictions of stability for the Euler method can be considered
    /// to be a conservative restriction for the RK4. Hence, this function uses the
    /// Euler method restrictions.
    ///
    /// Now, as explained in the [`build_thermal_network`] documentation, we are solving
    /// the following equation:
    ///
    /// ```math
    /// \dot{T} = \overline{C}^{-1} \overline{K}  T + \overline{C}^{-1} q
    /// ```
    ///
    /// And thus the stability of the numerical method will depend on the matrix:
    ///
    /// ```math
    /// \overline{K}^{\star} =\Delta t \overline{C}^{-1} \overline{K}
    /// ```
    ///
    /// Specifically, we don't want any of its [eigenvalues](https://en.wikipedia.org/wiki/Eigenvalues_and_eigenvectors)
    /// $`\xi_1, \xi_2,\xi_3, ...`$ to be outside of the Euler method's stability region. Since this
    /// matrix has onle Real eigenvalues, this is equivalent to saying:
    ///
    /// ```math
    /// -2 < \xi_i < 0 ; \forall i
    /// ```
    ///
    /// However, finding the eigenvalues for $`\overline{K}^{\star}`$ is far from trivial. So there is
    /// yet another heuristic I am using: I am treating the case of a wall with 1 layer that is subdivided
    /// into a single element as the limit case. I am not sure if this is correct, but most of the instabilities
    /// I identified through Trial and Error corresponded to this case.
    ///
    /// For this limit case:
    /// * $`R = \frac{\Delta x}{\lambda}`$
    /// * $`C = \rho  c_p  \Delta x`$
    ///
    /// thus the value of $`\overline{K}^{\star}`$ is:
    ///
    /// ```math
    /// \overline{K}^{\star}=\begin{bmatrix}
    /// -\frac{\Delta t}{C\times R} - \frac{\Delta t}{C\times R_s} & \frac{\Delta t}{C\times R} \\
    ///  \frac{\Delta t}{C\times R} & -\frac{\Delta t}{C\times R} - \frac{\Delta t}{C\times R_s}\\
    /// \end{bmatrix}   
    ///```
    /// Note that, in that equation, $`R_{si} = R_{so}`$. The reason for this is that this method
    /// does not know the real values of $`R_{si}`$ and $`R_{so}`$, so it simply using a placeholder
    /// value low enough to cover most general cases (`0.05`).
    ///
    /// Then, it can be found that the eigenvaues of this case—which we are treating as the limit case—are:
    /// ```math
    /// \xi_1 = -\frac{\Delta t} { R_s \rho c_p \Delta x }
    /// ```
    /// ```math
    /// \xi_2 = \xi_1 - 2 \frac{\Delta t \lambda}{ \rho  c_p  {\Delta x}^2}
    /// ```
    /// Both these values are always negative, so we don't need to worry about
    /// $`\xi_i < 0 `$. Also, it can be noticed that $`\xi_2 < \xi_1`$, meaning that
    /// what we need to comply with for Euler's stability criteria is:
    /// ```math
    /// -\frac{\Delta t} { R_s \rho c_p \Delta x } - 2 \frac{\Delta t \lambda}{ \rho  c_p  {\Delta x}^2} < -2
    /// ```
    ///
    /// Which means that the chosen $`\Delta x`$ must be greater than the (apparently only) positive
    /// solution to equation:
    ///
    /// ```math
    /// 0 = 2 {\Delta x}^2 - \left( \frac{\Delta t}{\rho c_p R_s} \right) \Delta x - \frac{2 \Delta t \lambda}{\rho c_p}
    /// ```
    ///
    /// So, this method will identify a combination of $`\Delta t`$ and $`\Delta x`$ that
    /// allows complying with this
    ///
    /// All that said, the value for $`\Delta t`$ actually used by the final model is
    /// actually half of what these equations use. This is because what we are using
    /// is a heuristic and I want to be safe... ish
    fn discretize_construction(
        construction: &Rc<Construction>,
        model_dt: Float,
        max_dx: Float,
        min_dt: Float,
    ) -> (usize, Vec<usize>) {
        // I could only think of how to make this recursively... so I did this.
        fn aux(
            construction: &Rc<Construction>,
            main_dt: Float,
            n: usize,
            max_dx: Float,
            min_dt: Float,
        ) -> (usize, Vec<usize>) {
            let dt = main_dt / (n as Float);

            // So, for each layer
            let n_layers = construction.materials.len();
            let mut n_elements: Vec<usize> = Vec::with_capacity(n_layers);
            const RS: Float = 0.05;

            for n_layer in 0..n_layers {
                let material = &construction.materials[n_layer];
                let substance = &material.substance;

                // Calculate the minimum_dx
                let thickness = material.thickness;
                let (k, rho, cp) = match substance {
                    Substance::Normal(s) => {
                        let k = s.thermal_conductivity().expect("Trying to discretize a construction that contains a Normal Substance without a 'thermal conductivity'");
                        let rho = s.density().expect("Trying to discretize a construction that contains a Normal Substance without a 'density'");
                        let cp = s.specific_heat_capacity().expect("Trying to discretize a construction that contains a Normal Substance without a 'specific heat capacity'");
                        (*k, *rho, *cp)
                    }
                    Substance::Gas(_) => {
                        n_elements.push(0);
                        continue;
                    }
                };

                let a_coef = 2.;
                let b_coef = -dt / (rho * cp * RS);
                let c_coef = -2. * dt * k / (rho * cp);
                let disc = b_coef * b_coef - 4. * a_coef * c_coef;
                // this should never happen...?
                debug_assert!(disc >= 0.);

                // One solution is apparently always negative...
                // i.e. it is meaningless
                debug_assert!((-b_coef - disc.sqrt()) / (2. * a_coef) < 0.);

                // The positive solution is the one we care about
                let min_dx = (-b_coef + disc.sqrt()) / (2. * a_coef);

                if min_dx > thickness {
                    // This means that this layer cannot comply with the
                    // given timestep because its thickness leads to a dx that
                    // does not ensure convergence...
                    // check if there is room for reducing dt (hence reducing min_dx)
                    let next_dt = main_dt / ((n + 1) as Float);
                    if next_dt > min_dt {
                        // If there is room for that, do it.
                        return aux(construction, main_dt, n + 1, max_dx, min_dt);
                    } else {
                        // otherwise, mark this layer as no-mass
                        n_elements.push(0);
                    }
                } else {
                    // subdivide the layer, making all the elements of equal thickness
                    let m = (thickness / min_dx).floor();
                    // this case belongs to the other branch of this if/else
                    debug_assert!(m as usize != 0);
                    let dx = thickness / m;
                    if dx > max_dx {
                        // If the found dx is larger than the max allowed d_x, try to change timestep
                        // check if there is room for reducing dt...
                        let next_dt = main_dt / ((n + 1) as Float);
                        if next_dt > min_dt {
                            // If there is room for that, do it.
                            return aux(construction, main_dt, n + 1, max_dx, min_dt);
                        } else {
                            // otherwise, mark this layer as no-mass
                            n_elements.push(0);
                        }
                    } else {
                        // "dx" is smaller than max_dx, and thus this works
                        // fine.
                        n_elements.push(m as usize)
                    }
                }
            }

            // Check stability requirements...
            // stability is assured by (alpha * dt / dx^2 <= 1/2 )
            #[cfg(debug_assertions)]
            {
                for (n_layer, _) in n_elements.iter().enumerate() {
                    let material = &construction.materials[n_layer];
                    let substance = &material.substance;

                    // Calculate the optimum_dx
                    let thickness = material.thickness;
                    let (k, rho, cp) = match substance {
                        Substance::Normal(s) => {
                            let k = s.thermal_conductivity().unwrap();
                            let rho = s.density().unwrap();
                            let cp = s.specific_heat_capacity().unwrap();
                            (*k, *rho, *cp)
                        }
                        Substance::Gas(_) => continue,
                    };
                    let dt = main_dt / n as Float;
                    let dx = thickness / n_elements[n_layer] as Float;

                    // assert!(alpha * dt / dx / dx <= 0.5);
                    let lambda1 = -dt / (RS * rho * cp * dx);
                    let r = dx / k;
                    let lambda2 = lambda1 - 2. * dt / (r * rho * cp * dx);
                    assert!(lambda1 >= -2.);
                    assert!(lambda1 <= 0.);
                    assert!(lambda2 >= -2.);
                    assert!(lambda2 <= 0.);
                }
            }

            // return
            (n, n_elements)
        }
        aux(construction, model_dt, 1, max_dx, min_dt)
    }

    /// Gets a squared [`Matrix`] containing the mass of the massive nodes
    /// and ignoring the no-mass ones
    fn calc_c_matrix(
        &self,         
    ) -> Vec<Float> {
        let masses : Vec<Float> = self.segments.iter().filter(|(mass,_r)|{ *mass > 0.0}).map(|(mass, _r)| *mass).collect();
        masses
        
    }

    pub fn r_between(&self,ini: usize, fin: usize, _temps: &[Float], gases: &[Gas])->Float{
        

        self.segments[ini..fin].iter().map(|(_mass, resistance)|{
            match resistance {
                Resistance::Solid(r)=>*r,
                Resistance::Cavity(i)=>{
                    let gas = &gases[*i];
                    todo!();
                },
                Resistance::Back => 0.0,//panic!("Todo: Cover Back in r_between"),
                Resistance::None => panic!("Found Resistance::None when calculating r_between"),

            }
        }).sum()

        
    }

    /// Calculates the `K` matrix (i.e., the thermal network) for massive constructions
    ///
    /// Constructions are assumed to be a sandwich where zero or more massive
    /// layers are located between two non-mass layers. These non-mass layers will always
    /// include the interior and exterior film convections coefficients, respectively (which is
    /// why they are called `full_rs_front` and `full_rs_back`, respectively). Additionally,
    /// they can also include some lightweight insulation or any other material of negligible
    /// thermal mass.
    ///
    /// This matrix is constructed based ona discretization of the layers of the
    /// construction; that is, each layer is subdivided into `n` elements all of equal
    /// thickness. One node is placed at the beginning and end of each elements. Each element can
    /// be represented by the 2x2 matrix
    ///
    /// ```math
    /// \overline{K}=\begin{bmatrix} -1/R & 1/R \\
    /// 1/R & -1/R
    /// \end{bmatrix}   ;   R=\frac{thickness}{\lambda}
    ///```
    ///
    /// Hence—ignoring external inputs—the K matrix for a construction subdivided into 3 elements (i.e., 4 nodes) can be written as follows:
    /// ```math
    /// \overline{K}=\begin{bmatrix}
    /// -1/R_{1\rightarrow2} & 1/R_{1\rightarrow2} & 0 & 0 \\
    /// 1/R_{1\rightarrow2} & -1/R_{1\rightarrow2} - 1/R_{2\rightarrow3} & 1/R_{2\rightarrow3} & 0 \\
    /// 0 & 1/R_{2\rightarrow3} & -1/R_{2\rightarrow3} - 1/R_{3\rightarrow4} & 1/R_{3\rightarrow4} \\
    /// 0 & 0 & 1/R_{3\rightarrow4} & -1/R_{3\rightarrow4} \\
    /// \end{bmatrix}   
    ///```
    ///
    /// Now, these nodes are also connected to an interior and an exterior temperatures through all the
    /// layers that do not have any thermal mass both in the interior and exterior (i.e., $`R_{si,full}`$ and $`R_{so,full}`$, respectively).
    /// This means that the Matrix $`\overline{K}`$ needs to become:
    /// ```math
    /// \overline{K}=\begin{bmatrix}
    /// -1/R_{1\rightarrow2} - 1/R_{si,full} & 1/R_{1\rightarrow2} & 0 & 0 \\
    /// 1/R_{1\rightarrow2} & -1/R_{1\rightarrow2} - 1/R_{2\rightarrow3} & 1/R_{2\rightarrow3} & 0 \\
    /// 0 & 1/R_{2\rightarrow3} & -1/R_{2\rightarrow3} - 1/R_{3\rightarrow4} & 1/R_{3\rightarrow4} \\
    /// 0 & 0 & 1/R_{3\rightarrow4} & -1/R_{3\rightarrow4}- 1/R_{so,full} \\
    /// \end{bmatrix}   
    ///```
    ///
    /// > **NOTE:** This method returns such a matrix, without the $`R_{si, full}`$ and $`R_{so, full}`$. They need to
    /// be added when marching (because these values change over time).
    fn calc_k_matrix(
        &self,
        temps: &[Float],
        gases: &[Gas],
        // c: &Rc<Construction>,
        // first_massive: usize,
        // last_massive: usize,
        // n_elements: &[usize],
        // all_nodes: usize,
        massive_only: bool,
    ) -> Matrix {
        // We need onte temperature per node
        assert_eq!(temps.len(), self.segments.len());

        
        let massives : Vec<usize> = if massive_only {

            self.segments.iter().enumerate().filter_map(|(i,(mass, _r))| {
                if *mass > 0. {
                    Some(i)
                } else {
                    None
                }}).collect()
        }else{
            self.segments.iter().enumerate().map(|(i,..)| i).collect()
        };
        let n_massive = massives.len();
        // initialize k_prime
        let mut k_matrix = Matrix::new(0.0, n_massive, n_massive);
        
        // * massive_node_i counts each row/col of the K matrix.
        // * ini and fin hold the position of the massive node among
        //   all nodes, including non-massive ones
        dbg!("Account for whatever is BEFOR the first massive node");
        for massive_node_i in 0..massives.len() {
            if massive_node_i == massives.len()-1{
                // We have reached the end... so, don't bother
                break
            }
            let ini = massives[massive_node_i];
            let fin = massives[massive_node_i+1];
            let thermal_resistance = self.r_between(ini, fin, temps, gases);

            // update values
            let u_value = 1. / thermal_resistance;
            // top left
            let i = massive_node_i;
            let old_value = k_matrix.get(i, i).unwrap();
            k_matrix.set(i, i, old_value - u_value).unwrap();
            // top right
            let old_value = k_matrix.get(i, i + 1).unwrap();
            k_matrix.set(i, i + 1, old_value + u_value).unwrap();
            // bottom left
            let old_value = k_matrix.get(i + 1, i).unwrap();
            k_matrix.set(i + 1, i, old_value + u_value).unwrap();
            // bottom right
            let old_value = k_matrix.get(i + 1, i + 1).unwrap();
            k_matrix
                .set(i + 1, i + 1, old_value - u_value)
                .unwrap();
        }
        

        k_matrix
        
    }

    /// Builds the necessary data for marching forward through time, solving the
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
    /// Note that—unless some layer of the Surface is generating heat—all the elements of $`q`$ are Zero
    /// exept the first one and the last one, which are $`\frac{T_{in}}{R_{si}}`$ and $`\frac{T_{out}}{R_{so}}`$,
    /// respectively.
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
    ///
    /// So, what this method does is to calculate what is needed in order to return a closure that
    /// calculates $`\Delta t \times f(t,T)`$; that is to say:
    /// ```math
    /// return = \Delta t \times f(t,T) = \Delta t \times \overline{C}^{-1} \overline{K}  T + \Delta t \times \overline{C}^{-1} q
    /// ```
    /// It is worth mentioning that, due to efficiency reasons, the following variables are defined
    /// within the code:
    /// * $`\overline{K}^{\star} = \Delta t \times \overline{C}^{-1}\overline{K}`$
    /// * $`\overline{C}^{\star} = \Delta t \times \overline{C}^{-1}`$
    fn build_static_massive_thermal_network(
        &mut self,
        gases: &[Gas],
        construction: &Rc<Construction>,
        // first_massive: usize,
        // last_massive: usize,
        dt: Float,
        // all_nodes: usize,
        // n_elements: &[usize],
        r_front: Float,
        r_back: Float,
    ) {
        // if this happens, we are trying to build the
        // thermal network for a non-massive wall... Which
        // does not make sense
        // debug_assert!(first_massive != last_massive);
        // debug_assert_eq!(calc_n_total_nodes(n_elements).unwrap(), all_nodes);

        // check coherence in input data
        // if n_elements.len() != construction.materials.len() {
        //     let err = format!("Mismatch between number of layers in construction ({}) and the number of elements in scheme ({})",construction.materials.len(),n_elements.len());
        //     return Err(err);
        // }


    
        // initialize k_prime as K... we will modify it later
        let mut k_prime = self.calc_k_matrix(
            &vec![0.; self.segments.len()], 
            gases,
            // first_massive,
            // last_massive,
            // n_elements,
            // all_nodes,
            true // massive onlye
        );

        // Calc the masses
        let c = self.calc_c_matrix(
            // construction, all_nodes, n_elements
        );
        let all_nodes = self.segments.len();

        // DIVIDE ONE BY THE OTHER to make K_prime
        let mut c_prime: Vec<Float> = Vec::with_capacity(all_nodes);
        for (i, mass) in c.iter().enumerate() {
            if *mass != 0.0 {
                // Multiply the whole column K by this.
                for j in 0..all_nodes {
                    let old_k_value = k_prime.get(i, j).unwrap();
                    k_prime.set(i, j, dt * old_k_value / mass).unwrap();
                }
                c_prime.push(dt / mass);
            } else {
                unreachable!()
            }
        }

        let front_mat = &construction.materials[0];
        let (front_thermal_absorbtance, front_solar_absorbtance) = match &front_mat.substance {
            Substance::Normal(s) => {
                let thermal = match s.thermal_absorbtance() {
                    Ok(v) => *v,
                    Err(_) => {
                        let v = 0.9;
                        eprintln!("Warning: Substance '{}' has no thermal absorbtance... assuming a value of {}", s.name, v);
                        v
                    }
                };
                let solar = match s.solar_absorbtance() {
                    Ok(v) => *v,
                    Err(_) => {
                        let v = 0.7;
                        eprintln!("Warning: Substance '{}' has no Solar absorbtance... assuming a value of {}", s.name, v);
                        v
                    }
                };
                (thermal, solar)
            }
            Substance::Gas(_) => {
                todo!()
            }
        };

        let back_mat = &construction.materials.last().unwrap(); // There should be at least one.
        let (back_thermal_absorbtance, back_solar_absorbtance) = match &back_mat.substance {
            Substance::Normal(s) => {
                let thermal = match s.thermal_absorbtance() {
                    Ok(v) => *v,
                    Err(_) => {
                        let v = 0.9;
                        eprintln!("Warning: Substance '{}' has no thermal absorbtance... assuming a value of {}", s.name, v);
                        v
                    }
                };
                let solar = match s.solar_absorbtance() {
                    Ok(v) => *v,
                    Err(_) => {
                        let v = 0.7;
                        eprintln!("Warning: Substance '{}' has no Solar absorbtance... assuming a value of {}", s.name, v);
                        v
                    }
                };
                (thermal, solar)
            }
            Substance::Gas(_) => {
                todo!()
            }
        };

        let func = move |
                                temperatures: &Matrix,
                                air_temp_front: Float,
                                air_temp_back: Float,
                                rs_front: Float,
                                rs_back: Float,
                                solar_irradiance_front: Float,
                                solar_irradiance_back: Float,
                                ir_irradiance_front: Float,
                                ir_irradiance_back: Float|
            -> Matrix {
            let full_rs_front = r_front + rs_front;
            let full_rs_back = r_back + rs_back;

            // Sol-air temperature
            let t_front = air_temp_front
                + (front_solar_absorbtance * solar_irradiance_front
                    + front_thermal_absorbtance * ir_irradiance_front)
                    * full_rs_front;
            let t_back = air_temp_back
                + (back_solar_absorbtance * solar_irradiance_back
                    + back_thermal_absorbtance * ir_irradiance_back)
                    * full_rs_back;

            let ts_front = temperatures.get(0, 0).unwrap();
            let ts_back = temperatures.get(all_nodes - 1, 0).unwrap();

            // Calculate: k_i = dt*inv(C) * q + h*inv(C)*k*T
            // But, dt*inv(C) = c_prime | h*inv(C)*k = k_prime
            // --> Calculate: k_i = c_prime * q + k_prime * T

            let mut k_i = k_prime.from_prod_n_diag(temperatures, 3).unwrap();

            // if we are generating heat in any layer (e.g., radiant floor) this
            // would need to change...
            let old_value = k_i.get(0, 0).unwrap();
            k_i.set(
                0,
                0,
                old_value
                /* Add RSFront */ - c_prime[0] * ts_front / full_rs_front
                /* And the heat flow*/ + c_prime[0] * t_front / full_rs_front,
            )
            .unwrap();

            // Ki is a column
            let old_value = k_i.get(all_nodes - 1, 0).unwrap();
            k_i.set(
                all_nodes - 1,
                0,
                old_value
                /* Add RSBack */ - c_prime[all_nodes - 1] * ts_back / full_rs_back
                /* And the heat flow*/ + c_prime[all_nodes - 1] * t_back / full_rs_back,
            )
            .unwrap();

            // return
            k_i
        };

        let march_closure = move |
                                temperatures: &mut Matrix,
                                air_temp_front: Float,
                                air_temp_back: Float,
                                rs_front: Float,
                                rs_back: Float,
                                solar_irradiance_front: Float,
                                solar_irradiance_back: Float,
                                ir_irradiance_front: Float,
                                ir_irradiance_back: Float|{
            // First
            let mut k1 = func(
                &temperatures,
                air_temp_front,
                air_temp_back,
                rs_front,
                rs_back,
                solar_irradiance_front,
                solar_irradiance_back,
                ir_irradiance_front,
                ir_irradiance_back,
            );

            // returning "temperatures + k1" is Euler... continuing is
            // Runge–Kutta 4th order

            // Second
            let mut aux = &k1 * 0.5;
            aux += &temperatures;
            let mut k2 = func(
                &aux,
                air_temp_front,
                air_temp_back,
                rs_front,
                rs_back,
                solar_irradiance_front,
                solar_irradiance_back,
                ir_irradiance_front,
                ir_irradiance_back,
            );

            // Third... put the result into `aux`
            k2.scale_into(0.5, &mut aux).unwrap(); //  aux = k2 /2
            aux += &temperatures;
            let mut k3 = func(
                &aux,
                air_temp_front,
                air_temp_back,
                rs_front,
                rs_back,
                solar_irradiance_front,
                solar_irradiance_back,
                ir_irradiance_front,
                ir_irradiance_back,
            );

            // Fourth... put the result into `aux`
            k3.scale_into(1., &mut aux).unwrap(); //  aux = k3
            aux += &temperatures;
            let mut k4 = func(
                &aux,
                air_temp_front,
                air_temp_back,
                rs_front,
                rs_back,
                solar_irradiance_front,
                solar_irradiance_back,
                ir_irradiance_front,
                ir_irradiance_back,
            );

            // Scale them and add them all up
            k1 /= 6.;
            k2 /= 3.;
            k3 /= 3.;
            k4 /= 6.;

            k1 += &k2;
            k1 += &k3;
            k1 += &k4;

            // Let's add it to the temperatures.
            // temperatures.add_to_this(&k1).unwrap();
            *temperatures += &k1;
        };
        self.march = Some(Box::new(march_closure));
        
    }


    /// Calculates the temperatures 
    fn build_static_no_mass_thermal_network(
        &mut self,
        gases: &[Gas],
        construction: &Rc<Construction>,
        // first_massive: usize,
        // last_massive: usize,
        dt: Float,
        // all_nodes: usize,
        // n_elements: &[usize],
        r_front: Float,
        r_back: Float,
    ) {
        let r = self.r_between(0, self.segments.len(), &vec![], gases);
        let u = 1./(r + r_front + r_back);

        let n_nodes = self.segments.len();
        
        let temps = vec![0.0; n_nodes];
        let k = self.calc_k_matrix(&temps, gases, false);// Not only massive ones
        // k *= -1.;
        println!("K={}", k);

        let march_closure = move |
        temperatures: &mut Matrix,
        air_temp_front: Float,
        air_temp_back: Float,
        rs_front: Float,
        rs_back: Float,
        _solar_irradiance_front: Float,
        _solar_irradiance_back: Float,
        _ir_irradiance_front: Float,
        _ir_irradiance_back: Float|{
            
            // Solve (A * temp = -q)
            let mut q = vec![0.0; n_nodes];
            q[0] -= air_temp_front/rs_front;
            q[n_nodes-1] -= air_temp_back/rs_back;
            let q = Matrix::from_data(n_nodes, 1, q);
            dbg!(rs_back);
            dbg!(rs_front);

            // Add the effect of convective heat transfer
            // let ts_front = temperatures.get(0, 0).unwrap();
            // let ts_back = temperatures.get(n_nodes - 1, 0).unwrap();

            // I don't like this... but it seems like I need this.
            let mut kp = k.clone();
            println!("K before = {}",kp);
            

            let old_value = kp.get(0, 0).unwrap();
            kp.set(
                0,
                0,
                old_value
                /* Add RSFront */ -  1. / rs_front
                // /* And the heat flow*/ + c_prime[0] * air_t_front / rs_front,
            )
            .unwrap();

            let old_value = kp.get(n_nodes - 1, n_nodes - 1).unwrap();
            kp.set(
                n_nodes - 1,
                n_nodes - 1,
                old_value
                /* Add RSBack */ -  1. / rs_back
                // /* And the heat flow*/ + c_prime[n_nodes - 1] * air_temp_back / rs_back,
            )
            .unwrap();



            dbg!(air_temp_front);
            dbg!(air_temp_back);
            println!("Q = {}",q);
            println!("K = {}",kp);
            
            // solve K * temp = -q
            let r = kp.gauss_seidel(&q, temperatures, 100, 0.002);
            if r.is_err(){                
                k.gauss_seidel(&q, temperatures, 100, 0.2).unwrap()                
            }
            

        };
        self.march = Some(Box::new(march_closure));
        
        
        
    }

    pub fn build_thermal_network(
        &mut self,
        gases: &[Gas],
        construction: &Rc<Construction>,
        // first_massive: usize,
        // last_massive: usize,
        dt: Float,
        // all_nodes: usize,
        // n_elements: &[usize],
        r_front: Float,
        r_back: Float,
    ) {

        // Is it Linear        
        if self.is_massive && self.is_static {
            // Normal surface.                        
            self.build_static_massive_thermal_network(gases, construction, dt, r_front, r_back);
            
        } else if !self.is_massive && self.is_static {
            self.build_static_no_mass_thermal_network(gases, construction, dt, r_front, r_back);
        };
    }

}



/***********/
/* TESTING */
/***********/

#[cfg(test)]
mod testing {
    use super::*;
    

    fn get_normal(thermal_cond: Float, density: Float, cp: Float, thickness: Float)->Rc<Construction>{
        
        let mut s = simple_model::substance::Normal::new("the substance".into());
        s.set_thermal_conductivity(thermal_cond)
            .set_density(density)
            .set_specific_heat_capacity(cp);
        let s = s.wrap();

        let material = simple_model::Material::new("the mat".into(), s, thickness);
        let material = Rc::new(material);
        let mut construction = simple_model::Construction::new("the construction".into());
        construction.materials.push(material);
        Rc::new(construction)
    }

    #[test]
    fn test_normal_march_closure(){
        let thermal_cond = 1.;
        let density = 2.1;
        let cp = 1.312;
        let thickness = 12.5/1000.;        
        let tstep_sub = 10;

        let gases = Vec::new();
        let construction = get_normal(thermal_cond, density, cp, thickness);
        let n_elements = vec![1]; 
        let mut d = Discretization::build(&construction, &gases, tstep_sub, &n_elements);
        
        let first_massive: usize = 0;
        let last_massive: usize = 1;
        let dt: Float = 0.1;
        let all_nodes: usize = 1;
        let r_front: Float = 0.2;
        let r_back: Float = 0.1;

        let f = d.build_thermal_network(&gases, &construction, dt, r_front, r_back);
    }
    
    #[test]
    fn test_build_k(){
        assert!(false)   
    }

    #[test]
    fn test_build_c(){
        assert!(false)   
        
    }

    #[test]
    fn test_build_normal_mass(){
        let thermal_cond = 1.;
        let density = 2.1;
        let cp = 1.312;
        let thickness = 12.5/1000.;        
        let tstep_sub = 10;

        let gases = Vec::new();
        let construction = get_normal(thermal_cond, density, cp, thickness);
        let d = Discretization::build(&construction, &gases, tstep_sub, &[1]);
        // normal --> linear
        assert!(d.is_static);
        assert_eq!(d.tstep_subdivision, tstep_sub);
        assert_eq!(d.segments.len(), 2);
        // mass of node 0
        let exp_mass = thickness * density * cp/2.;
        let mass = d.segments[0].0;
        assert!( (exp_mass - mass).abs() < 1e-17 , "Expecting mass to be {exp_mass}... found {mass}");
        if let Resistance::Solid(r) = d.segments[0].1{
            assert!( (r - thickness/thermal_cond).abs() < 1e-16 )
        }else{
            panic!("Expecting Solid!")
        }

        let mass = d.segments[1].0;
        assert!( (exp_mass - mass).abs() < 1e-17 , "Expecting mass to be {exp_mass}... found {mass}");
        if let Resistance::Back = d.segments[1].1{
            assert!(true)
        }else{
            panic!("Expecting Back!")
        }
    }

    #[test]
    fn test_build_normal_no_mass(){
        
        let thermal_cond = 1.;
        let density = 2.1;
        let cp = 1.312;
        let thickness = 12.5/1000.;        
        let tstep_sub = 10;
        
        let construction = get_normal(thermal_cond, density, cp, thickness);
        let gases = Vec::new();
        
        let d = Discretization::build(&construction, &gases, tstep_sub, &[0]);

        // normal --> linear
        assert!(d.is_static);
        assert_eq!(d.tstep_subdivision, tstep_sub);
        assert_eq!(d.segments.len(), 2);
        // mass of node 0
        let exp_mass = 0.0;
        let mass = d.segments[0].0;
        assert!( (exp_mass - mass).abs() < 1e-17 , "Expecting mass to be {exp_mass}... found {mass}");
        if let Resistance::Solid(r) = d.segments[0].1{
            assert!( (r - thickness/thermal_cond).abs() < 1e-16 )
        }else{
            panic!("Expecting Solid!")
        }

        let mass = d.segments[1].0;
        assert!( (exp_mass - mass).abs() < 1e-17 , "Expecting mass to be {exp_mass}... found {mass}");
        if let Resistance::Back = d.segments[1].1{
            assert!(true)
        }else{
            panic!("Expecting Back!")
        }
    }

    #[test]
    fn test_build_normal_gas_normal_mass(){
        let thermal_cond = 1.;
        let density = 2.1;
        let cp = 1.312;
        let thickness = 12.5/1000.;
        let mut gases = Vec::new();

        let tstep_sub = 10;
        let mut construction = simple_model::Construction::new("the construction".into());
        // add normal
        let mut normal = simple_model::substance::Normal::new("the substance".into());
        normal.set_thermal_conductivity(thermal_cond)
            .set_density(density)
            .set_specific_heat_capacity(cp);
        let normal = normal.wrap();
        let normal = simple_model::Material::new("the mat".into(), normal, thickness);
        let normal = Rc::new(normal);
        construction.materials.push(normal.clone());

        // add gas
        let mut gas = simple_model::substance::Gas::new("the gas".into());        
        gas.set_gas(simple_model::substance::gas::StandardGas::Air);
        gases.push(gas.clone().into());
        let gas = gas.wrap();
        let gas = simple_model::Material::new("the_gas".into(), gas, thickness);
        let gas = Rc::new(gas);
        construction.materials.push(gas);
        
        // add normal        
        construction.materials.push(normal);
        


        let construction = Rc::new(construction);
        let d = Discretization::build(&construction, &gases, tstep_sub, &[1, 1, 1]);

        // has gas --> linear
        assert!(!d.is_static);
        assert_eq!(d.tstep_subdivision, tstep_sub);
        assert_eq!(d.segments.len(), 4); // normal, gas, normal, back
        
        // node 0
        let exp_mass = thickness * density * cp/2.;
        let mass = d.segments[0].0;
        assert!( (exp_mass - mass).abs() < 1e-17 , "Expecting mass to be {exp_mass}... found {mass}");                
        if let Resistance::Solid(r) = d.segments[0].1 {
            assert!( (r - thickness/thermal_cond).abs() < 1e-16 )
        }else{
            panic!("Expecting Solid!")
        }

        // node 1        
        let mass = d.segments[1].0;
        assert!( (exp_mass - mass).abs() < 1e-17 , "Expecting mass to be {exp_mass}... found {mass}");
        if let Resistance::Cavity(i) = d.segments[1].1{
            // assert!( (r - thickness/thermal_cond).abs() < 1e-16 )
            assert_eq!(i, 0, "Expecting 0... found {i}");
        }else{
            panic!("Expecting Solid!")
        }

        // node 2
        let mass = d.segments[2].0;
        assert!( (exp_mass - mass).abs() < 1e-17 , "Expecting mass to be {exp_mass}... found {mass}");
        if let Resistance::Solid(r) = d.segments[2].1{
            assert!( (r - thickness/thermal_cond).abs() < 1e-16 )
        }else{
            panic!("Expecting Solid!")
        }

        // node 3
        let mass = d.segments[3].0;
        assert!( (exp_mass - mass).abs() < 1e-17 , "Expecting mass to be {exp_mass}... found {mass}");
        if let Resistance::Back = d.segments[3].1{
            assert!(true)
        }else{
            panic!("Expecting Solid!")
        }
    }

    #[test]
    fn test_build_normal_gas_normal_no_mass(){
        let thermal_cond = 1.;
        let density = 2.1;
        let cp = 1.312;
        let thickness = 12.5/1000.;
        let mut gases = Vec::new();

        let tstep_sub = 10;
        let mut construction = simple_model::Construction::new("the construction".into());
        // add normal
        let mut normal = simple_model::substance::Normal::new("the substance".into());
        normal.set_thermal_conductivity(thermal_cond)
            .set_density(density)
            .set_specific_heat_capacity(cp);
        let normal = normal.wrap();
        let normal = simple_model::Material::new("the mat".into(), normal, thickness);
        let normal = Rc::new(normal);
        construction.materials.push(normal.clone());

        // add gas
        let mut gas = simple_model::substance::Gas::new("the gas".into());        
        gas.set_gas(simple_model::substance::gas::StandardGas::Air);
        gases.push(gas.clone().into());
        let gas = gas.wrap();
        let gas = simple_model::Material::new("the_gas".into(), gas, thickness);
        let gas = Rc::new(gas);
        construction.materials.push(gas);
        
        // add normal        
        construction.materials.push(normal);
        


        let construction = Rc::new(construction);
        let d = Discretization::build(&construction, &gases,  tstep_sub, &[0,0,0]);

        // has gas --> linear
        assert!(!d.is_static);
        assert_eq!(d.tstep_subdivision, tstep_sub);
        assert_eq!(d.segments.len(), 4); // normal, gas, normal, back
        
        // node 0
        let exp_mass = 0.0;
        let mass = d.segments[0].0;
        assert!( (exp_mass - mass).abs() < 1e-17 , "Expecting mass to be {exp_mass}... found {mass}");                
        if let Resistance::Solid(r) = d.segments[0].1 {
            assert!( (r - thickness/thermal_cond).abs() < 1e-16 )
        }else{
            panic!("Expecting Solid!")
        }

        // node 1        
        let mass = d.segments[1].0;
        assert!( (exp_mass - mass).abs() < 1e-17 , "Expecting mass to be {exp_mass}... found {mass}");
        if let Resistance::Cavity(i) = d.segments[1].1{
            // assert!( (r - thickness/thermal_cond).abs() < 1e-16 )
            assert_eq!(i, 0, "Expecting 0... found {i}");
        }else{
            panic!("Expecting Solid!")
        }

        // node 2
        let mass = d.segments[2].0;
        assert!( (exp_mass - mass).abs() < 1e-17 , "Expecting mass to be {exp_mass}... found {mass}");
        if let Resistance::Solid(r) = d.segments[2].1{
            assert!( (r - thickness/thermal_cond).abs() < 1e-16 )
        }else{
            panic!("Expecting Solid!")
        }

        // node 3
        let mass = d.segments[3].0;
        assert!( (exp_mass - mass).abs() < 1e-17 , "Expecting mass to be {exp_mass}... found {mass}");
        if let Resistance::Back = d.segments[3].1{
            assert!(true)
        }else{
            panic!("Expecting Solid!")
        }
    }



    
    
    
}
