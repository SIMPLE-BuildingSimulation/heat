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


use crate::cavity::Cavity;
use crate::environment::Environment;
use crate::{ Float, SIGMA};
use matrix::Matrix;
use simple_model::{Construction, Substance};
use std::rc::Rc;


/// Represents a thermal connection in the thermal network.
/// It can be a Cavity, a Solid, or other.
#[derive(Debug, Clone)]
pub enum UValue {
    /// A normal (i.e., $`\lambda/\Delta x`$) U-value
    Solid(Float),

    /// A cavity, comprised of a gas
    Cavity(Box<Cavity>),

    /// The resistance is a surface coefficient.
    Back,

    /// Undefined yet
    None,
}

impl UValue {

    /// Gets the U-value of a `UValue` object
    pub fn u_value(&self, t_before: Float, t_after: Float) -> Float {
        match self {
            Self::Solid(u) => *u,
            Self::Cavity(c) => c.u_value(t_before, t_after),
            Self::Back => 0., // This should be calculated appart
            // Self::Front => {dbg!("Calc front Rs"); 0.1},
            Self::None => panic!("Attempting to get the u-value of None"),
        }
    }
}

impl std::default::Default for UValue {
    fn default() -> Self {
        UValue::None
    }
}




/// Represents the discretization of a [`Construction`] for heat transfer
/// calculation purposes. 
/// 
/// # Note
/// 
/// This object contains all the [`Cavity`] objects in it, which
/// contain information not only about their thickness but also their orientation.
/// This means that one `Discretization` should exist per `Surface`, not just by `Construction`.
pub struct Discretization {
    /// Contains the node mass and the `Resistance` of each segment
    pub segments: Vec<(Float, UValue)>,

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
    
}

impl Discretization {
    /// Creates a new `Discretization`.
    ///
    /// It first calculates the `tstep_subdivision` and the number of elements
    /// on each layer of Construction by calling `discretize_construction()`; and then builds the
    /// `Discretization` by calling `build()`.
    pub fn new(
        construction: &Rc<Construction>,
        model_dt: Float,
        max_dx: Float,
        min_dt: Float,
        height: Float,
        angle: Float,
    ) -> Result<Self, String> {
        let (tstep_subdivision, n_elements) =
            Self::discretize_construction(construction, model_dt, max_dx, min_dt);
        Self::build(construction, tstep_subdivision, &n_elements, height, angle)
    }

    

    /// Creates the `segments` of the `Discretization`.
    fn build(
        construction: &Rc<Construction>,
        tstep_subdivision: usize,
        n_elements: &[usize],
        height: Float,
        angle: Float,
    ) -> Result<Self, String> {
        debug_assert_eq!(n_elements.len(), construction.materials.len());

        // Let's start with an empty set of segments
        let mut n_nodes: usize = n_elements.iter().sum();
        n_nodes = n_nodes.max(construction.materials.len()); // At least one per layer... but Zero means  "no_mass wall"

        let mut segments: Vec<(Float, UValue)> = vec![(0.0, UValue::default()); n_nodes + 1];

        let mut n_segment = 0;
        for (n_layer, n) in n_elements.iter().enumerate() {
            let mut n = *n;
            let material = &construction.materials[n_layer];

            // get the mass of each segment.
            let mass = if n == 0 {
                0.0
            } else {
                match &material.substance {
                    Substance::Normal(s) => {
                        let dx = material.thickness / n as Float;
                        let rho = s.density().expect(
                            "Trying to calculate C_Matrix with a substance without 'density'",
                        );
                        let cp = s.specific_heat_capacity().expect("Trying to calculate C_Matrix with a substance without 'specific heat capacity'");
                        rho * cp * dx
                        
                    }
                    Substance::Gas(_s) => 0.,
                }
            };

            if n == 0 {
                n = 1;
            }
            // Now iterate all segments...
            for _ in 0..n {
                match &material.substance {
                    Substance::Normal(s) => {
                        // Add mass to this and next nodes (if it is NoMass, it is Zero)
                        segments[n_segment].0 += mass / 2.;
                        segments[n_segment + 1].0 += mass / 2.;

                        // Add resistance                        
                        let dx = material.thickness / n as Float;
                        let k = s.thermal_conductivity().unwrap_or_else(|_| panic!("Substance '{}' in material '{}' in Construction '{}' has no thermal conductivity, but we need it", s.name(), material.name(), construction.name()));
                        // Push U-value
                        segments[n_segment].1 = UValue::Solid(k / dx);
                    }
                    Substance::Gas(s) => {
                        dbg!("todo: Fill geometrical elements of Cavity properly");
                        // Search by name
                        let gas = match s.gas() {
                            Ok(simple_model::substance::gas::StandardGas::Air) => {
                                crate::gas::Gas::air()
                            }
                            Ok(simple_model::substance::gas::StandardGas::Argon) => {
                                crate::gas::Gas::argon()
                            }
                            Ok(simple_model::substance::gas::StandardGas::Xenon) => {
                                crate::gas::Gas::xenon()
                            }
                            Ok(simple_model::substance::gas::StandardGas::Krypton) => {
                                crate::gas::Gas::krypton()
                            }
                            _ => {
                                return Err(format!(
                                    "Substance '{}' does not have a standard gas.",
                                    &material.substance.name()
                                ))
                            }
                        };
                        if n_layer == 0 {
                            dbg!("This should be checked earlier.");
                            return Err(format!(
                                "Construction '{}' has a Gas as its first layer",
                                construction.name
                            ));
                        }
                        let prev_mat = construction.materials.get(n_layer - 1).unwrap(); // we already checked this
                        let next_mat = match construction.materials.get(n_layer + 1) {
                            Some(v) => v,
                            None => {
                                return Err(format!(
                                    "Construction '{}' has a Gas as its last layer",
                                    construction.name
                                ))
                            }
                        };

                        const DEFAULT_EM : Float = 0.84;
                        let ein = match &next_mat.substance{
                            Substance::Normal(s)=>match s.thermal_absorbtance(){
                                Ok(v)=>*v,
                                Err(_)=>{
                                    eprintln!("Substance '{}' has no thermal absorbtance... assuming {}", &construction.materials[0].substance.name(), DEFAULT_EM);
                                    DEFAULT_EM
                                }
                            },
                            Substance::Gas(_)=>return Err(format!("Construction '{}' has two gases without a solid layer between them", construction.name))
                        };

                        let eout = match &prev_mat.substance{
                            Substance::Normal(s)=>match s.thermal_absorbtance(){
                                Ok(v)=>*v,
                                Err(_)=>{
                                    eprintln!("Substance '{}' has no thermal absorbtance... assuming {}", &construction.materials[0].substance.name(), DEFAULT_EM);
                                    DEFAULT_EM
                                }
                            },
                            Substance::Gas(_)=>return Err(format!("Construction '{}' has two gases without a solid layer between them", construction.name))
                        };

                        let c = Cavity {
                            gas,
                            thickness: material.thickness,
                            height,
                            angle,
                            eout,
                            ein,
                        };
                        segments[n_segment].1 = UValue::Cavity(Box::new(c));
                    }
                }
                n_segment += 1;
            }
            // Process last one
            segments[n_nodes].1 = UValue::Back;
        }

        let is_static = segments
            .iter()
            .enumerate()
            .all(|(index, (_mass, resistance))| {
                // Last one does not count
                if index == n_nodes {
                    true
                } else if let UValue::Solid(_) = resistance {
                    true
                } else {
                    false
                }
            });
        let is_massive = segments.iter().any(|(mass, _resistance)| *mass > 0.0);

        Ok(Self {
            is_static,
            is_massive,
            segments,
            tstep_subdivision,            
        })
    }

    

    /// Calculates the R value of the whole system
    /// 
    /// # Panics
    /// Panics if the calculated R value is Zero (i.e., if there are no 
    /// layers or something like that)
    pub fn r_value(&self)->Float{
        let mut r = 0.0;

        for (_, u_value) in &self.segments{
            r += match u_value {
                UValue::Cavity(_c)=>todo!(), //c.u_value(t_front, t_back),
                UValue::Solid(v)=>1./v,
                UValue::Back => 0.0,
                UValue::None => unreachable!()
            }
        }

        assert!(r > 0.0, "Found Zero r-value");
        r
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
    /// We are solving the following equation:
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
            const MAX_RS: Float = 0.05;

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
                let b_coef = -dt / (rho * cp * MAX_RS);
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
                    let lambda1 = -dt / (MAX_RS * rho * cp * dx);
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

   
   
    /// Produces $`\overline{K}`$ and $`\vec{q}`$ (as in the equation $`\overline{C} \dot{\vec{T}} =  \overline{K} \vec{T} + \vec{q}`$),
    /// allowing to update the node temperatures in a surface. This method returns the matrix $`\overline{k}`$
    /// representing the thermal network, and the vector $`\vec{q}`$), accounting for the heat flows induced by the
    /// environment.
    ///
    /// # How it works
    ///
    /// This matrix is constructed based ona discretization of the layers of the
    /// construction; that is, each layer is subdivided into `n` elements all of equal
    /// thickness. One node is placed at the beginning and end of each elements. Each element can
    /// be represented by the 2x2 matrix
    ///
    /// ```math
    /// \overline{K}=\begin{bmatrix} -U & U \\
    /// U & -U
    /// \end{bmatrix}   
    ///```
    ///
    /// Then—ignoring all external inputs (i.e., a system disconnected from the environment—the K matrix
    /// for a construction subdivided into 3 elements (i.e., 4 nodes) can be written as follows:
    ///
    /// ```math
    /// \overline{K}=\begin{bmatrix}
    /// -U_{1\rightarrow2} & U_{1\rightarrow2} & 0 & 0 \\
    /// U_{1\rightarrow2} & -U_{1\rightarrow2} - U_{2\rightarrow3} & U_{2\rightarrow3} & 0 \\
    /// 0 & U_{2\rightarrow3} & -U_{2\rightarrow3} - U_{3\rightarrow4} & U_{3\rightarrow4} \\
    /// 0 & 0 & U_{3\rightarrow4} & -U_{3\rightarrow4} \\
    /// \end{bmatrix}   
    /// ```
    ///
    /// Now, these nodes are also connected to a front and back border conditions.
    /// This means that the Matrix $`\overline{K}`$ needs to become:
    ///
    /// ```math
    /// \overline{K}=\begin{bmatrix}
    /// -U_{1\rightarrow2} - h_{so} & U_{1\rightarrow2} & 0 & 0 \\
    /// U_{1\rightarrow2} & -U_{1\rightarrow2} - U_{2\rightarrow3} & U_{2\rightarrow3} & 0 \\
    /// 0 & U_{2\rightarrow3} & -U_{2\rightarrow3} - U_{3\rightarrow4} & U_{3\rightarrow4} \\
    /// 0 & 0 & U_{3\rightarrow4} & -U_{3\rightarrow4}- h_{si} \\
    /// \end{bmatrix}   
    /// ```
    ///
    /// Additionally, the vector $`\vec{q}`$ will have to account for the border conditions. How this is done depends
    /// on the border condition. If the border leads to a zone, a value of $` \epsilon_s + T_{env} h_s + E_{ir}\epsilon_s - \sigma {T_s}^4`$
    /// should be added to $`\vec{q}`$.  Note that $`T_{env}`$ is the temparture of the air in the environment,
    /// $`h_s`$ is the convection coefficient, $`E_{ir}`$ is the incident infrared radiation and $`\epsilon_s`$ is the
    /// emmisivity of the surface. On the contrary, if the border condition is a cavity, then a value of $`T_{pane} U_{cavity}`$
    /// should be added. $`T_{pane}`$ is the temperature of the surface before or after.
    ///         
    pub fn get_k_q(
        &self,
        ini: usize,
        fin: usize,
        temperatures: &Matrix,
        front_env: &Environment,
        front_emmisivity: Float,
        front_hs: Float,
        back_env: &Environment,
        back_emmisivity: Float,
        back_hs: Float
    ) -> (Matrix, Matrix) {
        let (nrows, ncols) = temperatures.size();
        assert_eq!(
            ncols, 1,
            "Expecting `temperatures` to be a Column matrix... but found {} columns",
            ncols
        );
        assert_eq!(
            nrows,
            self.segments.len(),
            "Expecting `temperatures` to have {} elements... found {}",
            self.segments.len(),
            nrows
        );
        let nnodes = fin - ini + 1;        

        let mut k = Matrix::new(0.0, nnodes, nnodes);
        let mut q = Matrix::new(0.0, nnodes, 1);


        // this is just quite helpful
        let get_t_after = |i: usize|->Float{
            match temperatures.get(i + 1, 0) {
                Ok(v) => v,
                Err(_) => {
                    let (.., s) = &self.segments[i];
                    assert!(!matches!(s, &UValue::Cavity(_)));
                    // If it is not Cavity, then it we can return whatever we want,
                    // as this u-value does not depend on tempeature
                    123.
                }
            }
        };


        for local_i in 0..nnodes-1 {
            let global_i = ini + local_i;            
            let (.., uvalue) = &self.segments[global_i];
            let t_before = temperatures.get(global_i, 0).unwrap(); // this should NEVER fail
            let t_after = get_t_after(global_i);
            let u = uvalue.u_value(t_before, t_after);

            // Top left... should be there            
            let old_value = k.get(local_i, local_i).unwrap();
            k.set(local_i, local_i, old_value - u).unwrap();            

            // Bottom right, only if within range.            
            if let Ok(old_value) = k.get(local_i + 1, local_i + 1){ 
                k.set(local_i + 1, local_i + 1, old_value - u).unwrap();
            }
            // Top right, only if within range.
            if let Ok(old_value) = k.get(local_i, local_i + 1) {
                k.set(local_i, local_i + 1, old_value + u).unwrap();
            }

            // Bottom left, only if within range.
            if let Ok(old_value) = k.get(local_i + 1, local_i) {
                k.set(local_i + 1, local_i, old_value + u).unwrap();
            }
        }

        // Add front border conditions
        let (hs_front, front_q) = if ini == 0 {            
            let ts = 273.15 + temperatures.get(0, 0).unwrap();
            let front_q = front_env.air_temperature * front_hs  // convection
                + front_env.ir_irrad * front_emmisivity // incident radiation
                - SIGMA * front_emmisivity * ts.powi(4); // outgoing radiation
                
            (front_hs, front_q)
        } else {
            let (.., resistance) = &self.segments[ini-1];
            let t_before = temperatures.get(ini-1, 0).unwrap(); // this should NEVER fail
            let t_after = get_t_after(fin);
            let u = resistance.u_value(t_before, t_after);

            (u, u*t_after)
        };        
        let old_v = q.get(0, 0).unwrap();
        q.set(0, 0, old_v + front_q).unwrap();

        let old_value = k.get(0, 0).unwrap();
        k.set(0, 0, old_value - hs_front).unwrap();

        // Add back border conditions
        let (hs_back, back_q) = if fin == nnodes - 1 {
            let ts = 273.15 + temperatures.get(fin, 0).unwrap();
                                    
            let back_q = back_env.air_temperature * back_hs  // convection
                + back_env.ir_irrad * back_emmisivity // incident radiation
                - SIGMA * back_emmisivity * ts.powi(4); // outgoing radiation

            (back_hs, back_q)
        } else {
            
            let (.., resistance) = &self.segments[fin];
            let t_before = temperatures.get(fin, 0).unwrap(); // this should NEVER fail
            let t_after = get_t_after(fin);
            let u = resistance.u_value(t_before, t_after);   
            dbg!(t_after - t_before, (t_after - t_before)*u, u);         

            (u, u*t_after)
        };
        let old_v = q.get(nnodes - 1, 0).unwrap();
        q.set(nnodes - 1, 0, old_v + back_q).unwrap();
        let old_value = k.get(nnodes - 1, nnodes - 1).unwrap();
        k.set(nnodes - 1, nnodes - 1, old_value - hs_back).unwrap();

        // return
        (k, q)
    }
}

/***********/
/* TESTING */
/***********/

#[cfg(test)]
mod testing {
    use super::*;
    use crate::gas::Gas;



    fn get_normal(
        thermal_cond: Float,
        density: Float,
        cp: Float,
        thickness: Float,
    ) -> Rc<Construction> {
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
    fn test_build_normal_mass() {
        let thermal_cond = 1.;
        let density = 2.1;
        let cp = 1.312;
        let thickness = 12.5 / 1000.;
        let tstep_sub = 10;

        let construction = get_normal(thermal_cond, density, cp, thickness);
        let d = Discretization::build(&construction, tstep_sub, &[1], 1., 0.).unwrap();
        // normal --> linear
        assert!(d.is_static);
        assert_eq!(d.tstep_subdivision, tstep_sub);
        assert_eq!(d.segments.len(), 2);
        // mass of node 0
        let exp_mass = thickness * density * cp / 2.;
        let mass = d.segments[0].0;
        assert!(
            (exp_mass - mass).abs() < 1e-17,
            "Expecting mass to be {exp_mass}... found {mass}"
        );
        if let UValue::Solid(u) = d.segments[0].1 {
            assert!((u -  thermal_cond/thickness).abs() < 1e-16)
        } else {
            panic!("Expecting Solid!")
        }

        let mass = d.segments[1].0;
        assert!(
            (exp_mass - mass).abs() < 1e-17,
            "Expecting mass to be {exp_mass}... found {mass}"
        );
        if let UValue::Back = d.segments[1].1 {
            assert!(true)
        } else {
            panic!("Expecting Back!")
        }
    }

    #[test]
    fn test_build_normal_no_mass() {
        let thermal_cond = 1.;
        let density = 2.1;
        let cp = 1.312;
        let thickness = 12.5 / 1000.;
        let tstep_sub = 10;

        let construction = get_normal(thermal_cond, density, cp, thickness);

        let d = Discretization::build(&construction, tstep_sub, &[0], 1., 0.).unwrap();

        // normal --> linear
        assert!(d.is_static);
        assert_eq!(d.tstep_subdivision, tstep_sub);
        assert_eq!(d.segments.len(), 2);
        // mass of node 0
        let exp_mass = 0.0;
        let mass = d.segments[0].0;
        assert!(
            (exp_mass - mass).abs() < 1e-17,
            "Expecting mass to be {exp_mass}... found {mass}"
        );
        if let UValue::Solid(u) = d.segments[0].1 {
            assert!((u -  thermal_cond/thickness).abs() < 1e-16)
        } else {
            panic!("Expecting Solid!")
        }

        let mass = d.segments[1].0;
        assert!(
            (exp_mass - mass).abs() < 1e-17,
            "Expecting mass to be {exp_mass}... found {mass}"
        );
        if let UValue::Back = d.segments[1].1 {
            assert!(true)
        } else {
            panic!("Expecting Back!")
        }
    }

    #[test]
    fn test_build_normal_gas_normal_mass() {
        let thermal_cond = 1.;
        let density = 2.1;
        let cp = 1.312;
        let thickness = 12.5 / 1000.;

        let tstep_sub = 10;
        let mut construction = simple_model::Construction::new("the construction".into());
        // add normal
        let mut normal = simple_model::substance::Normal::new("the substance".into());
        normal
            .set_thermal_conductivity(thermal_cond)
            .set_density(density)
            .set_thermal_absorbtance(0.9)
            .set_specific_heat_capacity(cp);
        let normal = normal.wrap();
        let normal = simple_model::Material::new("the mat".into(), normal, thickness);
        let normal = Rc::new(normal);
        construction.materials.push(normal.clone());

        // add gas
        let mut gas = simple_model::substance::Gas::new("the gas".into());
        gas.set_gas(simple_model::substance::gas::StandardGas::Air);

        let gas = gas.wrap();
        let gas = simple_model::Material::new("the_gas".into(), gas, thickness);
        let gas = Rc::new(gas);
        construction.materials.push(gas);

        // add normal
        construction.materials.push(normal);

        let construction = Rc::new(construction);
        let d = Discretization::build(&construction, tstep_sub, &[1, 1, 1], 1., 0.).unwrap();

        // has gas --> linear
        assert!(!d.is_static);
        assert_eq!(d.tstep_subdivision, tstep_sub);
        assert_eq!(d.segments.len(), 4); // normal, gas, normal, back

        // node 0
        let exp_mass = thickness * density * cp / 2.;
        let mass = d.segments[0].0;
        assert!(
            (exp_mass - mass).abs() < 1e-17,
            "Expecting mass to be {exp_mass}... found {mass}"
        );
        if let UValue::Solid(u) = d.segments[0].1 {
            assert!((u - thermal_cond/thickness).abs() < 1e-16)
        } else {
            panic!("Expecting Solid!")
        }

        // node 1
        let mass = d.segments[1].0;
        assert!(
            (exp_mass - mass).abs() < 1e-17,
            "Expecting mass to be {exp_mass}... found {mass}"
        );
        if let UValue::Cavity(c) = &d.segments[1].1 {            
            let u = c.u_value(10., 10.);
            dbg!(u);
        } else {
            panic!("Expecting a Cavity!")
        }

        // node 2
        let mass = d.segments[2].0;
        assert!(
            (exp_mass - mass).abs() < 1e-17,
            "Expecting mass to be {exp_mass}... found {mass}"
        );
        if let UValue::Solid(u) = d.segments[2].1 {
            assert!((u - thermal_cond/thickness).abs() < 1e-16)
        } else {
            panic!("Expecting Solid!")
        }

        // node 3
        let mass = d.segments[3].0;
        assert!(
            (exp_mass - mass).abs() < 1e-17,
            "Expecting mass to be {exp_mass}... found {mass}"
        );
        if let UValue::Back = d.segments[3].1 {
            assert!(true)
        } else {
            panic!("Expecting Solid!")
        }
    }

    #[test]
    fn test_build_normal_gas_normal_no_mass() {
        let thermal_cond = 1.;
        let density = 2.1;
        let cp = 1.312;
        let thickness = 12.5 / 1000.;

        let tstep_sub = 10;
        let mut construction = simple_model::Construction::new("the construction".into());
        // add normal
        let mut normal = simple_model::substance::Normal::new("the substance".into());
        normal
            .set_thermal_conductivity(thermal_cond)
            .set_density(density)
            .set_specific_heat_capacity(cp);
        let normal = normal.wrap();
        let normal = simple_model::Material::new("the mat".into(), normal, thickness);
        let normal = Rc::new(normal);
        construction.materials.push(normal.clone());

        // add gas
        let mut gas = simple_model::substance::Gas::new("the gas".into());
        gas.set_gas(simple_model::substance::gas::StandardGas::Air);

        let gas = gas.wrap();
        let gas = simple_model::Material::new("the_gas".into(), gas, thickness);
        let gas = Rc::new(gas);
        construction.materials.push(gas);

        // add normal
        construction.materials.push(normal);

        let construction = Rc::new(construction);
        let d = Discretization::build(&construction, tstep_sub, &[0, 0, 0], 1., 0.).unwrap();

        // has gas --> linear
        assert!(!d.is_static);
        assert_eq!(d.tstep_subdivision, tstep_sub);
        assert_eq!(d.segments.len(), 4); // normal, gas, normal, back

        // node 0
        let exp_mass = 0.0;
        let mass = d.segments[0].0;
        assert!(
            (exp_mass - mass).abs() < 1e-17,
            "Expecting mass to be {exp_mass}... found {mass}"
        );
        if let UValue::Solid(u) = d.segments[0].1 {
            assert!((u - thermal_cond/thickness).abs() < 1e-16)
        } else {
            panic!("Expecting Solid!")
        }

        // node 1
        let mass = d.segments[1].0;
        assert!(
            (exp_mass - mass).abs() < 1e-17,
            "Expecting mass to be {exp_mass}... found {mass}"
        );
        if let UValue::Cavity(c) = &d.segments[1].1 {
            // assert!( (r - thickness/thermal_cond).abs() < 1e-16 )
            let u = c.u_value(10., 20.);
            dbg!(u);
        } else {
            panic!("Expecting Solid!")
        }

        // node 2
        let mass = d.segments[2].0;
        assert!(
            (exp_mass - mass).abs() < 1e-17,
            "Expecting mass to be {exp_mass}... found {mass}"
        );
        if let UValue::Solid(u) = d.segments[2].1 {
            assert!((u -  thermal_cond/thickness).abs() < 1e-16, "u = {} | exp = {}", u, thermal_cond/thickness)
        } else {
            panic!("Expecting Solid!")
        }

        // node 3
        let mass = d.segments[3].0;
        assert!(
            (exp_mass - mass).abs() < 1e-17,
            "Expecting mass to be {exp_mass}... found {mass}"
        );
        if let UValue::Back = d.segments[3].1 {
            assert!(true)
        } else {
            panic!("Expecting Solid!")
        }
    }

    #[test]
    fn test_get_q_k_solid() {
        let thickness = 0.5;
        let n = 5;
        let dx = thickness / n as Float;
        let thermal_cond = 2.12;
        let u = thermal_cond / dx;

        let mut segments = Vec::with_capacity(n + 1);
        for _ in 0..n {
            segments.push((0.0, UValue::Solid(u)));
        }
        segments.push((0.0, UValue::Back));
        assert_eq!(segments.len(), n + 1);

        let d = Discretization {
            segments,
            tstep_subdivision: 1,
            is_static: true,
            is_massive: true,            
        };

        let front_env = Environment {
            solar_radiation: 0.0,
            air_temperature: 1.,
            air_speed: 0.,
            ir_irrad: SIGMA * (273.15 as Float).powi(4),
            .. Environment::default()
        };
        let front_emmisivity = 0.9;

        let back_env = Environment {
            solar_radiation: 0.0,
            air_temperature: 6.,
            air_speed: 0.,
            ir_irrad: SIGMA * (5. + 273.15 as Float).powi(4),
            .. Environment::default()
        };
        let back_emmisivity = 0.9;

        let temperatures = Matrix::from_data(n + 1, 1, vec![0., 1., 2., 3., 4., 5.]);
        let front_hs = front_env.get_hs();
        let back_hs = back_env.get_hs();

        let (k, q) = d.get_k_q(
            0,
            n,
            &temperatures,
            &front_env,
            front_emmisivity,
            front_hs,
            &back_env,
            back_emmisivity,
            back_hs,
        );
        println!("k = {}", k);
        println!("heat_flows = {}", q);

        

        for r in 0..n + 1 {
            // Check q
            let heat_flow = q.get(r, 0).unwrap();
            if r == 0 {
                assert!(heat_flow > 1e-5);
            } else if r == n {
                assert!(heat_flow > 1e-5);
            } else {
                assert!(heat_flow.abs() < 1e-29);
            }

            for c in 0..n + 1 {
                let v = k.get(r, c).unwrap();

                if c == r {
                    // Diagonal
                    if c == 0 {
                        // Exterior node
                        assert!((v - (-front_hs - u)).abs() < 1e-20);
                    } else if c == n {
                        // interior node
                        assert!((v - (-back_hs - u)).abs() < 1e-20);
                    } else {
                        // Any other node.
                        assert!((v - (-2. * u)).abs() < 1e-20);
                    }
                } else {
                    // Not diagonal
                    let r = r as i32;
                    let c = c as i32;
                    if (c - r).abs() <= 1 {
                        // if it is within tri-diagonal
                        assert!((v - u).abs() < 1e-20);
                    } else {
                        // if not
                        assert!(v.abs() < 1e-29);
                    }
                }
            }
        }
    }


    #[test]
    fn test_get_q_k_partial() {
        let thickness = 0.5;
        let n = 5;
        let dx = thickness / n as Float;
        let thermal_cond = 2.12;
        let u = thermal_cond / dx;

        let mut segments = Vec::with_capacity(n + 1);
        for _ in 0..n {
            segments.push((0.0, UValue::Solid(u)));
        }
        segments.push((0.0, UValue::Back));
        assert_eq!(segments.len(), n + 1);

        let d = Discretization {
            segments,
            tstep_subdivision: 1,
            is_static: true,
            is_massive: true,            
        };

        let front_env = Environment {
            air_temperature: 1.,
            air_speed: 0.,
            ir_irrad: SIGMA * (273.15 as Float).powi(4),
            solar_radiation: 0.0,
            .. Environment::default()
        };
        let front_emmisivity = 0.9;

        let back_env = Environment {
            air_temperature: 6.,
            air_speed: 0.,
            ir_irrad: SIGMA * (5. + 273.15 as Float).powi(4),
            solar_radiation: 0.0,
            .. Environment::default()
        };
        let back_emmisivity = 0.9;

        let temperatures = Matrix::from_data(n + 1, 1, vec![0., 1., 2., 3., 4., 5.]);
        let front_hs = 0.1;
        let back_hs = 0.1;
        let (k, q) = d.get_k_q(
            1,
            n-1,
            &temperatures,
            &front_env,
            front_emmisivity,
            front_hs,
            &back_env,
            back_emmisivity,
            back_hs,
        );
        println!("k = {}", k);
        println!("heat_flows = {}", q);
        

        for r in 0..n-1 {
            // Check q
            let heat_flow = q.get(r, 0).unwrap();
            if r == 0 {
                assert!(heat_flow > 1e-5);
            } else if r == n-2 {
                assert!(heat_flow > 1e-5);
            } else {
                assert!(heat_flow.abs() < 1e-29);
            }

            for c in 0..n-1  {
                let v = k.get(r, c).unwrap();

                if c == r {
                    // Diagonal
                    if c == 0 {
                        // Exterior node
                        assert!((v + 2.* u).abs() < 1e-20, "v = {} | found = {}", v, 2.*u);
                    } else if c == n-2 {
                        // interior node
                        assert!((v - (-2. * u)).abs() < 1e-20);
                    } else {
                        // Any other node.
                        assert!((v - (-2. * u)).abs() < 1e-20);
                    }
                } else {
                    // Not diagonal
                    let r = r as i32;
                    let c = c as i32;
                    if (c - r).abs() <= 1 {
                        // if it is within tri-diagonal
                        assert!((v - u).abs() < 1e-20);
                    } else {
                        // if not
                        assert!(v.abs() < 1e-29);
                    }
                }
            }
        }
    }


    #[test]
    fn test_u(){
        // https://github.com/LBNL-ETA/Windows-CalcEngine/blob/main/src/Tarcog/tst/units/DoubleClear_UValueEnvironment.unit.cpp

        
        let gap_thickness = 0.0127;
        let solid_thickness = 0.003048;   // [m]
        let solid_conductance = 1.0; 

        let gap = Cavity{
            thickness: gap_thickness,
            height: 1.,
            gas: Gas::air(),
            eout: 0.84, 
            ein: 0.84,
            angle: 0.*crate::PI/2.,
        };

        let mut segments = Vec::with_capacity(4);
        // layer 1.
        segments.push((0., UValue::Solid(solid_conductance/solid_thickness)));

        // Layer 2: Gap
        segments.push((0.0, UValue::Cavity(Box::new(gap))));

        // layer 3.        
        segments.push((0., UValue::Solid(solid_conductance/solid_thickness)));

        // Last node
        segments.push((0.0, UValue::Back));
        let d = Discretization {
            segments,
            tstep_subdivision: 1,
            is_static: true,
            is_massive: true,            
        };

        // Borders
        let air_temperature = 255.15 - 273.15;   // Kelvins into C
        let air_speed = 5.5;            // meters per second
        let t_sky: Float = 255.15;             // Kelvins into C
        let solar_radiation = 789.0;
        let front_env = Environment{
            air_speed,
            air_temperature,
            solar_radiation,
            ir_irrad: SIGMA * (t_sky.powi(4)),
            .. Environment::default()
        };

        let back_env = Environment{
            air_speed: 0.0,
            air_temperature: 294.15 - 273.15, // K into C
            solar_radiation : 0.,
            ir_irrad: SIGMA * (294.15 as Float).powi(4),
            .. Environment::default()
        };

        let kelvin = Matrix::new(273.14, 4, 1);


        // First, u-value        
        let exp_temps = vec![258.791640 , 259.116115, 279.323983, 279.648458];
        let mut temperatures = Matrix::from_data(4, 1, exp_temps.clone());
        temperatures -= &kelvin;
        
        let front_hs = 36.34359273; // these were found by analyzing WINDOW's response
        let back_hs = 5.;

        let (k, mut q) = d.get_k_q(0, 3, &temperatures, &front_env, 0.9, front_hs, &back_env, 0.9, back_hs);
        q *= -1.;
        // let mut temps = Matrix::new(0.0, 4, 1);
        
        println!("K = {}", k);
        let keep_k = k.clone();
        println!("q = {}", q); 
        let keep_q = q.clone();       
        let t = k.mut_n_diag_gaussian(q, 3).unwrap();        
        
        println!("T = {}", &t + &kelvin);
        
        for (i,exp) in exp_temps.iter().enumerate() {
            let found = t.get(i, 0).unwrap() + 273.15;
            assert!( (exp - found).abs() < 0.17, "Expecting {}, found {}... delta is {}", exp, found, (exp - found).abs() );
        }        
        let mut check =  &keep_k*&t; // this should be equals to keep_q
        check -= &keep_q;
        println!("check = {}", check);


        // Then, with sun
        let exp_temps = vec![261.920088, 262.408524, 284.752662, 285.038190];
        let mut temperatures = Matrix::from_data(4, 1, exp_temps.clone());
        temperatures -= &kelvin; // into C

        let front_hs = 32.6; // these were found by analyzing WINDOW's response
        let back_hs = 6.8;
        let (k, mut q) = d.get_k_q(0, 3, &temperatures, &front_env, 0.9, front_hs, &back_env, 0.9, back_hs);
        
        // add sun?
        let old = q.get(0, 0).unwrap();
        q.set(0, 0, old  + 0.096489921212/2. * solar_radiation).unwrap();
        let old = q.get(1, 0).unwrap();
        q.set(1, 0, old  + 0.096489921212/2. * solar_radiation).unwrap();

        let old = q.get(2, 0).unwrap();
        q.set(2, 0, old  + 0.072256758809/2. * solar_radiation).unwrap();
        let old = q.get(3, 0).unwrap();
        q.set(3, 0, old  + 0.072256758809/2. * solar_radiation).unwrap();


        q *= -1.;

        let keep_k = k.clone();
        // let keep_q = q.clone();
        
        
        println!("K = {}", k);
        println!("q = {}", q);
        
        let t = k.mut_n_diag_gaussian(q, 3).unwrap();
                
        println!("T = {}", &t - &kelvin);
        
        for (i,exp) in exp_temps.iter().enumerate() {
            let found = t.get(i, 0).unwrap() + 273.15;
            assert!( (exp - found).abs() < 0.15, "Expecting {}, found {}... delta is {}", exp, found, (exp - found).abs() );
        }        
        let mut check =  &keep_k*&t; // this should be equals to keep_q
        check -= &t;
        println!("check = {}", check);
    }
}
