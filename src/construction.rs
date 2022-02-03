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

use crate::Float;

use std::rc::Rc;

use matrix::Matrix;

use simple_model::{Construction, Substance};

/// Given a Maximum element thickness ($`\Delta x_{max}`$) and a minimum timestep ($`\Delta t_{min}`$), this function
/// will find an arguibly good (i.e., stable and accurate) combination of $\\Delta t\$ and number of elements in each
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
pub fn discretize_construction(
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
            let (k, rho, cp) = match substance{
                Substance::Normal(s)=>{                    
                    let k = s.thermal_conductivity().expect("Trying to discretize a construction that contains a Normal Substance without a 'thermal conductivity'");            
                    let rho = s.density().expect("Trying to discretize a construction that contains a Normal Substance without a 'density'");            
                    let cp = s.specific_heat_capacity().expect("Trying to discretize a construction that contains a Normal Substance without a 'specific heat capacity'");            
                    (*k, *rho, *cp)
                },
                // Potentially... Substance::Gas(_)=>{n_elements.push(0); continue}
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
                let (k, rho, cp) = match substance{
                    Substance::Normal(s)=>{                    
                        let k = s.thermal_conductivity().unwrap();
                        let rho = s.density().unwrap();
                        let cp = s.specific_heat_capacity().unwrap();
                        (*k, *rho, *cp)
                    },
                    // Potentially... Substance::Gas(_)=>{panic! ?}
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

/// In a discretization scheme, this function finds the
/// first and the last massive layers... the no-mass layers
/// before and after the first and last massive ones, respectively,
/// are just bulked with the values of RSi and RSo.
pub fn get_first_and_last_massive_elements(n_elements: &[usize]) -> Result<(usize, usize), String> {
    // We need somethig to process!
    if n_elements.is_empty() {
        return Err(
            "Impossible to check first and last massive layers in empty discretization scheme"
                .to_string(),
        );
    }
    let n_layers = n_elements.len();

    // Border cases: one layer
    if n_layers == 1 {
        if n_elements[0] == 0 {
            // if the layer is no-mass, then return 0,0
            return Ok((0, 0));
        } else {
            // if the layer has mas,
            return Ok((0, 1));
        }
    }

    // find the first and last massive layers
    let mut first_massive = n_layers;
    let mut last_massive = 0;
    let mut any_massive = false;
    for (i, nodes) in n_elements.iter().enumerate() {
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
    if any_massive {
        Ok((first_massive, last_massive + 1))
    } else {
        Ok((0, 0))
    }
}

/// This function calculates the number of nodes that result
/// from the Construction's discretization.
///
/// For instance, if a construction has only one layer with one
/// element, then that construction has two nodes (one interior and one
/// exterior one). Similarly, if there are two layers but one has no mass,
/// then the first one will add only one node because it will be merged
/// with the exterior film heat transfer coefficient.
///
/// Nodes are placed in the joints between elements, and
/// the nodes for the Rsi and Rso layers are considered.
pub fn calc_n_total_nodes(n_elements: &[usize]) -> Result<usize, String> {
    // We need somethig to process!
    if n_elements.is_empty() {
        return Err("Wrong discretization scheme... it contains no elements".to_string());
    }

    // Border case: Only one element.
    if n_elements.len() == 1 {
        if n_elements[0] == 0 {
            return Ok(2);
        } else {
            return Ok(n_elements[0] + 1);
        }
    }

    // Otherwise, let's do this.
    let mut n: usize = 1; // the end node.
    let (first_massive, last_massive) = get_first_and_last_massive_elements(n_elements)?;

    // No massive layer at all in the construction...
    if first_massive == 0 && last_massive == 0 {
        return Ok(2);
    }

    // Now, process from the first massive to the last.
    // interior and exterior no-mass layers do not contribute
    // nodes.
    let mut in_nomass_layer: bool = false;
    for nodes in &n_elements[first_massive..last_massive] {
        if *nodes > 0 {
            // Material with mass.
            in_nomass_layer = false;
            n += nodes;
        } else {
            // material with no-mass
            if !in_nomass_layer {
                n += 1;
                in_nomass_layer = true;
            }
        }
    }

    // return
    Ok(n)
}

pub fn calc_r_front(c: &Rc<Construction>, first_massive: usize) -> Float {
    let mut r_front = 0.;
    for i in 0..first_massive {        
        let material = &c.materials[i];                                         
        match &material.substance{
            Substance::Normal(s)=>{
                r_front += material.thickness / s.thermal_conductivity().expect("calc_r_front() requires that Normal substances have Thermal Conductivity");

            }
        }
    }
    r_front
}

pub fn calc_r_back(c: &Rc<Construction>, last_massive: usize) -> Float {
    let mut r_back = 0.;
    for i in last_massive..c.materials.len() {
        let material = &c.materials[i];
        match &material.substance{
            Substance::Normal(s)=>{r_back += material.thickness / s.thermal_conductivity().expect("calc_r_back requires that Normal substances have Thermal Conductivity");}
        }
        
    }
    r_back
}

fn calc_c_matrix(
    construction: &Rc<Construction>,
    all_nodes: usize,
    n_elements: &[usize],
) -> Vec<Float> {
    let mut c_matrix: Vec<Float> = vec![0.0; all_nodes];
    let mut node = 0;

    for n_layer in 0..construction.materials.len() {
        // let layer_index = c.get_layer_index(n_layer).unwrap();
        let material = &construction.materials[n_layer]; //model.get_material(layer_index).unwrap();

        let m = n_elements[n_layer];

        if m != 0 {
            // if has mass
            for _ in 0..m {
                // Calc mass
                // let substance_index = material.get_substance_index().unwrap();
                let substance = &material.substance; //model.get_substance(substance_index).unwrap();

                // let rho = substance.density().unwrap();
                // let cp = substance.specific_heat_capacity().unwrap();
                let (rho, cp) = match substance{
                    Substance::Normal(s)=>{                                            
                        let rho = s.density().expect("Trying to calculate C_Matrix with a substance without 'density'");            
                        let cp = s.specific_heat_capacity().expect("Trying to calculate C_Matrix with a substance without 'specific heat capacity'");            
                        (*rho, *cp)
                    },
                    // Potentially... Substance::Gas(_)=>{n_elements.push(0); continue}
                };
                let dx = material.thickness / (m as Float);
                let m = rho * cp * dx; // dt;

                // Nodes are in the joints between
                // finite elements... add half mass from each
                c_matrix[node] += m / 2.0;
                c_matrix[node + 1] += m / 2.0;

                // advance node.
                node += 1;
            }
        } else {
            // else, they are zero already...
            // still, we need to move forward one
            // node.
            if n_layer > 0 && n_elements[n_layer - 1] > 0 {
                // Only do it if the previous
                // layer was massive
                node += 1;
            }
        }
    }
    c_matrix
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
/// This method returns such a matrix, without the $`R_{si, full}`$ and $`R_{so, full}`$. They need to
/// be added when marching (because these values change over time).
fn calc_k_matrix(
    c: &Rc<Construction>,
    first_massive: usize,
    last_massive: usize,
    n_elements: &[usize],
    all_nodes: usize,
) -> Matrix {
    // initialize k_prime
    let mut k_matrix = Matrix::new(0.0, all_nodes, all_nodes);

    // Calculate what is in between the massive layers
    let mut node: usize = 0;
    let mut n_layer: usize = first_massive;

    while n_layer < last_massive {
        // let material_index = c.get_layer_index(n_layer).unwrap();
        let material = &c.materials[n_layer]; //model.get_material(material_index).unwrap();

        let m = n_elements[n_layer];
        if m == 0 {
            // no-mass material
            // add up all the R of the no-mass layers that
            // are together
            let mut thermal_resistance = 0.0; // if the material is no mass, then the first value
            while n_layer < last_massive && n_elements[n_layer] == 0 {
                // let material_index = c.get_layer_index(n_layer).unwrap();
                let material = &c.materials[n_layer]; //model.get_material(material_index).unwrap();
                let dx = material.thickness; //().unwrap();

                // let substance_index = material.get_substance_index().unwrap();                
                let k = match &material.substance {
                    Substance::Normal(s)=>{
                        let k = s.thermal_conductivity().expect("Trying to calc K-matrix of a construciton containing a Normal Substance without 'Thermal conductivity'");
                        k
                    }
                };

                thermal_resistance += dx / k;
                n_layer += 1;
            }

            // update values
            let u_value = 1. / thermal_resistance;
            // top left
            let old_value = k_matrix.get(node, node).unwrap();
            k_matrix.set(node, node, old_value - u_value).unwrap();
            // top right
            let old_value = k_matrix.get(node, node + 1).unwrap();
            k_matrix.set(node, node + 1, old_value + u_value).unwrap();
            // bottom left
            let old_value = k_matrix.get(node + 1, node).unwrap();
            k_matrix.set(node + 1, node, old_value + u_value).unwrap();
            // bottom right
            let old_value = k_matrix.get(node + 1, node + 1).unwrap();
            k_matrix
                .set(node + 1, node + 1, old_value - u_value)
                .unwrap();

            // Move one node ahead
            node += 1;
        } else {
            // calc U value
            
            let k = match &material.substance {
                Substance::Normal(s)=>{
                    let k = s.thermal_conductivity().expect("Trying to calc K-matrix of a construciton containing a Normal Substance without 'Thermal conductivity'");
                    k
                }
            };
            let dx = material.thickness / (m as Float);
            let u = k / dx;

            for _ in 0..m {
                // top left
                let old_value = k_matrix.get(node, node).unwrap();
                k_matrix.set(node, node, old_value - u).unwrap();
                // top right
                let old_value = k_matrix.get(node, node + 1).unwrap();
                k_matrix.set(node, node + 1, old_value + u).unwrap();
                // bottom left
                let old_value = k_matrix.get(node + 1, node).unwrap();
                k_matrix.set(node + 1, node, old_value + u).unwrap();
                // bottom right
                let old_value = k_matrix.get(node + 1, node + 1).unwrap();
                k_matrix.set(node + 1, node + 1, old_value - u).unwrap();

                // advance node.
                node += 1;
            }
            n_layer += 1;
        }
    }

    // // ADD RSI AND RSO
    // // add r_si to first node
    // if r_front > 0.0 {
    //     let old_value = k.get(0, 0).unwrap();
    //     k.set(0, 0, old_value - 1.0 / r_front).unwrap();
    // }

    // // rs_o to the last node
    // if r_back > 0.0 {
    //     let old_value = k.get(all_nodes - 1, all_nodes - 1).unwrap();
    //     k.set(all_nodes - 1, all_nodes - 1, old_value - 1.0 / r_back)
    //         .unwrap();
    // }
    // return
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
pub fn build_thermal_network(
    construction: &Rc<Construction>,
    first_massive: usize,
    last_massive: usize,
    dt: Float,
    all_nodes: usize,
    n_elements: &[usize],
    r_front: Float,
    r_back: Float,
) -> Result<impl Fn(&Matrix, Float, Float, Float, Float) -> Matrix, String> {
    // if this happens, we are trying to build the
    // thermal network for a non-massive wall... Which
    // does not make sense
    debug_assert!(first_massive != last_massive);
    debug_assert_eq!(calc_n_total_nodes(n_elements).unwrap(), all_nodes);

    // check coherence in input data
    if n_elements.len() != construction.materials.len() {
        let err = format!("Mismatch between number of layers in construction ({}) and the number of elements in scheme ({})",construction.materials.len(),n_elements.len());
        return Err(err);
    }

    // initialize k_prime as K... we will modify it later    
    let mut k_prime = calc_k_matrix(
        construction,
        first_massive,
        last_massive,
        n_elements,
        all_nodes,
    );

    // Calc the masses
    let c = calc_c_matrix(construction, all_nodes, n_elements);

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

    // SET C_I AND C_O
    // *c_i = c[0]/dt;
    // *c_o = c[all_nodes - 1]/dt;

    // *given_k_prime = k_prime;

    // let r_front = calc_r_front(construction, first_massive);
    // let r_back = calc_r_back(construction, last_massive);

    let march_closure = move |
            nodes_temps: &Matrix,
            t_front: Float,
            t_back: Float,
            rs_front: Float,
            rs_back: Float|
          -> Matrix {
        let full_rs_front = r_front + rs_front;
        let full_rs_back = r_back + rs_back;
        let ts_front = nodes_temps.get(0, 0).unwrap();
        let ts_back = nodes_temps.get(all_nodes - 1, 0).unwrap();

        // Calculate: k_i = dt*inv(C) * q + h*inv(C)*k*T
        // But, dt*inv(C) = c_prime | h*inv(C)*k = k_prime
        // --> Calculate: k_i = c_prime * q + k_prime * T

        let mut k_i = k_prime.from_prod_n_diag(nodes_temps, 3).unwrap();

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

    Ok(march_closure)
    // Ok(())
}
