#![allow(clippy::too_many_arguments)]

use matrix::Matrix;

//use building_model::construction::Construction;
use building_model::building::Building;
use building_model::construction::Construction;

/// Calculates the R-value of the Construction.
/// # Parameters:
/// * materials: A reference to the vector containing the materials from which this object was built from
pub fn r_value(building: &Building, construction: &Construction) -> Result<f64, String> {
    //let construction = building.get_construction(construction_index).unwrap();

    //let materials = building.get_materials();
    //let substances = building.get_substances();

    let mut r = 0.0;

    for material_index in construction.layers() {
        let material = building.get_material(*material_index).unwrap();
        let substance_index = material.get_substance_index().unwrap();
        let substance = building.get_substance(substance_index).unwrap();
        let lambda = substance.thermal_conductivity().unwrap();

        r += material.thickness().unwrap() / lambda;
    }

    Ok(r)
}

/// Given a Maximum thickness (max_dx) and a minimum timestep (max_dt), this function
/// will find a good combination of dt and number of elements in each
/// layer of the construction.
///
/// This function recursively increases N in order to reduce dt to numbers
/// that divide main_dt as a while (e.g. main_dt/1, main_dt/2, main_dt/3...)
pub fn discretize_construction(
    building: &Building,
    c: &Construction,
    main_dt: f64,
    max_dx: f64,
    min_dt: f64,
) -> (usize, Vec<usize>) {
    // I could only think of how to make this recursively... so I did this.
    fn aux(
        building: &Building,
        c: &Construction,
        main_dt: f64,
        n: usize,
        max_dx: f64,
        min_dt: f64,
    ) -> (usize, Vec<usize>) {
        let dt = main_dt / (n as f64);
        let safety = 1.5;

        // Choose a dx so that dt allows convergence.
        // stability is assured by (alpha * dt / dx^2 <= 1/2 )
        // meaning, we need to satisfy dx >= sqrt(2 * alpha * dt)

        // So, for each layer
        let mut n_elements: Vec<usize> = Vec::with_capacity(c.n_layers());

        for n_layer in 0..c.n_layers() {
            let material_index = c.get_layer_index(n_layer).unwrap();
            let material = building.get_material(material_index).unwrap();
            let substance_index = material.get_substance_index().unwrap();
            let substance = building.get_substance(substance_index).unwrap();

            // Calculate the optimum_dx
            let thickness = material.thickness().unwrap();
            let alpha = substance.thermal_diffusivity().unwrap();
            let optimum_dx = (safety * 2.0 * alpha * dt).sqrt();

            // make all the elements of equal thickness
            let m = (thickness / optimum_dx).floor();

            let dx: f64;
            if m == 0. {
                // The given timestep allows for simulating the whole
                // layer as a single element (i.e. optimum_dx > dx).
                // However, dx cannot be greater than thickness, so limit it.
                dx = thickness;
            } else {
                // Layer must be subdivided in m equal parts.
                dx = thickness / m;
            }

            if dx > max_dx {
                // If the found dx is larger than the max allowed d_x, try to change timestep
                // check if there is room for reducing dt...
                let next_dt = main_dt / ((n + 1) as f64);
                if next_dt > min_dt {
                    // If there is room for that, do it.
                    return aux(building, c, main_dt, n + 1, max_dx, min_dt);
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

        // return
        (n, n_elements)
    }

    aux(building, c, main_dt, 1, max_dx, min_dt)
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
        return Err("Wrong discretization scheme... no elements considered".to_string());
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

/// Constructions are assumed to be a sandwich where zero or more massive
/// layers are located between two non-mass layers. These non-mass layers will always
/// include the interior and exterior film convections coefficients, respectively (which is
/// why they are called `full_rsi` and `full_rso`, respectively). Additionally,
/// they can also include some lightweight insulation or any other material of negligible
/// thermal mass.
///
/// The nodes are placed in between elements. Each element can be represented by the
/// 2x2 matrix
///
/// ```math
///
/// K=\begin{bmatrix} -1/R & 1/R \\
/// 1/R & -1/R
/// \end{bmatrix}   ;   R=\frac{thickness}{\lambda}
///
///```
///
///
/// The equation to solve is the following:
///
/// ```math
/// \overline{C}  \dot{T} + \overline{K}  T = q
/// ```
///
/// Where $`\overline{C}`$ and $`\overline{K}`$ are (more or less) matrices representing the
/// thermal mass and resistance of each node, respectively; and where `T` and `q` are
/// vectors representing the Temperature and the heat flow into each node, respectively.
/// $`\overline{C}`$ and $`\overline{K}`$ are the result of finite difference method.
///
/// This model also uses finite differences to march through time. When doing this,
/// the previous equation can be represented as follows:
///
/// ```math
/// \overline{C}  \frac{ (T_{i+1}- T_{i}) } {\Delta t} - \overline{K} T_i = q
/// ```
///
/// Hence, the temperatures of one step ($`T_{i+1}`$) can be calculated based on
/// the temperatures of the step before $`T_{i}`$ as follows.
///
/// ```math
/// {T_{i+1}}  =   T_i +  \Delta t  \overline{C}^{-1}\overline{K} T_i + \Delta t\overline{C}^{-1}q
/// ```
/// Where $`I_n`$ is the Identity Matrix. Based on that equation, it is convenient to define $`\overline{K}'`$ as follows:
///
/// ```math
/// \overline{K}' = \Delta t  \overline{C}^{-1}\overline{K}
/// ```
///
/// The value of $`\overline{K}`$ is never stored, as only $`\overline{K}'`$ is really needed for marching
/// through time.
///
/// Also, note that—unless some layer of the Surface is generating heat—all the elements of $`q`$ are Zero
/// exept the first one and the last one, which are $`\frac{T_{in}}{R_{si}}`$ and $`\frac{T_{out}}{R_{so}}`$,
/// respectively. These values are calculated when marching forward in time. However,
/// the value of $`\Delta t\overline{C}^{-1}`$ that accompanies $`q`$ in the previous
/// equation are stored in the values `c_i` and `c_o`, respectively.
///
///
/// # Aguments
/// * building: the Building containing the construction and all other data.
/// * construction: the construction being discretized
/// * dt: the timestep for the model (relevant for pre-processing some elements that will be used when marching forward in time)
/// * n_elements: The number of elements in each layer of the construction
/// * rsi: the interior film coefficient
/// * rso: the exterior film coefficient
///
/// # Outputs
///
/// * k_prime: the Matrix representing the heat transfer between the nodes... will be filled by this function, and the input should be full of zeros
/// * full_rsi: R_si + the thermal resistance of all the no-mass layers before the first massive element.
/// * full_rso: R_so + the thermal resistance of all the no-mass layers after the last massive element.
/// * c_i: The thermal mass of the most interior node, divided by the timestep
/// * c_o: The thermal mass of the most exterior node, divided by the timestep
pub fn build_thermal_network(
    building: &Building,
    c: &Construction,
    dt: f64,
    n_elements: &[usize],
    rs_i: f64,
    rs_o: f64,
    k_prime: &mut Matrix,
    full_rsi: &mut f64,
    full_rso: &mut f64,
    c_i: &mut f64,
    c_o: &mut f64,
) -> Result<(), String> {
    // check coherence in input data
    if n_elements.len() != c.n_layers() {
        let err = format!("Mismatch between number of layers in construction ({}) and the number of elements in scheme ({})",c.n_layers(),n_elements.len());
        return Err(err);
    }

    // Calculate number of nodes
    let all_nodes = calc_n_total_nodes(&n_elements)?;

    // Check the size of k_prime
    let (rows, cols) = k_prime.size();
    if rows != all_nodes || cols != all_nodes {
        let err = format!(
            "Unexpected size of given matrix - found ({},{}) and was expecting ({},{})",
            rows, cols, all_nodes, all_nodes
        );
        return Err(err);
    }

    #[cfg(debug_assertions)]
    {
        let (nrows, ncols) = k_prime.size();
        for row in 0..nrows {
            for col in 0..ncols {
                debug_assert!(k_prime.get(row, col).unwrap().abs() < f64::EPSILON);
            }
        }
    }

    // NOW, PROCESS
    ////////////////
    let (first_massive, last_massive) = get_first_and_last_massive_elements(n_elements)?;

    if first_massive == 0 && last_massive == 0 {
        // no massive layers at all in construction.
        // Simple case... return Zero mass and an R value of 1/R
        let r = r_value(building, c).unwrap() + rs_i + rs_o;
        *full_rsi = r;
        *full_rso = r;
        k_prime.set(0, 0, -1.0 / r).unwrap();
        k_prime.set(0, 1, 1.0 / r).unwrap();
        k_prime.set(1, 0, 1.0 / r).unwrap();
        k_prime.set(1, 1, -1.0 / r).unwrap();
        *c_i = 0.0;
        *c_o = 0.0;
        return Ok(());
    } else {
        // Calculate inner and outer surface Resistances

        // Everything before the first massive layer
        *full_rsi = rs_i;
        for i in 0..first_massive {
            let material_index = c.get_layer_index(i).unwrap();
            let material = building.get_material(material_index).unwrap();
            let substance_index = material.get_substance_index().unwrap();
            let substance = building.get_substance(substance_index).unwrap();

            *full_rsi += material.thickness().unwrap() / substance.thermal_conductivity().unwrap();
        }

        // Everything after the last massive layer
        *full_rso = rs_o;
        for i in last_massive..c.n_layers() {
            let material_index = c.get_layer_index(i).unwrap();
            let material = building.get_material(material_index).unwrap();
            let substance_index = material.get_substance_index().unwrap();
            let substance = building.get_substance(substance_index).unwrap();
            *full_rso += material.thickness().unwrap() / substance.thermal_conductivity().unwrap();
        }
    }

    // Calculate what is in between the massive layers
    let mut node: usize = 0;
    let mut n_layer: usize = first_massive;

    while n_layer < last_massive {
        let material_index = c.get_layer_index(n_layer).unwrap();
        let material = building.get_material(material_index).unwrap();

        let m = n_elements[n_layer];
        if m == 0 {
            // nomass material
            // add up all the R of the no-mass layers that
            // are together
            let mut r = 0.0; // if the material is no mass, then the first value
            while n_layer < last_massive && n_elements[n_layer] == 0 {
                let material_index = c.get_layer_index(n_layer).unwrap();
                let material = building.get_material(material_index).unwrap();
                let dx = material.thickness().unwrap();

                let substance_index = material.get_substance_index().unwrap();
                let substance = building.get_substance(substance_index).unwrap();
                let k = substance.thermal_conductivity().unwrap();

                r += dx / k;
                n_layer += 1;
            }

            // update values
            let u = 1. / r;
            // top left
            let old_value = k_prime.get(node, node).unwrap();
            k_prime.set(node, node, old_value - u).unwrap();
            // top right
            let old_value = k_prime.get(node, node + 1).unwrap();
            k_prime.set(node, node + 1, old_value + u).unwrap();
            // bottom left
            let old_value = k_prime.get(node + 1, node).unwrap();
            k_prime.set(node + 1, node, old_value + u).unwrap();
            // bottom right
            let old_value = k_prime.get(node + 1, node + 1).unwrap();
            k_prime.set(node + 1, node + 1, old_value - u).unwrap();

            // Move one node ahead
            node += 1;
        } else {
            // calc U value
            let substance_index = material.get_substance_index().unwrap();
            let substance = building.get_substance(substance_index).unwrap();

            let k = substance.thermal_conductivity().unwrap();
            let dx = material.thickness().unwrap() / (m as f64);
            let u: f64 = k / dx;

            for _ in 0..m {
                // top left
                let old_value = k_prime.get(node, node).unwrap();
                k_prime.set(node, node, old_value - u).unwrap();
                // top right
                let old_value = k_prime.get(node, node + 1).unwrap();
                k_prime.set(node, node + 1, old_value + u).unwrap();
                // bottom left
                let old_value = k_prime.get(node + 1, node).unwrap();
                k_prime.set(node + 1, node, old_value + u).unwrap();
                // bottom right
                let old_value = k_prime.get(node + 1, node + 1).unwrap();
                k_prime.set(node + 1, node + 1, old_value - u).unwrap();

                // advance node.
                node += 1;
            }
            n_layer += 1;
        }
    }

    // ADD RSI AND RSO
    // add r_si to first node
    if *full_rsi > 0.0 {
        let old_value = k_prime.get(0, 0).unwrap();
        k_prime.set(0, 0, old_value - 1.0 / (*full_rsi)).unwrap();
    }

    if *full_rso > 0.0 {
        let old_value = k_prime.get(all_nodes - 1, all_nodes - 1).unwrap();
        k_prime
            .set(all_nodes - 1, all_nodes - 1, old_value - 1.0 / (*full_rso))
            .unwrap();
    }

    // CALCULATE MASSES
    let mut left_side: Vec<f64> = vec![0.0; all_nodes];
    node = 0;

    for n_layer in 0..c.n_layers() {
        let layer_index = c.get_layer_index(n_layer).unwrap();
        let material = building.get_material(layer_index).unwrap();

        let m = n_elements[n_layer];

        if m != 0 {
            // if has mass
            for _ in 0..m {
                // Calc mass
                let substance_index = material.get_substance_index().unwrap();
                let substance = building.get_substance(substance_index).unwrap();

                let rho = substance.density().unwrap();
                let cp = substance.specific_heat_capacity().unwrap();
                let dx = material.thickness().unwrap() / (m as f64);
                let m = rho * cp * dx / dt;

                // Nodes are in the joints between
                // finite elements... add half mass from each
                left_side[node] += m / 2.0;
                left_side[node + 1] += m / 2.0;

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

    // DIVIDE ONE BY THE OTHER
    for (i, mass) in left_side.iter().enumerate() {
        if *mass != 0.0 {
            // Multiply the whole column by this.
            for j in 0..all_nodes {
                let old_k_value = k_prime.get(i, j).unwrap();
                k_prime.set(i, j, old_k_value / mass).unwrap();
            }
        } // Else, masses are already Zero.
    }

    // SET C_I AND C_O
    *c_i = left_side[0];
    *c_o = left_side[all_nodes - 1];

    Ok(())
}
