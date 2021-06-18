#![allow(clippy::too_many_arguments)]

use matrix::Matrix;

use building_model::boundary::Boundary;
use building_model::building::Building;
use building_model::fenestration::Fenestration;
use building_model::object_trait::ObjectTrait;
use building_model::surface::Surface;
use convection::*;
use simulation_state::simulation_state::SimulationState;
use simulation_state::simulation_state_element::SimulationStateElement;

use crate::construction::*;

/// This is a Surface from the point of view of our thermal solver.
/// Since this module only calculate heat transfer (and not short-wave solar
/// radiation, e.g., light), both building_model::Fenestration and building_model::Surface
/// are treated in the same way.
pub struct ThermalSurface {
    /// The index of this surface within
    /// the Thermal Model surfaces array
    index: usize,

    // A reference to the original Surface in the
    // building model
    surface_index: usize,

    /// The front side convection coefficient
    pub rs_front: f64,

    /// The back side convection coefficient
    pub rs_back: f64,

    /// The interior (i.e. front side) resistance before
    /// any layer with mass. It includes the r_si and also
    /// any light-weight material at the front of the construction.
    /// If the first layer in the construction has mass, then this
    /// value will be equal to r_si
    pub full_rs_front: f64,

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
    pub full_rs_back: f64,

    /// The matrix that represents the thermal network
    k_prime: matrix::Matrix,

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
    front_boundary: Boundary,

    /// The location of the back boundary zone in the
    /// Zones array of the Thermal Model    
    back_boundary: Boundary,

    /// The area of the Surface
    area: f64,

    /// Is this a Fenestration or a Surface in the original Building model?
    pub is_fenestration: bool,
}

impl ThermalSurface {
    /// Constructs a new ThermalSurface object.
    ///
    /// # Parameters
    /// * building: the Building object containing the Construction
    /// * state: the SimulationState element that contains the results
    /// * surface_index: The index of the Surface represented by this surface in the Building object
    /// * construction_index: The index of the Construction of this surface in the Building object
    /// * dt: the timestep for the model
    /// * area: area of the Surface,
    /// * rs_front: the interior film coefficient,
    /// * rs_back: the exterior film coefficient,
    /// * n_elements: A reference to a Vec<usize> containing the number of nodes in each layer of the construction
    /// * index: The index assigned to this surface within the Thermal Model surfaces array
    /// * is_fenestration: Is a fenestration?    
    fn new(
        building: &Building,
        state: &mut SimulationState,
        surface_index: usize,
        construction_index: usize,
        dt: f64,
        area: f64,
        rs_front: f64,
        rs_back: f64,
        n_elements: &[usize],
        index: usize,
        is_fenestration: bool,
    ) -> Result<Self, String> {
        // Borrow the construction from the Building object
        let construction = building.get_construction(construction_index).unwrap();

        // Calculate number of nodes in that construction.
        let n_nodes = calc_n_total_nodes(&n_elements)?;

        // Push elements to the SimulationState
        for i in 0..n_nodes {
            if is_fenestration {
                state.push(SimulationStateElement::FenestrationNodeTemperature(
                    surface_index,
                    i,
                    22.0,
                ));
            } else {
                state.push(SimulationStateElement::SurfaceNodeTemperature(
                    surface_index,
                    i,
                    22.0,
                ));
            }
        }

        // Build ThermalSurface with some placeholders
        let mut ret = ThermalSurface {
            surface_index,
            rs_front,
            rs_back,
            full_rs_back: 0.0,  // filled when building thermal network
            full_rs_front: 0.0, // filled when building thermal network
            c_o: 0.0,           // filled when building thermal network
            c_i: 0.0,           // filled when building thermal network
            k_prime: Matrix::new(0.0, n_nodes, n_nodes), // filled when building thermal network
            n_nodes,
            massive: true,                  // filled after building the thermal network
            front_boundary: Boundary::None, // filled when setting boundary
            back_boundary: Boundary::None,
            index,
            area,
            is_fenestration,
        };

        // Build the thermal network for this surface
        build_thermal_network(
            building,
            &construction,
            dt,
            &n_elements,
            rs_front,
            rs_back,
            &mut ret.k_prime,
            &mut ret.full_rs_front,
            &mut ret.full_rs_back,
            &mut ret.c_i,
            &mut ret.c_o,
        )
        .unwrap();
        ret.massive = !(ret.c_o == 0.0 && ret.c_i == 0.);

        // return
        Ok(ret)
    }

    pub fn from_surface(
        building: &Building,
        state: &mut SimulationState,
        surface: &Surface,
        dt: f64,
        n_elements: &[usize],
        index: usize,
    ) -> Result<Self, String> {
        let surface_index = ObjectTrait::index(surface);

        // Check if Surface is valid... or else
        ObjectTrait::is_full(surface).unwrap();

        // this (should not fail because the surface is valid)
        let construction_index = surface.get_construction_index().unwrap();

        let area = surface.area().unwrap(); // should not fail because surface is full

        let (rs_front, rs_back) = calc_convection_coefficients(surface);

        ThermalSurface::new(
            building,
            state,
            surface_index,
            construction_index,
            dt,
            area,
            rs_front,
            rs_back,
            n_elements,
            index,
            false, // not a fenestration
        )
    }

    pub fn from_fenestration(
        building: &Building,
        state: &mut SimulationState,
        fenestration: &Fenestration,
        dt: f64,
        n_elements: &[usize],
        index: usize,
    ) -> Result<Self, String> {
        //let fenestration = building.get_fenestration(fenestration_index);
        let fenestration_index = ObjectTrait::index(fenestration);

        // Check if Surface is valid... or else
        ObjectTrait::is_full(fenestration).unwrap();

        // this (should not fail because the surface is valid)
        let construction_index = fenestration.get_construction_index().unwrap();

        let area = fenestration.area().unwrap(); // should not fail because surface is full

        let (rs_front, rs_back) = calc_convection_coefficients_for_fenestration(fenestration);

        ThermalSurface::new(
            building,
            state,
            fenestration_index,
            construction_index,
            dt,
            area,
            rs_front,
            rs_back,
            n_elements,
            index,
            true, // it is a fenestration
        )
    }

    pub fn area(&self) -> f64 {
        self.area
    }

    pub fn rs_front(&self) -> f64 {
        self.rs_front
    }

    pub fn rs_back(&self) -> f64 {
        self.rs_back
    }

    /// Calculates the heat flow out of the layer, based
    /// on the inside and outside temperatures
    fn calc_heat_flow(
        &self,
        building: &Building,
        state: &SimulationState,
        t_front: f64,
        t_back: f64,
    ) -> (f64, f64) {
        // Positive is going out of the layer.
        let q_in;
        let q_out;
        let t_si = self.front_temperature(building, state);
        let t_so = self.back_temperature(building, state);
        if self.massive {
            q_in = (t_si - t_front) / self.full_rs_front;
            q_out = (t_so - t_back) / self.full_rs_back;
        } else {
            q_in = (t_si - t_front) / self.rs_front;
            q_out = (t_so - t_back) / self.rs_back;
        }
        // return
        (q_in, q_out)
    }

    /// Checks whether a wall has thermal mass
    pub fn is_massive(&self) -> bool {
        self.massive
    }

    fn surface_front_temperature(&self, building: &Building, state: &SimulationState) -> f64 {
        let surf = building.get_surface(self.surface_index).unwrap();

        let i = surf.get_first_node_temperature_index().unwrap();

        if let SimulationStateElement::SurfaceNodeTemperature(surf_index, node_index, temperature) =
            state[i]
        {
            if surf_index != self.index {
                panic!(
                    "Incorrect index allocated for Temperature of SurfaceNode of Surface '{}'",
                    self.index
                );
            }
            if node_index != 0 {
                panic!(
                    "Incorrect index allocated for Front Temperature of of Surface '{}'",
                    self.index
                );
            }
            // all Good here... return
            temperature
        } else {
            panic!("Incorrect StateElement kind allocated for Temperature of SurfaceNode of Surface '{}'", self.index);
        }
    }

    fn fenestration_front_temperature(&self, building: &Building, state: &SimulationState) -> f64 {
        let surf = building.get_fenestration(self.surface_index).unwrap();

        let i = surf.get_first_node_temperature_index().unwrap();

        if let SimulationStateElement::FenestrationNodeTemperature(
            surf_index,
            node_index,
            temperature,
        ) = state[i]
        {
            if surf_index != self.index {
                panic!(
                    "Incorrect index allocated for Temperature of FenestrationNodeTemperature of Surface '{}'",
                    self.index
                );
            }
            if node_index != 0 {
                panic!(
                    "Incorrect index allocated for Front Temperature of of Fenestration '{}'",
                    self.index
                );
            }
            // all Good here... return
            temperature
        } else {
            panic!("Incorrect StateElement kind allocated for Temperature of FenestrationNodeTemperature of Fenestration '{}'", self.index);
        }
    }

    /// Gets the Front temperature
    pub fn front_temperature(&self, building: &Building, state: &SimulationState) -> f64 {
        if self.is_fenestration {
            self.fenestration_front_temperature(building, state)
        } else {
            self.surface_front_temperature(building, state)
        }
    }

    /// Gets the Back temperature
    pub fn back_temperature(&self, building: &Building, state: &SimulationState) -> f64 {
        if self.is_fenestration {
            self.fenestration_back_temperature(building, state)
        } else {
            self.surface_back_temperature(building, state)
        }
    }

    /// Gets the Back temperature of a surface
    fn surface_back_temperature(&self, building: &Building, state: &SimulationState) -> f64 {
        let surf = building.get_surface(self.surface_index).unwrap();

        let i = surf.get_last_node_temperature_index().unwrap();

        if let SimulationStateElement::SurfaceNodeTemperature(surf_index, node_index, temperature) =
            state[i]
        {
            if surf_index != self.index {
                panic!(
                    "Incorrect index allocated for Back Temperature of SurfaceNode of Surface '{}'",
                    self.index
                );
            }
            if node_index != self.n_nodes - 1 {
                panic!(
                    "Incorrect index allocated for Back Temperature of Surface '{}' ({})... expected {}, found {}",
                    self.index,
                    if self.is_fenestration {"fenestration"}else{"non-fenestration"},
                    self.n_nodes - 1,
                    node_index
                );
            }
            // all Good here... return
            temperature
        } else {
            panic!("Incorrect StateElement kind allocated for Temperature of SurfaceNode of Surface '{}'", self.index);
        }
    }

    /// Gets the Back temperature of a surface
    fn fenestration_back_temperature(&self, building: &Building, state: &SimulationState) -> f64 {
        let surf = building.get_fenestration(self.surface_index).unwrap();

        let i = surf.get_last_node_temperature_index().unwrap();

        if let SimulationStateElement::FenestrationNodeTemperature(
            surf_index,
            node_index,
            temperature,
        ) = state[i]
        {
            if surf_index != self.index {
                panic!(
                    "Incorrect index allocated for Back Temperature of FenestrationNodeTemperature of Fenestration '{}'",
                    self.index
                );
            }
            if node_index != self.n_nodes - 1 {
                panic!(
                    "Incorrect index allocated for Back Temperature of Fenestration '{}' ({})... expected {}, found {}",
                    self.index,
                    if self.is_fenestration {"fenestration"}else{"non-fenestration"},
                    self.n_nodes - 1,
                    node_index
                );
            }
            // all Good here... return
            temperature
        } else {
            panic!("Incorrect StateElement kind allocated for Temperature of FenestrationNodeTemperature of Fenestration '{}'", self.index);
        }
    }

    /// Sets the node temperatures in the state
    fn set_node_temperatures(
        &self,
        building: &Building,
        state: &mut SimulationState,
        matrix: &Matrix,
    ) {
        if self.is_fenestration {
            self.set_fenestration_node_temperatures(building, state, matrix)
        } else {
            self.set_surface_node_temperatures(building, state, matrix)
        }
    }

    fn set_surface_node_temperatures(
        &self,
        building: &Building,
        state: &mut SimulationState,
        matrix: &Matrix,
    ) {
        let surf = building.get_surface(self.surface_index).unwrap();
        let ini = surf.get_first_node_temperature_index().unwrap();
        let fin = ini + self.n_nodes;

        for i in ini..fin {
            if let SimulationStateElement::SurfaceNodeTemperature(surf_index, node_index, _) =
                state[i]
            {
                if surf_index != self.index {
                    panic!(
                        "Incorrect index allocated for Temperature of SurfaceNode of Surface '{}'",
                        self.index
                    );
                }
                // all Good here
                let new_t = matrix.get(node_index, 0).unwrap();
                state[i] =
                    SimulationStateElement::SurfaceNodeTemperature(surf_index, node_index, new_t);
            } else {
                panic!("Incorrect StateElement kind allocated for Temperature of SurfaceNode of Surface '{}'", self.index);
            }
        }
    }

    fn set_fenestration_node_temperatures(
        &self,
        building: &Building,
        state: &mut SimulationState,
        matrix: &Matrix,
    ) {
        let surf = building.get_fenestration(self.surface_index).unwrap();
        let ini = surf.get_first_node_temperature_index().unwrap();
        let fin = ini + self.n_nodes;

        for i in ini..fin {
            if let SimulationStateElement::FenestrationNodeTemperature(surf_index, node_index, _) =
                state[i]
            {
                if surf_index != self.index {
                    panic!(
                        "Incorrect index allocated for Temperature of SurfaceNode of Surface '{}'",
                        self.index
                    );
                }
                // all Good here
                let new_t = matrix.get(node_index, 0).unwrap();
                state[i] = SimulationStateElement::FenestrationNodeTemperature(
                    surf_index, node_index, new_t,
                );
            } else {
                panic!("Incorrect StateElement kind allocated for Temperature of SurfaceNode of Surface '{}'", self.index);
            }
        }
    }

    /// Retrieves the state of the Surface as a Matrix
    /// object.
    pub fn get_node_temperatures(&self, building: &Building, state: &SimulationState) -> Matrix {
        if self.is_fenestration {
            self.get_fenestration_node_temperatures(building, state)
        } else {
            self.get_surface_node_temperatures(building, state)
        }
    }

    fn get_fenestration_node_temperatures(
        &self,
        building: &Building,
        state: &SimulationState,
    ) -> Matrix {
        let mut ret = Matrix::new(0.0, self.n_nodes, 1);
        let surf = building.get_fenestration(self.surface_index).unwrap();
        let ini = surf.get_first_node_temperature_index().unwrap();
        let fin = ini + self.n_nodes;
        for i in ini..fin {
            if let SimulationStateElement::FenestrationNodeTemperature(
                surf_index,
                node_index,
                temperature,
            ) = state[i]
            {
                if surf_index != self.index {
                    panic!("Incorrect index allocated for Temperature of FenestrationNode... surf_index = {}, self.index = {}", surf_index, self.index);
                }
                // all Good here
                ret.set(node_index, 0, temperature).unwrap();
            } else {
                panic!("Incorrect StateElement kind allocated for Temperature of FenestrationNode of Fenestration '{}'... found '{}'", self.index, state[i].to_string());
            }
        }
        ret
    }

    fn get_surface_node_temperatures(
        &self,
        building: &Building,
        state: &SimulationState,
    ) -> Matrix {
        let mut ret = Matrix::new(0.0, self.n_nodes, 1);
        let surf = building.get_surface(self.surface_index).unwrap();
        let ini = surf.get_first_node_temperature_index().unwrap();
        let fin = ini + self.n_nodes;
        for i in ini..fin {
            if let SimulationStateElement::SurfaceNodeTemperature(
                surf_index,
                node_index,
                temperature,
            ) = state[i]
            {
                if surf_index != self.index {
                    panic!("Incorrect index allocated for Temperature of SurfaceNode... surf_index = {}, self.index = {}", surf_index, self.index);
                }
                // all Good here
                ret.set(node_index, 0, temperature).unwrap();
            } else {
                panic!("Incorrect StateElement kind allocated for Temperature of SurfaceNode of Surface '{}'... found '{}'", self.index, state[i].to_string());
            }
        }
        ret
    }

    /// Marches one timestep. Returns front and back heat flow    
    pub fn march(
        &self,
        building: &Building,
        state: &mut SimulationState,
        t_front: f64,
        t_back: f64,
    ) -> (f64, f64) {
        let mut temperatures = self.get_node_temperatures(building, state);

        // println!("T front = {} | T back = {} ", t_front, t_back);
        // println!(" OLD TEMPS = {}", temperatures);

        if self.massive {
            // Update temperatures... T_i+1 = Ti + K_prime*Ti + {t_front/full_rs_front/C_i ... 0,0,0... t_back/full_rs_back/C_o}
            //                                     ^^^^^^^^^                    ^^^^^^^^^^^^
            //                        lets call this vector 'a'         These are the F components
            let mut a = self.k_prime.from_prod_n_diag(&temperatures, 3).unwrap(); // k_prime is tri-diagonal
                                                                                  // ... 'a' should be a vector

            // println!(" A_before = {}", a);

            // Let's add the F components
            let old_value = a.get(0, 0).unwrap();
            a.set(0, 0, old_value + t_front / self.full_rs_front / self.c_i)
                .unwrap();

            let old_value = a.get(self.n_nodes - 1, 0).unwrap();
            a.set(
                self.n_nodes - 1,
                0,
                old_value + t_back / self.full_rs_back / self.c_o,
            )
            .unwrap();

            //println!(" A_after = {}", a);

            // Let's add a to the temperatures.
            temperatures.add_to_this(&a).unwrap();
        } else {
            // full_rs_front is, indeed, the whole R
            let q = (t_back - t_front) / self.full_rs_front;

            let ts_front = t_front + q * self.rs_front;
            let ts_back = t_back - q * self.rs_back;

            temperatures.set(0, 0, ts_front).unwrap();
            temperatures.set(self.n_nodes - 1, 0, ts_back).unwrap();
        }
        // println!(" Temperatures_after = {}", temperatures);
        // Set state
        self.set_node_temperatures(building, state, &temperatures);

        // return        
        self.calc_heat_flow(building, state, t_front, t_back)
    }

    pub fn set_front_boundary(&mut self, b: Boundary) {
        self.front_boundary = b;
    }

    pub fn set_back_boundary(&mut self, b: Boundary) {
        self.back_boundary = b;
    }

    pub fn front_boundary(&self) -> Boundary {
        self.front_boundary
    }

    pub fn back_boundary(&self) -> Boundary {
        self.back_boundary
    }
} // END OF SURFACE IMPL

/***********/
/* TESTING */
/***********/

#[cfg(test)]
mod testing {
    use super::*;
    use building_model::material::MaterialProperties;
    use building_model::substance::SubstanceProperties;
    use geometry3d::loop3d::Loop3D;
    use geometry3d::point3d::Point3D;
    use geometry3d::polygon3d::Polygon3D;

    fn add_polyurethane(building: &mut Building) -> usize {
        let poly_index = building.add_substance("polyurethane".to_string());
        building
            .set_substance_properties(
                poly_index,
                SubstanceProperties {
                    density: 17.5,                 // kg/m3... reverse engineered from paper
                    specific_heat_capacity: 2400., // J/kg.K
                    thermal_conductivity: 0.0252,  // W/m.K
                },
            )
            .unwrap();
        {
            let poly = building.get_substance(poly_index).unwrap();
            assert_eq!(poly.thermal_diffusivity().unwrap(), 0.6E-6);
        }

        poly_index
    }

    fn add_brickwork(building: &mut Building) -> usize {
        let brickwork_index = building.add_substance("brickwork".to_string());

        building
            .set_substance_properties(
                brickwork_index,
                SubstanceProperties {
                    density: 1700.,               // kg/m3... reverse engineered from paper
                    specific_heat_capacity: 800., // J/kg.K
                    thermal_conductivity: 0.816,  // W/m.K
                },
            )
            .unwrap();
        {
            let brickwork = building.get_substance(brickwork_index).unwrap();
            assert_eq!(brickwork.thermal_diffusivity().unwrap(), 0.6E-6);
        }
        brickwork_index
    }

    fn add_layer(building: &mut Building, substance_index: usize, thickness: f64) -> usize {
        let m0_index = building.add_material("123123".to_string());
        building
            .set_material_substance(m0_index, substance_index)
            .unwrap();
        building
            .set_material_properties(
                m0_index,
                MaterialProperties {
                    thickness: thickness,
                },
            )
            .unwrap();

        m0_index
    }

    fn get_wall_1() -> Building {
        let mut building = Building::new("The Building".to_string());

        /* SUBSTANCES */

        // Add polyurethane Substance
        let poly_index = add_polyurethane(&mut building);

        /* MATERIAL */
        let m0_thickness = 200.0 / 1000. as f64;
        let m0_index = add_layer(&mut building, poly_index, m0_thickness);

        /* WALL 1 */
        let c0_index = building.add_construction("Wall 1".to_string());
        building
            .add_material_to_construction(c0_index, m0_index)
            .unwrap();

        return building;
    }

    fn get_wall_2() -> Building {
        let mut building = Building::new("The Building".to_string());

        /* SUBSTANCES */
        // Add polyurethane Substance
        let poly_index = add_polyurethane(&mut building);

        // Add brickwork
        let brickwork_index = add_brickwork(&mut building);

        /* MATERIALS */
        let m0_thickness = 200.0 / 1000. as f64;
        let m0_index = add_layer(&mut building, poly_index, m0_thickness);

        let m1_thickness = 110.0 / 1000. as f64;
        let m1_index = add_layer(&mut building, brickwork_index, m1_thickness);

        /* WALL 2 */

        let c0_index = building.add_construction("Wall 2".to_string());
        building
            .add_material_to_construction(c0_index, m0_index)
            .unwrap();
        building
            .add_material_to_construction(c0_index, m1_index)
            .unwrap();

        return building;
    }

    fn get_wall_3() -> Building {
        let mut building = Building::new("The Building".to_string());

        /* SUBSTANCES */
        // Add polyurethane Substance
        let poly_index = add_polyurethane(&mut building);

        // Add brickwork
        let brickwork_index = add_brickwork(&mut building);

        /* MATERIALS */
        let poly_mat_thickness = 20.0 / 1000. as f64;
        let poly_mat_index = add_layer(&mut building, poly_index, poly_mat_thickness);

        let brickwork_mat_thickness = 220.0 / 1000. as f64;
        let brickwork_mat_index =
            add_layer(&mut building, brickwork_index, brickwork_mat_thickness);

        /* WALL 3 */

        let c0_index = building.add_construction("Wall 3".to_string());
        building
            .add_material_to_construction(c0_index, poly_mat_index)
            .unwrap();
        building
            .add_material_to_construction(c0_index, brickwork_mat_index)
            .unwrap();
        building
            .add_material_to_construction(c0_index, poly_mat_index)
            .unwrap();

        return building;
    }

    #[test]
    fn test_new() {
        let mut building = Building::new("A building".to_string());

        /* SUBSTANCES */
        // Add polyurethane Substance
        let poly_index = add_polyurethane(&mut building);

        let m0_thickness = 200.0 / 1000. as f64;
        let m0_index = add_layer(&mut building, poly_index, m0_thickness);

        /* CONSTRUCTION */
        let c0_index = building.add_construction("Wall".to_string());
        building
            .add_material_to_construction(c0_index, m0_index)
            .unwrap();

        /* GEOMETRY */
        let mut the_loop = Loop3D::new();
        let l = 1. as f64;
        the_loop.push(Point3D::new(-l, -l, 0.)).unwrap();
        the_loop.push(Point3D::new(l, -l, 0.)).unwrap();
        the_loop.push(Point3D::new(l, l, 0.)).unwrap();
        the_loop.push(Point3D::new(-l, l, 0.)).unwrap();
        the_loop.close().unwrap();
        let p = Polygon3D::new(the_loop).unwrap();

        /* SURFACE */
        let s_index = building.add_surface("Surface".to_string());
        building.set_surface_polygon(s_index, p).unwrap();
        building
            .set_surface_construction(s_index, c0_index)
            .unwrap();

        let m0 = building.get_material(m0_index).unwrap();
        let c0 = building.get_construction(c0_index).unwrap();
        let surface = building.get_surface(s_index).unwrap();

        // This should result in a dt of 300, with 8 layers of 0.25m thick.
        let main_dt = 300.0;
        let max_dx = m0.thickness().unwrap() / 4.0;
        let min_dt = 1.;
        let (n_subdivisions, nodes) =
            discretize_construction(&building, &c0, main_dt, max_dx, min_dt);
        let dt = main_dt / n_subdivisions as f64;
        let mut state: SimulationState = SimulationState::new();
        let ts =
            ThermalSurface::from_surface(&building, &mut state, &surface, dt, &nodes, 0).unwrap();

        let (rs_front, rs_back) = calc_convection_coefficients(&surface);
        assert!(ts.massive);
        assert_eq!(ts.n_nodes, 9);
        assert_eq!(ts.rs_front, rs_front);
        assert_eq!(ts.rs_back, rs_back);
        assert_eq!(ts.area, 4.0);
    }

    #[test]
    fn test_get_dt_and_n_nodes_wall_1() {
        let building = get_wall_1();

        let c0_index = 0;
        let c0 = building.get_construction(c0_index).unwrap();

        let m0_index = 0;
        let m0 = building.get_material(m0_index).unwrap();
        let m0_thickness = m0.thickness().unwrap();

        // This should result in a dt of 300, with 8 layers of 0.25m thick.
        let main_dt = 300.;
        let (n_subdivisions, n_nodes) =
            discretize_construction(&building, &c0, main_dt, m0_thickness / 4., 1.);
        let dt = main_dt / n_subdivisions as f64;

        assert_eq!(n_nodes.len(), c0.n_layers());
        assert_eq!(dt, 300.);
        assert_eq!(n_nodes[0], 8);

        // This should result in a dt of 300, with 8 layers of 0.25m thick.
        let main_dt = 300.;
        let (n_subdivisions, n_nodes) =
            discretize_construction(&building, &c0, main_dt, m0_thickness / 8., 1.);
        let dt = main_dt / n_subdivisions as f64;
        assert_eq!(n_nodes.len(), c0.n_layers());
        assert_eq!(dt, 300.);
        assert_eq!(n_nodes[0], 8);

        // This should result in a dt of 300/2=150, with 12 layers of 0.0166667m thick.
        let main_dt = 300.;
        let (n_subdivisions, n_nodes) =
            discretize_construction(&building, &c0, main_dt, m0_thickness / 9., 1.);
        let dt = main_dt / n_subdivisions as f64;
        assert_eq!(n_nodes.len(), c0.n_layers());
        assert_eq!(dt, 150.);
        assert_eq!(n_nodes[0], 12);

        // This should result in a dt of 300/4=75, with 17 layers of 0.011...m thick.
        let main_dt = 300.;
        let (n_subdivisions, n_nodes) =
            discretize_construction(&building, &c0, main_dt, m0_thickness / 15., 1.);
        let dt = main_dt / n_subdivisions as f64;
        assert_eq!(n_nodes.len(), c0.n_layers());
        assert_eq!(dt, 75.);
        assert_eq!(n_nodes[0], 17);

        // We imposed a min_dt of 80, meaning that this should result in a no-mass layer
        let main_dt = 300.;
        let (n_subdivisions, n_nodes) =
            discretize_construction(&building, &c0, main_dt, m0_thickness / 15., 80.);
        let dt = main_dt / n_subdivisions as f64;
        assert_eq!(n_nodes.len(), c0.n_layers());
        assert_eq!(dt, 100.);
        assert_eq!(n_nodes[0], 0);
    }

    #[test]
    fn test_get_dt_and_n_nodes_wall_2() {
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
        let (n_subdivisions, n_nodes) = discretize_construction(&building, &c0, main_dt, 0.03, 1.);
        let dt = main_dt / n_subdivisions as f64;
        assert_eq!(n_nodes.len(), 2);
        assert_eq!(dt, 300.);
        assert_eq!(n_nodes[0], 8);
        assert_eq!(n_nodes[1], 4);

        // This should result in a dt of 300/2=150, with 12 and 6 layers.
        let main_dt = 300.;
        let (n_subdivisions, n_nodes) = discretize_construction(&building, &c0, main_dt, 0.02, 1.);
        let dt = main_dt / n_subdivisions as f64;
        assert_eq!(n_nodes.len(), c0.n_layers());
        assert_eq!(dt, 150.);
        assert_eq!(n_nodes[0], 12);
        assert_eq!(n_nodes[1], 6);

        // This should result in a dt of 300/3=100, with 14 and 8
        let main_dt = 300.;
        let (n_subdivisions, n_nodes) = discretize_construction(&building, &c0, main_dt, 0.015, 1.);
        let dt = main_dt / n_subdivisions as f64;
        assert_eq!(n_nodes.len(), c0.n_layers());
        assert_eq!(dt, 100.);
        assert_eq!(n_nodes[0], 14);
        assert_eq!(n_nodes[1], 8);

        // We imposed a min_dt of 100, meaning that this should result in a no-mass layer
        //let main_dt = 300.;
        let (n_subdivisions, n_nodes) =
            discretize_construction(&building, &c0, main_dt, 0.015, 110.);
        let dt = main_dt / n_subdivisions as f64;
        assert_eq!(n_nodes.len(), c0.n_layers());
        assert_eq!(dt, 150.);
        assert_eq!(n_nodes[0], 0);
        assert_eq!(n_nodes[1], 0);
    }

    #[test]
    fn test_get_dt_and_n_nodes_wall_3() {
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
        let (n_subdivisions, n_nodes) = discretize_construction(&building, &c0, main_dt, 0.03, 1.);

        let dt = main_dt / n_subdivisions as f64;
        assert_eq!(n_nodes.len(), c0.n_layers());
        assert_eq!(dt, 300.);
        assert_eq!(n_nodes[0], 0);
        assert_eq!(n_nodes[1], 9);
        assert_eq!(n_nodes[2], 0);

        // This should result in a dt of 300/6=50, with [2,23,2] layers
        let main_dt = 300.;
        let (n_subdivisions, n_nodes) = discretize_construction(&building, &c0, main_dt, 0.015, 1.);
        let dt = main_dt / n_subdivisions as f64;
        assert_eq!(n_nodes.len(), c0.n_layers());
        assert_eq!(dt, 50.);
        assert_eq!(n_nodes[0], 2);
        assert_eq!(n_nodes[1], 23);
        assert_eq!(n_nodes[2], 2);

        // Limit min_time to 65... This should result in a dt of 300/4=75, with [0, 18, 0] layers
        let main_dt = 300.;
        let (n_subdivisions, n_nodes) =
            discretize_construction(&building, &c0, main_dt, 0.015, 65.);
        let dt = main_dt / n_subdivisions as f64;
        assert_eq!(n_nodes.len(), c0.n_layers());
        assert_eq!(dt, 75.);
        assert_eq!(n_nodes[0], 0);
        assert_eq!(n_nodes[1], 18);
        assert_eq!(n_nodes[2], 0);
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
    fn test_wall_brickwork() {
        let mut building = Building::new("The building".to_string());

        // Add brickwork
        let brickwork_index = add_brickwork(&mut building);

        /* WALL 1: Single layer, with mass. */
        let brickwork_200_thickness = 200.0 / 1000. as f64;
        let brickwork_200_index =
            add_layer(&mut building, brickwork_index, brickwork_200_thickness);

        let wall_1_index = building.add_construction("Wall 1".to_string());
        building
            .add_material_to_construction(wall_1_index, brickwork_200_index)
            .unwrap();

        let wall_1 = building.get_construction(wall_1_index).unwrap();

        let dt = 156.0;
        let n_elements: Vec<usize> = vec![4];
        let all_nodes = calc_n_total_nodes(&n_elements).unwrap();

        let mut k_prime = Matrix::new(0.0, all_nodes, all_nodes);
        let mut full_rs_front = 0.0;
        let mut full_rs_back = 0.0;
        let mut c_i = 0.0;
        let mut c_o = 0.0;
        build_thermal_network(
            &building,
            &wall_1,
            dt,
            &n_elements,
            0., // rs_front
            0., // rs_back
            &mut k_prime,
            &mut full_rs_front,
            &mut full_rs_back,
            &mut c_i,
            &mut c_o,
        )
        .unwrap();

        let substance = building.get_substance(brickwork_index).unwrap();
        let rho = substance.density().unwrap();
        let cp = substance.specific_heat_capacity().unwrap();
        let k = substance.thermal_conductivity().unwrap();

        let n_layer = 0;
        let m = n_elements[n_layer];
        let dx = brickwork_200_thickness / (m as f64);
        let substance = building.get_substance(brickwork_index).unwrap();
        let alpha = substance.thermal_diffusivity().unwrap();
        let mass = rho * cp * dx / dt;

        // check coefficients
        assert_eq!(full_rs_front, 0.0); // 2.0*dt/(rho*cp*dx));
        assert_eq!(full_rs_back, 0.0); // 2.0*dt/(rho*cp*dx));
        assert_eq!(c_i, rho * cp * dx / (2.0 * dt));
        assert_eq!(c_o, rho * cp * dx / (2.0 * dt));

        // Check first node
        let node = 0;

        let found = k_prime.get(node, node).unwrap();
        let exp = -2.0 * alpha * dt / dx / dx;

        //assert!( (found-exp).abs() < 1E-10 );
        assert_eq!(exp, found);

        let found = k_prime.get(node, node + 1).unwrap();
        let exp = 2.0 * k / dx / mass;
        assert!((found - exp).abs() < 1E-10);

        // check middle nodes
        for node in 1..all_nodes - 1 {
            let found = k_prime.get(node, node).unwrap();
            let exp = -2.0 * alpha * dt / dx / dx;
            //assert!( (found-exp).abs() < 1E-10 );
            assert_eq!(found, exp);

            let found = k_prime.get(node, node + 1).unwrap();
            let exp = k / dx / mass;
            assert!((found - exp).abs() < 1E-10);

            let found = k_prime.get(node, node - 1).unwrap();
            let exp = k / dx / mass;
            assert!((found - exp).abs() < 1E-10);
        }
        // check end node
        let node = all_nodes - 1;

        let found = k_prime.get(node, node).unwrap();
        let exp = -2.0 * alpha * dt / dx / dx;
        assert!((found - exp).abs() < 1E-10);

        let found = k_prime.get(node, node - 1).unwrap();
        let exp = 2.0 * k / dx / mass;
        assert!((found - exp).abs() < 1E-10);
    }

    #[test]
    fn test_wall_2() {
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
        let n_nodes: Vec<usize> = vec![3, 3];
        let all_nodes = calc_n_total_nodes(&n_nodes).unwrap();

        let cp0 = m0_substance.specific_heat_capacity().unwrap();
        let rho0 = m0_substance.density().unwrap();
        let dx0 = m0.thickness().unwrap() / 3.;
        let mass0 = rho0 * cp0 * dx0 / dt;

        let cp1 = m1_substance.specific_heat_capacity().unwrap();
        let rho1 = m1_substance.density().unwrap();
        let dx1 = m1.thickness().unwrap() / 3.;
        let mass1 = rho1 * cp1 * dx1 / dt;

        let mut k_prime = Matrix::new(0.0, all_nodes, all_nodes);
        let mut full_rs_front = 0.0;
        let mut full_rs_back = 0.0;
        let mut c_i = 0.0;
        let mut c_o = 0.0;
        build_thermal_network(
            &building,
            &c0,
            dt,
            &n_nodes,
            0.,
            0.,
            &mut k_prime,
            &mut full_rs_front,
            &mut full_rs_back,
            &mut c_i,
            &mut c_o,
        )
        .unwrap();

        // check coefficients
        assert_eq!(full_rs_front, 0.0); // 2.0*dt/(rho*cp*dx));
        assert_eq!(full_rs_back, 0.0); // 2.0*dt/(rho*cp*dx));
        assert_eq!(c_i, rho0 * cp0 * dx0 / (2.0 * dt));
        assert_eq!(c_o, rho1 * cp1 * dx1 / (2.0 * dt));

        // FIRST LAYER
        //////////////
        let n_layer = 0;
        //let rho = m1.substance.density;
        //let cp = m1.substance.heat_capacity;
        let m = n_nodes[n_layer];
        let dx = m0.thickness().unwrap() / (m as f64);
        let k = m0_substance.thermal_conductivity().unwrap();
        let alpha = m0_substance.thermal_diffusivity().unwrap();

        // Check first node in layer
        let node = 0;

        let found = k_prime.get(node, node).unwrap();
        let exp = -2.0 * alpha * dt / dx / dx;
        assert!((found - exp).abs() < 1E-10);

        let found = k_prime.get(node, node + 1).unwrap();
        let exp = 2.0 * k / dx / mass0;
        assert!((found - exp).abs() < 1E-10);

        // check middle nodes
        for node in 1..3 {
            let found = k_prime.get(node, node).unwrap();
            let exp = -2.0 * alpha * dt / dx / dx;
            assert!((found - exp).abs() < 1E-10);

            let found = k_prime.get(node, node + 1).unwrap();
            let exp = k / dx / mass0;
            assert!((found - exp).abs() < 1E-10);

            let found = k_prime.get(node, node - 1).unwrap();
            let exp = k / dx / mass0;
            assert!((found - exp).abs() < 1E-10);
        }
        // check middle node (e.g. node 3), which
        // is also the beginning of Layer 2
        let node = 3;
        let mass_0 = m0_substance.density().unwrap()
            * m0_substance.specific_heat_capacity().unwrap()
            * m0.thickness().unwrap()
            / (6.0 * dt);
        let mass_1 = m1_substance.density().unwrap()
            * m1_substance.specific_heat_capacity().unwrap()
            * m1.thickness().unwrap()
            / (6.0 * dt);

        let found = k_prime.get(node, node).unwrap();
        let u0 = 3.0 * m0_substance.thermal_conductivity().unwrap() / m0.thickness().unwrap();
        let u1 = 3.0 * m1_substance.thermal_conductivity().unwrap() / m1.thickness().unwrap();
        let exp = -(u0 + u1) / (mass_0 + mass_1);
        assert_eq!(found, exp);

        let found = k_prime.get(node, node - 1).unwrap();
        assert_eq!(found, u0 / (mass0 / 2. + mass1 / 2.0));

        let found = k_prime.get(node, node + 1).unwrap();
        assert_eq!(found, u1 / (mass0 / 2.0 + mass1 / 2.0));

        // SECOND LAYER
        //////////////
        let n_layer = 1;
        //let rho = m2.substance.density;
        //let cp = m2.substance.heat_capacity;
        let m = n_nodes[n_layer];
        let dx = m1.thickness().unwrap() / (m as f64);
        let k = m1_substance.thermal_conductivity().unwrap();
        let alpha = m1_substance.thermal_diffusivity().unwrap();

        // check middle nodes in layer 2
        for node in 4..6 {
            let found = k_prime.get(node, node).unwrap();
            let exp = -2.0 * alpha * dt / dx / dx;
            assert!((found - exp).abs() < 1E-10);

            let found = k_prime.get(node, node + 1).unwrap();
            let exp = k / dx / mass1;
            assert!((found - exp).abs() < 1E-10);

            let found = k_prime.get(node, node - 1).unwrap();
            let exp = k / dx / mass1;
            assert!((found - exp).abs() < 1E-10);
        }
        // check end node
        let node = 6;

        let found = k_prime.get(node, node).unwrap();
        let u2 = 3.0 * m1_substance.thermal_conductivity().unwrap() / m1.thickness().unwrap();
        let exp = -2.0 * u2 / (mass1);
        assert_eq!(found, exp);

        let found = k_prime.get(node, node - 1).unwrap();
        assert_eq!(found, 2.0 * u1 / mass1);
    }

    #[test]
    fn test_wall_1() {
        let building = get_wall_1();

        let c0_index = 0;
        let c0 = building.get_construction(c0_index).unwrap();

        let m0_index = 0;
        let m0 = building.get_material(m0_index).unwrap();
        let m0_substance_index = m0.get_substance_index().unwrap();
        let m0_substance = building.get_substance(m0_substance_index).unwrap();

        let dt = 156.0;
        let n_nodes: Vec<usize> = vec![0];
        let all_nodes = calc_n_total_nodes(&n_nodes).unwrap();
        let mut k_prime = Matrix::new(0.0, all_nodes, all_nodes);

        let mut full_rs_front = 0.0;
        let mut full_rs_back = 0.0;
        let mut c_i = 0.0;
        let mut c_o = 0.0;
        build_thermal_network(
            &building,
            &c0,
            dt,
            &n_nodes,
            0.,
            0.,
            &mut k_prime,
            &mut full_rs_front,
            &mut full_rs_back,
            &mut c_i,
            &mut c_o,
        )
        .unwrap();

        // check coefficients
        assert_eq!(
            full_rs_front,
            m0.thickness().unwrap() / m0_substance.thermal_conductivity().unwrap()
        ); // 2.0*dt/(rho*cp*dx));
        assert_eq!(
            full_rs_back,
            m0.thickness().unwrap() / m0_substance.thermal_conductivity().unwrap()
        ); // 2.0*dt/(rho*cp*dx));
        assert_eq!(c_i, 0.0);
        assert_eq!(c_o, 0.0);

        for row in 0..all_nodes {
            for col in 0..all_nodes {
                let exp = m0_substance.thermal_conductivity().unwrap() / m0.thickness().unwrap();
                let found = k_prime.get(row, col).unwrap();
                if row == col {
                    assert_eq!(-exp, found);
                } else {
                    assert_eq!(exp, found);
                }
            }
        }
    }

    #[test]
    fn test_double_wall_1() {
        let mut building = get_wall_1();

        let c0_index = 0;
        let m0_index = 0;

        // add second layer
        building
            .add_material_to_construction(c0_index, m0_index)
            .unwrap();

        let c0 = building.get_construction(c0_index).unwrap();
        let m0 = building.get_material(m0_index).unwrap();

        let m0_substance_index = m0.get_substance_index().unwrap();
        let m0_substance = building.get_substance(m0_substance_index).unwrap();

        let dt = 156.0;
        let n_nodes: Vec<usize> = vec![0, 0];
        let all_nodes = calc_n_total_nodes(&n_nodes).unwrap();
        let r = 2.0 * m0.thickness().unwrap() / m0_substance.thermal_conductivity().unwrap();
        let mut k_prime = Matrix::new(0.0, all_nodes, all_nodes);
        let mut full_rs_front = 0.0;
        let mut full_rs_back = 0.0;
        let mut c_i = 0.0;
        let mut c_o = 0.0;
        build_thermal_network(
            &building,
            &c0,
            dt,
            &n_nodes,
            0.,
            0.,
            &mut k_prime,
            &mut full_rs_front,
            &mut full_rs_back,
            &mut c_i,
            &mut c_o,
        )
        .unwrap();

        // check coefficients
        assert_eq!(full_rs_front, r);
        assert_eq!(full_rs_back, r);
        assert_eq!(c_i, 0.);
        assert_eq!(c_o, 0.);

        for row in 0..all_nodes {
            for col in 0..all_nodes {
                let exp = 1.0 / r;
                let found = k_prime.get(row, col).unwrap();
                if row == col {
                    assert_eq!(-exp, found);
                } else {
                    assert_eq!(exp, found);
                }
            }
        }
    }

    #[test]
    fn test_wall_5() {
        let mut building = Building::new("A building".to_string());

        /* SUBSTANCES */
        let poly_index = add_polyurethane(&mut building);
        let brickwork_index = add_brickwork(&mut building);

        /* MATERIALS */
        let m0_thickness = 200.0 / 1000. as f64;
        let m0_index = add_layer(&mut building, poly_index, m0_thickness);

        let m1_thickness = 200.0 / 1000. as f64;
        let m1_index = add_layer(&mut building, brickwork_index, m1_thickness);

        /* WALL */

        let c0_index = building.add_construction("Wall 2".to_string());
        building
            .add_material_to_construction(c0_index, m0_index)
            .unwrap();
        building
            .add_material_to_construction(c0_index, m1_index)
            .unwrap();

        let m0 = building.get_material(m0_index).unwrap();
        let m0_substance_index = m0.get_substance_index().unwrap();
        let m0_substance = building.get_substance(m0_substance_index).unwrap();

        let m1 = building.get_material(m1_index).unwrap();
        let m1_substance_index = m1.get_substance_index().unwrap();
        let m1_substance = building.get_substance(m1_substance_index).unwrap();

        let c0 = building.get_construction(c0_index).unwrap();

        let dt = 156.0;
        let n_nodes: Vec<usize> = vec![1, 0];
        let all_nodes = calc_n_total_nodes(&n_nodes).unwrap();
        assert_eq!(all_nodes, 2);

        let cp = m0_substance.specific_heat_capacity().unwrap();
        let rho = m0_substance.density().unwrap();
        let dx = m0.thickness().unwrap();

        let mut k_prime = Matrix::new(0.0, all_nodes, all_nodes);
        let mut full_rs_front = 0.0;
        let mut full_rs_back = 0.0;
        let mut c_i = 0.0;
        let mut c_o = 0.0;
        build_thermal_network(
            &building,
            &c0,
            dt,
            &n_nodes,
            0.,
            0.,
            &mut k_prime,
            &mut full_rs_front,
            &mut full_rs_back,
            &mut c_i,
            &mut c_o,
        )
        .unwrap();

        // check coefficients
        assert_eq!(full_rs_front, 0.0); // 2.0*dt/(rho*cp*dx));
        assert_eq!(
            full_rs_back,
            m1.thickness().unwrap() / m1_substance.thermal_conductivity().unwrap()
        ); // 2.0*dt/(rho*cp*dx));
        assert_eq!(c_i, rho * cp * dx / (2.0 * dt));
        assert_eq!(c_o, rho * cp * dx / (2.0 * dt));

        let half_mass = m0_substance.density().unwrap()
            * m0_substance.specific_heat_capacity().unwrap()
            * m0.thickness().unwrap()
            / (2.0 * dt);
        let u = m0_substance.thermal_conductivity().unwrap() / m0.thickness().unwrap();

        // Check first node
        let found = k_prime.get(0, 0).unwrap();
        let exp = -u / half_mass;
        assert_eq!(exp, found);

        // Check last node
        let rso = m1.thickness().unwrap() / m1_substance.thermal_conductivity().unwrap();
        let found = k_prime.get(1, 1).unwrap();
        let exp = -(u + 1.0 / rso) / half_mass;
        assert!((exp - found).abs() < 1E-10);
    }

    #[test]
    fn test_wall_6() {
        let mut building = Building::new("A building".to_string());

        /* SUBSTANCES */

        let poly_index = add_polyurethane(&mut building);
        let brickwork_index = add_brickwork(&mut building);

        /* MATERIALS */

        let m0_thickness = 200.0 / 1000. as f64;
        let m0_index = add_layer(&mut building, brickwork_index, m0_thickness);

        let m1_thickness = 200.0 / 1000. as f64;
        let m1_index = add_layer(&mut building, poly_index, m1_thickness);

        /* WALL */

        let c0_index = building.add_construction("Wall".to_string());
        building
            .add_material_to_construction(c0_index, m1_index)
            .unwrap();
        building
            .add_material_to_construction(c0_index, m0_index)
            .unwrap();

        let c0 = building.get_construction(c0_index).unwrap();

        let m0 = building.get_material(m0_index).unwrap();
        let m0_substance_index = m0.get_substance_index().unwrap();
        let m0_substance = building.get_substance(m0_substance_index).unwrap();

        let m1 = building.get_material(m1_index).unwrap();
        let m1_substance_index = m1.get_substance_index().unwrap();
        let m1_substance = building.get_substance(m1_substance_index).unwrap();

        /* TESTS */
        let dt = 156.0;
        let n_nodes: Vec<usize> = vec![0, 1];
        let all_nodes = calc_n_total_nodes(&n_nodes).unwrap();
        assert_eq!(all_nodes, 2);

        let cp = m0_substance.specific_heat_capacity().unwrap();
        let rho = m0_substance.density().unwrap();
        let dx = m0.thickness().unwrap();

        let mut k_prime = Matrix::new(0.0, all_nodes, all_nodes);
        let mut full_rs_front = 0.0;
        let mut full_rs_back = 0.0;
        let mut c_i = 0.0;
        let mut c_o = 0.0;
        build_thermal_network(
            &building,
            &c0,
            dt,
            &n_nodes,
            0.,
            0.,
            &mut k_prime,
            &mut full_rs_front,
            &mut full_rs_back,
            &mut c_i,
            &mut c_o,
        )
        .unwrap();

        let rsi = m1.thickness().unwrap() / m1_substance.thermal_conductivity().unwrap();

        // check coefficients
        assert_eq!(full_rs_front, rsi); // 2.0*dt/(rho*cp*dx));
        assert_eq!(full_rs_back, 0.0); // 2.0*dt/(rho*cp*dx));
        assert_eq!(c_i, rho * cp * dx / (2.0 * dt));
        assert_eq!(c_o, rho * cp * dx / (2.0 * dt));

        let half_mass = m0_substance.density().unwrap()
            * m0_substance.specific_heat_capacity().unwrap()
            * m0.thickness().unwrap()
            / (2.0 * dt);
        let u = m0_substance.thermal_conductivity().unwrap() / m0.thickness().unwrap();

        // Check first node

        let found = k_prime.get(0, 0).unwrap();
        let exp = -(u + 1.0 / rsi) / half_mass;
        assert!((exp - found).abs() < 1E-10);

        // Check last node
        let found = k_prime.get(1, 1).unwrap();
        let exp = -u / half_mass;
        assert!((exp - found).abs() < 1E-10);
    }

    #[test]
    fn test_wall_7() {
        let mut building = Building::new("A building".to_string());

        /* SUBSTANCES */

        let poly_index = add_polyurethane(&mut building);
        let brickwork_index = add_brickwork(&mut building);

        /* MATERIALS */

        let m0_thickness = 200.0 / 1000. as f64;
        let m0_index = add_layer(&mut building, brickwork_index, m0_thickness);

        let m1_thickness = 200.0 / 1000. as f64;
        let m1_index = add_layer(&mut building, poly_index, m1_thickness);

        /* WALL */

        let c0_index = building.add_construction("Wall".to_string());
        building
            .add_material_to_construction(c0_index, m0_index)
            .unwrap();
        building
            .add_material_to_construction(c0_index, m1_index)
            .unwrap();
        building
            .add_material_to_construction(c0_index, m1_index)
            .unwrap();
        building
            .add_material_to_construction(c0_index, m0_index)
            .unwrap();

        let c0 = building.get_construction(c0_index).unwrap();

        let m0 = building.get_material(m0_index).unwrap();
        let m0_substance_index = m0.get_substance_index().unwrap();
        let m0_substance = building.get_substance(m0_substance_index).unwrap();

        let m1 = building.get_material(m1_index).unwrap();
        let m1_substance_index = m1.get_substance_index().unwrap();
        let m1_substance = building.get_substance(m1_substance_index).unwrap();

        let cp = m0_substance.specific_heat_capacity().unwrap();
        let rho = m0_substance.density().unwrap();
        let dx = m0.thickness().unwrap();

        let dt = 156.0;
        let n_nodes: Vec<usize> = vec![1, 0, 0, 1];
        let all_nodes = calc_n_total_nodes(&n_nodes).unwrap();
        assert_eq!(all_nodes, 4);

        let mut k_prime = Matrix::new(0.0, all_nodes, all_nodes);
        let mut full_rs_front = 0.0;
        let mut full_rs_back = 0.0;
        let mut c_i = 0.0;
        let mut c_o = 0.0;
        build_thermal_network(
            &building,
            &c0,
            dt,
            &n_nodes,
            0.,
            0.,
            &mut k_prime,
            &mut full_rs_front,
            &mut full_rs_back,
            &mut c_i,
            &mut c_o,
        )
        .unwrap();

        // check coefficients
        assert_eq!(full_rs_front, 0.0); // 2.0*dt/(rho*cp*dx));
        assert_eq!(full_rs_back, 0.0); // 2.0*dt/(rho*cp*dx));
        assert_eq!(c_i, rho * cp * dx / (2.0 * dt));
        assert_eq!(c_o, rho * cp * dx / (2.0 * dt));

        let half_mass = m0_substance.density().unwrap()
            * m0_substance.specific_heat_capacity().unwrap()
            * m0.thickness().unwrap()
            / (2.0 * dt);
        let u = m0_substance.thermal_conductivity().unwrap() / m0.thickness().unwrap();
        let insulation_r =
            m1.thickness().unwrap() * 2.0 / m1_substance.thermal_conductivity().unwrap();

        // Check first node
        let found = k_prime.get(0, 0).unwrap();
        let exp = -u / half_mass;
        assert!((exp - found).abs() < 1E-10);

        // Check second node
        let found = k_prime.get(1, 1).unwrap();
        let exp = -(u + 1.0 / insulation_r) / half_mass;
        assert!((exp - found).abs() < 1E-10);

        // Check third node
        let found = k_prime.get(2, 2).unwrap();
        let exp = -(u + 1.0 / insulation_r) / half_mass;
        assert!((exp - found).abs() < 1E-6);

        // Check last node
        let found = k_prime.get(3, 3).unwrap();
        let exp = -u / half_mass;
        assert!((exp - found).abs() < 1E-10);
    }

    #[test]
    fn test_calc_heat_flow_with_mass() {
        let mut building = Building::new("A building".to_string());

        /* SUBSTANCES */

        let poly_index = add_polyurethane(&mut building);

        /* MATERIALS */

        let m0_thickness = 200.0 / 1000. as f64;
        let m0_index = add_layer(&mut building, poly_index, m0_thickness);

        /* CONSTRUCTION */
        let c0_index = building.add_construction("Wall".to_string());
        building
            .add_material_to_construction(c0_index, m0_index)
            .unwrap();

        /* GEOMETRY */
        let mut the_loop = Loop3D::new();
        let l = 1. as f64;
        the_loop.push(Point3D::new(-l, -l, 0.)).unwrap();
        the_loop.push(Point3D::new(l, -l, 0.)).unwrap();
        the_loop.push(Point3D::new(l, l, 0.)).unwrap();
        the_loop.push(Point3D::new(-l, l, 0.)).unwrap();
        the_loop.close().unwrap();
        let p = Polygon3D::new(the_loop).unwrap();

        /* SURFACE */
        let s_index = building.add_surface("Surface".to_string());
        building.set_surface_polygon(s_index, p).unwrap();
        building
            .set_surface_construction(s_index, c0_index)
            .unwrap();

        /* TEST */
        let surface = building.get_surface(s_index).unwrap();
        let m0 = building.get_material(m0_index).unwrap();
        let c0 = building.get_construction(c0_index).unwrap();

        let main_dt = 300.0;
        let max_dx = m0.thickness().unwrap() / 4.0;
        let min_dt = 1.0;
        let (n_subdivisions, nodes) =
            discretize_construction(&building, &c0, main_dt, max_dx, min_dt);

        let dt = main_dt / n_subdivisions as f64;
        let mut state: SimulationState = SimulationState::new();
        let ts =
            ThermalSurface::from_surface(&building, &mut state, &surface, dt, &nodes, 0).unwrap();
        assert!(ts.massive);

        // MAP THE STATE
        building.map_simulation_state(&mut state).unwrap();

        // TEST
        let temperatures = ts.get_node_temperatures(&building, &state);

        assert_eq!(22.0, temperatures.get(0, 0).unwrap());
        assert_eq!(22.0, temperatures.get(ts.n_nodes - 1, 0).unwrap());

        let t_front = 10.0;
        let q_in_ref = (22.0 - t_front) / ts.rs_front;

        let t_back = 10.0;
        let q_out_ref = (22.0 - t_back) / ts.rs_back;
        let (q_in, q_out) = ts.calc_heat_flow(&building, &state, t_front, t_back);

        assert_eq!(q_in, q_in_ref);
        assert_eq!(q_out, q_out_ref);
    }

    #[test]
    fn test_calc_heat_flow_no_mass() {
        let mut building = Building::new("A building".to_string());

        /* SUBSTANCES */

        let poly_index = add_polyurethane(&mut building);

        /* MATERIALS */
        let m0_thickness = 200.0 / 1000. as f64;
        let m0_index = add_layer(&mut building, poly_index, m0_thickness);

        /* CONSTRUCTION */
        let c0_index = building.add_construction("Wall".to_string());
        building
            .add_material_to_construction(c0_index, m0_index)
            .unwrap();

        /* GEOMETRY */
        let mut the_loop = Loop3D::new();
        let l = 1. as f64;
        the_loop.push(Point3D::new(-l, -l, 0.)).unwrap();
        the_loop.push(Point3D::new(l, -l, 0.)).unwrap();
        the_loop.push(Point3D::new(l, l, 0.)).unwrap();
        the_loop.push(Point3D::new(-l, l, 0.)).unwrap();
        the_loop.close().unwrap();
        let p = Polygon3D::new(the_loop).unwrap();

        /* SURFACE */
        let s_index = building.add_surface("Surface".to_string());
        building.set_surface_polygon(s_index, p).unwrap();
        building
            .set_surface_construction(s_index, c0_index)
            .unwrap();

        /* TEST */
        let surface = building.get_surface(s_index).unwrap();
        let m0 = building.get_material(m0_index).unwrap();

        let main_dt = 300.0;
        let max_dx = m0.thickness().unwrap() / 15.;
        let min_dt = 80.;
        let (n_subdivisions, nodes) = {
            let c0 = building.get_construction(c0_index).unwrap();
            discretize_construction(&building, &c0, main_dt, max_dx, min_dt)
        };

        let dt = main_dt / n_subdivisions as f64;
        let mut state: SimulationState = SimulationState::new();
        let ts =
            ThermalSurface::from_surface(&building, &mut state, &surface, dt, &nodes, 0).unwrap();

        // MAP THE STATE
        building.map_simulation_state(&mut state).unwrap();

        // TEST

        assert!(!ts.massive);
        let temperatures = ts.get_node_temperatures(&building, &state);
        assert_eq!(22.0, temperatures.get(0, 0).unwrap());
        assert_eq!(22.0, temperatures.get(ts.n_nodes - 1, 0).unwrap());

        let t_front = 10.0;
        let q_in_ref = (22.0 - t_front) / ts.rs_front;

        let t_back = 10.0;
        let q_out_ref = (22.0 - t_back) / ts.rs_back;
        let (q_in, q_out) = ts.calc_heat_flow(&building, &state, t_front, t_back);

        assert_eq!(q_in, q_in_ref);
        assert_eq!(q_out, q_out_ref);
    }

    #[test]
    fn test_calc_heat_flow_mixed_mass() {
        // THIRD TEST -- WITH ONE NO-MASS LAYER IN THE EXTERIOR AND ONE IN THE INTERIOR

        let mut building = Building::new("Some building".to_string());

        /* SUBSTANCES */
        let poly_index = add_polyurethane(&mut building);
        let brickwork_index = add_brickwork(&mut building);

        /* MATERIALS */
        let m1_index = add_layer(&mut building, poly_index, 20. / 1000.);
        let m2_index = add_layer(&mut building, brickwork_index, 220. / 1000.);

        /* CONSTRUCTION */

        let c_index = building.add_construction("construction".to_string());
        building
            .add_material_to_construction(c_index, m1_index)
            .unwrap();
        building
            .add_material_to_construction(c_index, m2_index)
            .unwrap();
        building
            .add_material_to_construction(c_index, m1_index)
            .unwrap();

        /* GEOMETRY */

        let mut the_loop = Loop3D::new();
        let l = 1. as f64;
        the_loop.push(Point3D::new(-l, -l, 0.)).unwrap();
        the_loop.push(Point3D::new(l, -l, 0.)).unwrap();
        the_loop.push(Point3D::new(l, l, 0.)).unwrap();
        the_loop.push(Point3D::new(-l, l, 0.)).unwrap();
        the_loop.close().unwrap();
        let p = Polygon3D::new(the_loop).unwrap();

        /* SURFACE */

        let s_index = building.add_surface("Surface 1".to_string());
        building.set_surface_construction(s_index, c_index).unwrap();
        building.set_surface_polygon(s_index, p).unwrap();

        /* TEST */

        let surface = building.get_surface(s_index).unwrap();

        let main_dt = 300.0;
        let max_dx = 0.015;
        let min_dt = 65.;
        let (n_subdivisions, nodes) = {
            let c = building.get_construction(c_index).unwrap();
            discretize_construction(&building, &c, main_dt, max_dx, min_dt)
        };

        let dt = main_dt / n_subdivisions as f64;
        let mut state: SimulationState = SimulationState::new();
        let ts =
            ThermalSurface::from_surface(&building, &mut state, &surface, dt, &nodes, 0).unwrap();

        // MAP THE STATE
        building.map_simulation_state(&mut state).unwrap();

        // TEST

        assert!(ts.massive);
        let temperatures = ts.get_node_temperatures(&building, &state);
        assert_eq!(22.0, temperatures.get(0, 0).unwrap());
        assert_eq!(22.0, temperatures.get(ts.n_nodes - 1, 0).unwrap());

        let t_front = 10.0;
        let q_in_ref = (22.0 - t_front) / ts.full_rs_front;

        let t_back = 10.0;
        let q_out_ref = (22.0 - t_back) / ts.full_rs_back;
        let (q_in, q_out) = ts.calc_heat_flow(&building, &state, t_front, t_back);

        assert_eq!(q_in, q_in_ref);
        assert_eq!(q_out, q_out_ref);
    }

    #[test]
    fn test_march_massive() {
        let mut building = Building::new("Some building".to_string());

        /* SUBSTANCES */
        let brickwork_index = add_brickwork(&mut building);

        /* MATERIALS */
        let m1_index = add_layer(&mut building, brickwork_index, 200. / 1000.);

        /* CONSTRUCTION */

        let c_index = building.add_construction("construction".to_string());
        building
            .add_material_to_construction(c_index, m1_index)
            .unwrap();

        /* GEOMETRY */
        let mut the_loop = Loop3D::new();
        let l = 1. as f64;
        the_loop.push(Point3D::new(-l, -l, 0.)).unwrap();
        the_loop.push(Point3D::new(l, -l, 0.)).unwrap();
        the_loop.push(Point3D::new(l, l, 0.)).unwrap();
        the_loop.push(Point3D::new(-l, l, 0.)).unwrap();
        the_loop.close().unwrap();
        let p = Polygon3D::new(the_loop).unwrap();

        /* SURFACE */
        let s_index = building.add_surface("Surface 1".to_string());
        building.set_surface_polygon(s_index, p).unwrap();
        building.set_surface_construction(s_index, c_index).unwrap();

        let m1 = building.get_material(m1_index).unwrap();
        let surface = building.get_surface(s_index).unwrap();

        // FIRST TEST -- 10 degrees on each side
        let main_dt = 300.0;
        let max_dx = m1.thickness().unwrap() / 2.0;
        let min_dt = 1.0;
        let (n_subdivisions, nodes) = {
            let c = building.get_construction(c_index).unwrap();
            discretize_construction(&building, &c, main_dt, max_dx, min_dt)
        };
        let dt = main_dt / n_subdivisions as f64;

        let mut state: SimulationState = SimulationState::new();
        let ts =
            ThermalSurface::from_surface(&building, &mut state, surface, dt, &nodes, 0).unwrap();
        assert!(ts.massive);

        // MAP THE STATE
        building.map_simulation_state(&mut state).unwrap();

        // TEST

        // Try marching until q_in and q_out are zero.
        let mut q: f64 = 9000009.0;
        let c = building.get_construction(c_index).unwrap();
        let mut counter: usize = 0;
        while q.abs() > 1E-5 {
            let (q_in, q_out) = ts.march(&building, &mut state, 10.0, 10.0);
            // the same amount of heat needs to leave in each direction
            assert!((q_in - q_out).abs() < 1E-5);

            // q_front is positive
            assert!(q_in >= 0.);
            assert!(q_out >= 0.);

            // q_in needs to be getting smaller
            assert!(q_in < q);
            q = q_in;

            counter += 1;
            if counter > 99999 {
                panic!("Exceded number of iterations")
            }
        }

        // all nodes should be at 10.0 now.
        let temperatures = ts.get_node_temperatures(&building, &state);
        for i in 0..ts.n_nodes {
            let t = temperatures.get(i, 0).unwrap();
            assert!((t - 10.0).abs() < 1E-5);
        }

        // SECOND TEST -- 10 degrees in and 30 out.
        // We expect the final q to be (30-10)/R from
        // outside to inside. I.E. q_in = (30-10)/R,
        // q_out = -(30-10)/R

        // March until q converges
        let mut change: f64 = 99.0;
        let mut counter: usize = 0;
        let mut previous_q: f64 = -125.0;
        let mut final_qin: f64 = -12312.;
        let mut final_qout: f64 = 123123123.;
        while change.abs() > 1E-10 {
            let (q_in, q_out) = ts.march(&building, &mut state, 10.0, 30.0);
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
        assert!((final_qin + final_qout).abs() < 1E-6);

        let r = r_value(&building, c).unwrap();

        let exp_q = (30.0 - 10.0) / (r + ts.rs_front + ts.rs_back);
        assert!((exp_q - final_qin).abs() < 1E-4);
        assert!((exp_q + final_qout).abs() < 1E-4);
    }

    #[test]
    fn test_march_nomass() {
        let mut building = Building::new("A building".to_string());

        /* SUBSTANCE */
        let brickwork_index = add_brickwork(&mut building);

        /* MATERIAL */
        let m1_index = add_layer(&mut building, brickwork_index, 3. / 1000.);

        /* CONSTRUCTION */
        let c_index = building.add_construction("Construction".to_string());
        building
            .add_material_to_construction(c_index, m1_index)
            .unwrap();

        /* GEOMETRY */
        let mut the_loop = Loop3D::new();
        let l = 1. as f64;
        the_loop.push(Point3D::new(-l, -l, 0.)).unwrap();
        the_loop.push(Point3D::new(l, -l, 0.)).unwrap();
        the_loop.push(Point3D::new(l, l, 0.)).unwrap();
        the_loop.push(Point3D::new(-l, l, 0.)).unwrap();
        the_loop.close().unwrap();
        let p = Polygon3D::new(the_loop).unwrap();

        /* SURFACE */
        let s_index = building.add_surface("WALL".to_string());
        building.set_surface_construction(s_index, c_index).unwrap();
        building.set_surface_polygon(s_index, p).unwrap();

        let mut state: SimulationState = SimulationState::new();

        /* TEST */

        let m1 = building.get_material(m1_index).unwrap();
        let surface = building.get_surface(s_index).unwrap();

        // FIRST TEST -- 10 degrees on each side
        let main_dt = 300.0;
        let max_dx = m1.thickness().unwrap() / 2.0;
        let min_dt = 100.0;
        let (n_subdivisions, nodes) = {
            let c = building.get_construction(c_index).unwrap();
            discretize_construction(&building, c, main_dt, max_dx, min_dt)
        };

        let dt = main_dt / n_subdivisions as f64;

        let ts =
            ThermalSurface::from_surface(&building, &mut state, surface, dt, &nodes, 0).unwrap();
        assert!(!ts.massive);

        // MAP THE STATE
        building.map_simulation_state(&mut state).unwrap();

        // TEST

        let c = building.get_construction(c_index).unwrap();
        // Try marching until q_in and q_out are zero.

        let (q_in, q_out) = ts.march(&building, &mut state, 10.0, 10.0);

        // this should show instantaneous update. So,
        let temperatures = ts.get_node_temperatures(&building, &state);
        assert_eq!(temperatures.get(0, 0).unwrap(), 10.0);
        assert_eq!(temperatures.get(1, 0).unwrap(), 10.0);
        assert_eq!(q_in, 0.0);
        assert_eq!(q_out, 0.0);

        // SECOND TEST -- 10 degrees in and 30 out.
        // We expect the final q to be (30-10)/R from
        // outside to inside. I.E. q_in = (30-10)/R,
        // q_out = -(30-10)/R

        let (q_in, q_out) = ts.march(&building, &mut state, 10.0, 30.0);

        // Expecting
        assert!(q_in > 0.0);
        assert!(q_out < 0.0);
        assert!((q_in + q_out).abs() < 1E-6);

        let exp_q = (30.0 - 10.0) / (r_value(&building, c).unwrap() + ts.rs_front + ts.rs_back);
        assert!((exp_q - q_in).abs() < 1E-4);
        assert!((exp_q + q_out).abs() < 1E-4);
    }
}
