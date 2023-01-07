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

use crate::convection::ConvectionParams;
use crate::discretization::Discretization;
use crate::glazing::Glazing;
use crate::surface_trait::SurfaceTrait;
use crate::Float;
use geometry3d::Vector3D;
use matrix::Matrix;
use simple_model::{
    Boundary, Construction, Fenestration, SimpleModel, SimulationStateHeader, Substance, Surface,
    TerrainClass,
};
use simple_model::{SimulationState, SiteDetails};
use std::rc::Rc;

/// Calculates whether a surface is facing the wind direction
/// **wind_direction in Radians**
pub fn is_windward(wind_direction: Float, cos_tilt: Float, normal: Vector3D) -> bool {
    if cos_tilt.abs() < 0.98 {
        // tilted
        let wind_direction = Vector3D::new(wind_direction.sin(), wind_direction.cos(), 0.0);
        normal * wind_direction > 0.0
    } else {
        // if it is horizontal
        true
    }
}


/// The memory needed to simulate the marching forward
/// of a massive chunk
pub struct ChunkMemory {
    /// memory for a matrix
    pub aux: Matrix,
    /// memory for a matrix
    pub k: Matrix,
    /// memory for a matrix
    pub q: Matrix,
    /// memory for a matrix
    pub k1: Matrix,
    /// memory for a matrix
    pub k2: Matrix,
    /// memory for a matrix
    pub k3: Matrix,
    /// memory for a matrix
    pub k4: Matrix,
}

impl ChunkMemory {
    fn new(ini: usize, fin: usize)->Self{
        let n = fin - ini - 1;
        ChunkMemory {
            aux: Matrix::new(0.0, n+1, 1),
            k: Matrix::new(0.0, n+1, n+1),
            q: Matrix::new(0.0, n+1, 1),
            k1: Matrix::new(0.0, n+1, 1),
            k2: Matrix::new(0.0, n+1, 1),
            k3: Matrix::new(0.0, n+1, 1),
            k4: Matrix::new(0.0, n+1, 1),
        }
    }
}

/// The memory needed to simulate the marching of 
/// a surface
pub struct SurfaceMemory {
    /// Memory for each massive chunk
    pub massive_chunks: Vec<ChunkMemory>,
    /// Memory for each no-mass chunk
    pub nomass_chunks: Vec<ChunkMemory>,
    /// The temperatures
    pub temperatures: Matrix,

    /// The solar absorption on each node
    pub q: Matrix
}

/// Calculates a surface's wind speed modifier; that is to say, the value by which
/// the weather file wind speed needs to be multiplied in order to estimate the wind
/// speed next to the window
///
/// This is a rip off from EnergyPlus' Engineering Reference, where they explain that
/// the corrected wind speed ($`V_z`$ in $`m/s`$) at a certain altitude $`z`$ in $`m`$
/// (e.g., the height of the window's centroid) can be estimated through an equation that
/// relates the measurements at the meteorological station and those in the site.
///
/// Specifically, this equation depends on the altitude
/// at which the wind speed was measured at the meteorological station ($`z_{met}`$,
/// assumed to be $`10m`$), the so-called "wind speed profile boundary layer" at the
/// weather station ($`\delta_{met}`$, assumed to be $`240m`$) and the "wind speed profile
/// exponent" at the meteorological station $`\alpha_{met}`$. Also, it depends on the
/// "wind speed profile boundary layer" at the site ($`\delta`$) and the "wind speed profile
/// exponent" $`\alpha`$.
///
/// ```math
/// V_z = \left(\frac{\delta_{met}}{z_{met}}\right)^{\alpha_{met}}\left(\frac{z}{\delta} \right)^{\alpha}
/// ```
/// The values for $`\alpha`$ and $`\delta`$ depend on the kind of terrain.
///
/// | Terrain Class | $`\alpha`$ | $`\delta`$ |
/// |---------------|------------|------------|
/// | Country       | 0.14       | 270        |
/// | Suburbs       | 0.22       | 370        |
/// | City          | 0.33       | 460        |
/// | Ocean         | 0.10       | 210        |
/// | Urban         | 0.22       | 370        |
///
/// > Note: if height is Zero, then we assume the wind speed to be Zero
pub fn wind_speed_modifier(height: Float, site_details: &Option<SiteDetails>) -> Float {
    // Surface touching the ground... no wind
    if height < 1e-5 {
        return 0.0;
    }
    // this bit does not change... it is assumed constant for all meterological stations

    let mut alpha = 0.0;
    let mut delta = 0.0;

    if let Some(details) = site_details {
        // raise all surfaces, if needed
        // if let Ok(z_terrain) = details.altitude(){
        //     height += z_terrain;
        // }
        if let Ok(terrain) = details.terrain() {
            (alpha, delta) = match terrain {
                TerrainClass::Country => (0.14, 270.),
                TerrainClass::Suburbs => (0.22, 370.),
                TerrainClass::City => (0.33, 460.),
                TerrainClass::Ocean => (0.10, 210.),
                TerrainClass::Urban => (0.22, 370.),
            }
        }
    } else {
        // default to Urban
        alpha = 0.22;
        delta = 370.;
    }

    (270. / 10. as Float).powf(0.14) * (height / delta).powf(alpha)
}

fn rearrange_k(dt: Float, c: &Matrix, memory: &mut ChunkMemory) -> Result<(), String> {
    let (crows, ..) = c.size();
    // Rearrenge into dT = (dt/C) * K + (dt/C)*q
    for nrow in 0..crows {
        let v = dt / c.get(nrow, nrow)?;
        // transform k into k_prime (i.e., k * dt/C)
        let ini = if nrow == 0 { 0 } else { nrow - 1 };
        let fin = if nrow == crows - 1 {
            crows - 1
        } else {
            nrow + 1
        };
        // for ncol in 0..kcols{
        for ncol in ini..=fin {
            memory.k.scale_element(nrow, ncol, v)?;
        }
        memory.q.scale_element(nrow, 0, v)?;
    }
    Ok(())
}

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
fn rk4(
    // dt: Float,
    c: &Matrix,
    /*mut k: Matrix, mut q: Matrix,*/
    memory: &mut ChunkMemory,
    t: &mut Matrix,
) -> Result<(), String> {
    let (krows, kcols) = memory.k.size();
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
    let (qrows, qcols) = memory.q.size();
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

    // // Rearrenge into dT = (dt/C) * K + (dt/C)*q
    // for nrow in 0..krows {
    //     let v = dt / c.get(nrow, nrow)?;
    //     // transfrom q into q_prime (i.e., q * dt/C)
    //     memory.q.scale_element(nrow, 0, v)?;
    // }

    // I am not sure why I need to clean... I thought this was not necessary.
    memory.k1 *= 0.0;
    memory.k2 *= 0.0;
    memory.k3 *= 0.0;
    memory.k4 *= 0.0;
    memory.aux *= 0.0;

    // get k1
    memory.k.prod_tri_diag_into(t, &mut memory.k1)?;
    memory.k1 += &memory.q;
    
    // returning "temperatures + k1" is Euler... continuing is
    // Runge–Kutta 4th order
    // *t += &memory.k1;    
    // return Ok(());

    memory.k1.scale_into(0.5, &mut memory.aux)?;
    memory.aux += t;    

    // k2
    memory.k.prod_tri_diag_into(&memory.aux, &mut memory.k2)?;
    memory.k2 += &memory.q;    

    // k3
    memory.k2.scale_into(0.5, &mut memory.aux)?;
    memory.aux += t;
    memory.k.prod_tri_diag_into(&memory.aux, &mut memory.k3)?;
    memory.k3 += &memory.q;    

    // k4
    memory.aux.copy_from(&memory.k3);
    memory.aux += t;
    memory.k.prod_tri_diag_into(&memory.aux, &mut memory.k4)?;
    memory.k4 += &memory.q;    

    // Scale them and add them all up
    memory.k1 /= 6.;
    memory.k2 /= 3.;
    memory.k3 /= 3.;
    memory.k4 /= 6.;

    // Let's add it all and return
    *t += &memory.k1;
    *t += &memory.k2;
    *t += &memory.k3;
    *t += &memory.k4;

    Ok(())
}

/// This is a Surface from the point of view of our thermal solver.
/// Since this module only calculate heat transfer (and not short-wave solar
/// radiation, e.g., light), both simple_model::Fenestration and simple_model::Surface
/// are treated in the same way.
pub struct ThermalSurfaceData<T: SurfaceTrait> {
    /// A reference to the element in the [`SimpleModel`] which this struct represents
    pub parent: Rc<T>,

    /// The [`Discretization`] that represents this `ThermalSurfaceData`
    pub discretization: Discretization,

    /// The back boundary
    pub front_boundary: Option<Boundary>,

    /// The index of the space in front of this, if any
    pub front_space_index: Option<usize>,

    /// The front boundary
    pub back_boundary: Option<Boundary>,

    /// The index of the space at the back of this, if any
    pub back_space_index: Option<usize>,

    /// The thermal absorbtance on the front side (from 0 to 1)
    pub front_emissivity: Float,

    /// The thermal absorbtance on the back side (from 0 to 1)
    pub back_emissivity: Float,

    /// The area of the Surface
    pub area: Float,

    /// The perimeter of the surface
    pub perimeter: Float,

    /// The normal of the surface
    pub normal: Vector3D,

    /// The wind velocity changes with altitude. This field
    /// stores the factor with which the wind velocity of the weather
    /// file needs to be multiplied in order to estimate the wind speed
    /// at the exterior of the surface.
    pub wind_speed_modifier: Float,

    /// The cosine of the tilt angle (normal * Vector3D(0., 0., 1.))
    pub cos_tilt: Float,

    /// The chunks of nodes that have mass
    pub massive_chunks: Vec<(usize, usize)>,

    /// The chunks of nodes that have nomass
    pub nomass_chunks: Vec<(usize, usize)>,

    /// The absorbtances of each node in the system, proportional
    /// to the front incident radiation (i.e., they do not add up to 1.0)
    pub front_alphas: Matrix,

    /// The absorbtances of each node in the system, proportional
    /// to the back incident radiation (i.e., they do not add up to 1.0)
    pub back_alphas: Matrix,

    /// [**Only available during testing**] this allows setting a fixed convection
    /// coefficient
    #[cfg(debug_assertions)]
    pub front_hs: Option<Float>,

    /// [**Only available during testing**] this allows setting a fixed convection
    /// coefficient
    #[cfg(debug_assertions)]
    pub back_hs: Option<Float>,
}

impl<T: SurfaceTrait> ThermalSurfaceData<T> {
    /// Allocates memory for the simulation
    pub fn allocate_memory(&self) -> SurfaceMemory {
        let massive_chunks = self.massive_chunks.iter().map(|(ini, fin)|{
            ChunkMemory::new(*ini,*fin)
        }).collect();

        let nomass_chunks = self.nomass_chunks.iter().map(|(ini, fin)|{
            ChunkMemory::new(*ini,*fin)
        }).collect();

        let ini = self.parent.first_node_temperature_index();
        let fin = self.parent.last_node_temperature_index() + 1;
        let n_nodes = fin - ini;
        let q = Matrix::new(0.0, n_nodes, 1);
        let temperatures = Matrix::new(0.0, n_nodes, 1);

        SurfaceMemory{
            massive_chunks,
            nomass_chunks,
            temperatures,
            q,
        }
    }

    /// Creates a new [`ThermalSurfaceData`]
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        state: &mut SimulationStateHeader,
        model: &SimpleModel,
        site_details: &Option<SiteDetails>,
        ref_surface_index: usize,
        parent: &Rc<T>,
        area: Float,
        perimeter: Float,
        height: Float,
        normal: Vector3D,
        construction: &Rc<Construction>,
        discretization: Discretization,
    ) -> Result<ThermalSurfaceData<T>, String> {
        // Set Front and Back state
        parent.add_front_convection_state(state, ref_surface_index)?;
        parent.add_back_convection_state(state, ref_surface_index)?;

        parent.add_front_convective_heatflow_state(state, ref_surface_index)?;
        parent.add_back_convective_heatflow_state(state, ref_surface_index)?;

        parent.add_front_solar_irradiance_state(state, ref_surface_index)?;
        parent.add_back_solar_irradiance_state(state, ref_surface_index)?;

        parent.add_front_ir_irradiance_state(state, ref_surface_index)?;
        parent.add_back_ir_irradiance_state(state, ref_surface_index)?;

        // Add node data.
        let n_nodes = discretization.segments.len();
        parent.add_node_temperature_states(state, ref_surface_index, n_nodes)?;

        let front_substance = model.get_material_substance(&construction.materials[0])?;

        let back_substance =
            model.get_material_substance(construction.materials.last().unwrap())?;

        const DEFAULT_EM: Float = 0.84;
        let front_emissivity = match &front_substance {
            Substance::Normal(s) => {
                s.front_thermal_absorbtance_or(crate::model::MODULE_NAME, DEFAULT_EM)
            }
            _ => panic!("Front Emissivity not available for this particular kind of Substance"),
        };
        let back_emissivity = match &&back_substance {
            Substance::Normal(s) => {
                s.back_thermal_absorbtance_or(crate::model::MODULE_NAME, DEFAULT_EM)
            }
            _ => panic!("Front Emissivity not available for this particular kind of Substance"),
        };

        let (massive_chunks, nomass_chunks) = discretization.get_chunks();

        // Calculate solar absoption
        let front_glazing = Glazing::get_front_glazing_system(construction, model)?;
        let back_glazing = Glazing::get_back_glazing_system(construction, model)?;
        // These two are the absorbtion of each glazing layer. We need the absorption of each node
        let front_alphas_prev = Glazing::alphas(&front_glazing);
        if front_alphas_prev.len() != 1 && front_alphas_prev.len() != construction.materials.len() {
            panic!("Construction '{}' seems to have a mixture of transparent and opaque layers. This is not currently supported.", construction.name());
        }
        let n_nodes = discretization.segments.len();
        let n_layers = construction.materials.len();

        let mut front_alphas = Matrix::new(0.0, n_nodes, 1);
        let mut global_i = 0;
        for (alpha_i, alpha) in front_alphas_prev.iter().enumerate() {
            let layer_index = 2 * alpha_i; // We need to skip cavities
            let n = if discretization.n_elements[layer_index] == 0 {
                1
            } else {
                discretization.n_elements[layer_index]
            };
            let substance = model.get_material_substance(&construction.materials[layer_index])?;
            if let Substance::Normal(sub) = &substance {
                let tau = *sub.solar_transmittance().unwrap_or(&0.0);
                if tau > 0.0 {
                    // Distribute across all the nodes
                    for local_i in 0..=n {
                        front_alphas
                            .add_to_element(global_i + local_i, 0, *alpha / (n + 1) as Float)
                            .unwrap();
                    }
                } else {
                    // Add only to the first node
                    front_alphas.add_to_element(global_i, 0, *alpha).unwrap();
                }
            } else {
                unreachable!()
            }
            global_i += n + 1;
        }

        let back_alphas_prev = Glazing::alphas(&back_glazing);
        if back_alphas_prev.len() != 1 && back_alphas_prev.len() != construction.materials.len() {
            panic!("Construction '{}' seems to have a mixture of transparent and opaque layers. This is not currently supported.", construction.name());
        }
        let mut back_alphas = Matrix::new(0.0, n_nodes, 1);
        let mut global_i = n_nodes;
        for (alpha_i, alpha) in back_alphas_prev.iter().enumerate() {
            let layer_index = n_layers - 2 * alpha_i - 1; // We need to skip cavities
            let n = if discretization.n_elements[layer_index] == 0 {
                1
            } else {
                discretization.n_elements[layer_index]
            };
            let substance = model.get_material_substance(&construction.materials[layer_index])?;
            if let Substance::Normal(sub) = substance {
                let tau = *sub.solar_transmittance().unwrap_or(&0.0);
                if tau > 0.0 {
                    // Distribute across all the nodes
                    for local_i in 0..=n {
                        let row = global_i - local_i - 1;
                        back_alphas
                            .add_to_element(row, 0, *alpha / (n + 1) as Float)
                            .unwrap();
                    }
                } else {
                    // Add only to the last node
                    back_alphas.add_to_element(global_i - 1, 0, *alpha).unwrap();
                }
            } else {
                unreachable!()
            }
            global_i -= n + 1;
        }

        let cos_tilt = normal * Vector3D::new(0., 0., 1.);
        let wind_speed_modifier = wind_speed_modifier(height, site_details);

        // Build resulting
        Ok(ThermalSurfaceData {
            parent: parent.clone(),
            area,
            perimeter,
            normal,
            cos_tilt,
            discretization,
            front_boundary: None,
            back_boundary: None,
            front_space_index: None,
            back_space_index: None,
            front_emissivity,
            back_emissivity,
            wind_speed_modifier,
            front_alphas,
            back_alphas,
            massive_chunks,
            nomass_chunks,
            #[cfg(debug_assertions)]
            front_hs: None,
            #[cfg(debug_assertions)]
            back_hs: None,
        })
    }

    /// Sets the front boundary
    pub fn set_front_boundary(&mut self, b: &Boundary, model: &SimpleModel) {
        self.front_boundary = Some(b.clone());
        if let Boundary::Space { space } = b {
            for (i, s) in model.spaces.iter().enumerate() {
                if s.name() == space {
                    self.front_space_index = Some(i);
                    break;
                }
            }
        }
    }

    /// Sets the back boundary
    pub fn set_back_boundary(&mut self, b: &Boundary, model: &SimpleModel) {
        self.back_boundary = Some(b.clone());
        if let Boundary::Space { space } = b {
            for (i, s) in model.spaces.iter().enumerate() {
                if s.name() == space {
                    self.back_space_index = Some(i);
                    break;
                }
            }
        }
    }

    fn calc_border_conditions(
        &self,
        state: &SimulationState,
        t_front: Float,
        t_back: Float,
        wind_direction: Float,
        wind_speed: Float,
    ) -> (ConvectionParams, ConvectionParams, Float, Float) {
        // Calculate and set Front and Back IR Irradiance
        let ir_front = self.parent.front_infrared_irradiance(state);
        let ir_back = self.parent.back_infrared_irradiance(state);

        let windward = is_windward(wind_direction, self.cos_tilt, self.normal);

        // TODO: There is something to do here if we are talking about windows
        let (front_env, front_hs) = if let Some(b) = &self.front_boundary {
            match &b {
                Boundary::Space { .. } => {
                    let front_env = ConvectionParams {
                        air_temperature: t_front,
                        air_speed: 0.0,
                        rad_temperature: t_front,
                        surface_temperature: self.parent.front_temperature(state),
                        roughness_index: 1,
                        cos_surface_tilt: self.cos_tilt,
                    };

                    (
                        front_env,
                        front_env.get_tarp_natural_convection_coefficient(),
                    )
                }
                Boundary::AmbientTemperature { temperature } => {
                    let front_env = ConvectionParams {
                        air_temperature: *temperature,
                        air_speed: 0.0,
                        rad_temperature: t_front,
                        surface_temperature: self.parent.front_temperature(state),
                        roughness_index: 1,
                        cos_surface_tilt: self.cos_tilt,
                    };

                    (
                        front_env,
                        front_env.get_tarp_natural_convection_coefficient(),
                    )
                }
                Boundary::Ground => unreachable!(),
            }
        } else {
            let mut front_env = ConvectionParams {
                air_temperature: t_front,
                air_speed: wind_speed * self.wind_speed_modifier,
                rad_temperature: (ir_front / crate::SIGMA).powf(0.25) - 273.15,
                surface_temperature: self.parent.front_temperature(state),
                roughness_index: 1,
                cos_surface_tilt: self.cos_tilt,
            };
            front_env.cos_surface_tilt = -self.cos_tilt;
            (
                front_env,
                front_env.get_tarp_convection_coefficient(self.area, self.perimeter, windward),
            )
        };

        let (back_env, back_hs) = if let Some(b) = &self.back_boundary {
            match &b {
                Boundary::Space { .. } => {
                    let back_env = ConvectionParams {
                        air_temperature: t_back,
                        air_speed: 0.0,
                        rad_temperature: t_back, //self.parent.back_temperature(state),//(ir_back/crate::SIGMA).powf(0.25) - 273.15,
                        surface_temperature: self.parent.back_temperature(state),
                        roughness_index: 1,
                        cos_surface_tilt: self.cos_tilt,
                    };
                    (back_env, back_env.get_tarp_natural_convection_coefficient())
                }
                Boundary::AmbientTemperature { temperature } => {
                    let front_env = ConvectionParams {
                        air_temperature: *temperature,
                        air_speed: 0.0,
                        rad_temperature: t_front,
                        surface_temperature: self.parent.front_temperature(state),
                        roughness_index: 1,
                        cos_surface_tilt: self.cos_tilt,
                    };

                    (
                        front_env,
                        front_env.get_tarp_natural_convection_coefficient(),
                    )
                }
                Boundary::Ground => unreachable!(),
            }
        } else {
            let back_env = ConvectionParams {
                air_temperature: t_back,
                air_speed: wind_speed * self.wind_speed_modifier,
                rad_temperature: (ir_back / crate::SIGMA).powf(0.25) - 273.15,
                surface_temperature: self.parent.back_temperature(state),
                roughness_index: 1,
                cos_surface_tilt: self.cos_tilt,
            };
            (
                back_env,
                back_env.get_tarp_convection_coefficient(self.area, self.perimeter, windward),
            )
        };

        assert!(
            !front_hs.is_nan() && !back_hs.is_nan(),
            "Found NaN convection coefficients: Front={front_hs} | back={back_hs}"
        );
        #[cfg(debug_assertions)]
        return (
            front_env,
            back_env,
            self.front_hs.unwrap_or(front_hs),
            self.back_hs.unwrap_or(back_hs),
        );
        #[cfg(not(debug_assertions))]
        (front_env, back_env, front_hs, back_hs)
    }

    fn march_mass(&self,
        global_temperatures: &mut Matrix, 
        solar_radiation: &Matrix, 
        dt: Float,
        t_front: Float,
        t_back: Float,
        front_rad_hs: Float, 
        back_rad_hs: Float, 
        wind_direction: Float,
        wind_speed: Float,
        ini: usize, 
        fin: usize, 
        memory: &mut ChunkMemory,
        state: &SimulationState,
    )->Result<(),String>{

        

        
        let (front_env, back_env, front_hs, back_hs) =
            self.calc_border_conditions(state, t_front, t_back, wind_direction, wind_speed);

        self.discretization.get_k_q(
            ini,
            fin,
            &global_temperatures,
            &front_env,
            front_hs,
            front_rad_hs,
            &back_env,
            back_hs,
            back_rad_hs,
            memory,
        )?;

        let c = self
            .discretization
            .segments
            .iter()
            .skip(ini)
            .take(fin - ini)
            .map(|(mass, _)| *mass)
            .collect();
        let c = Matrix::diag(c);
        
        
        // ... here we add solar gains
        for (local_i, global_i) in (ini..fin).into_iter().enumerate() {
            let v = solar_radiation.get(global_i, 0).unwrap();
            memory.q.add_to_element(local_i, 0, v).unwrap();
        }
        
        rearrange_k(dt, &c, memory)?;

        // Use RT4 for updating temperatures of massive nodes.
        let mut local_temps = Matrix::new(0.0, fin - ini, 1);
        for (local_i, global_i) in (ini..fin).into_iter().enumerate() {
            let v = global_temperatures.get(global_i, 0).unwrap();
            local_temps.set(local_i, 0, v).unwrap();
        }

        
        rk4( &c, memory, &mut local_temps)?;

        for (local_i, global_i) in (ini..fin).into_iter().enumerate() {
            let v = local_temps.get(local_i, 0).unwrap();
            global_temperatures.set(global_i, 0, v).unwrap();
        }
        Ok(())
    }

    fn march_nomass(&self, 
        global_temperatures: &mut Matrix, 
        solar_radiation: &Matrix, 
        t_front: Float,
        t_back: Float,
        front_rad_hs: Float, 
        back_rad_hs: Float, 
        wind_direction: Float,
        wind_speed: Float,
        ini: usize, 
        fin: usize, 
        memory: &mut ChunkMemory,
        state: &SimulationState,
    )->Result<(),String>{
        let mut old_err = 99999.;
        let mut count = 0;
        
       

        loop {

            // Update convection coefficients
            let (front_env, back_env, front_hs, back_hs) =
                self.calc_border_conditions(state, t_front, t_back, wind_direction, wind_speed);

            
            // Calculate q based on heat transfer (convection, IR radiation)
            self.discretization.get_k_q(
                ini,
                fin,
                &global_temperatures,
                &front_env,
                front_hs,
                front_rad_hs,
                &back_env,
                back_hs,
                back_rad_hs,
                memory,
            )?;
            
            // add solar gains
            for (local_i, i) in (ini..fin).into_iter().enumerate() {
                let v = solar_radiation.get(i, 0)?;
                memory.q.add_to_element(local_i, 0, v)?;
            }
            memory.q *= -1.;

            let temps = memory.k.clone().mut_n_diag_gaussian(memory.q.clone(), 3)?; // and just like that, q is the new temperatures

            let mut err = 0.0;
            for (local_i, i) in (ini..fin).into_iter().enumerate() {
                let local_temp = temps.get(local_i, 0).unwrap();
                let global_temp = global_temperatures.get(i, 0)?;
                err += (local_temp - global_temp).abs();
            }
            if err > old_err {
                #[cfg(debug_assertions)]
                if count > 100 {
                    eprintln!("Breaking after {} iterations... because BAD!", count);
                }
                break;
            }

            assert!(
                !err.is_nan(),
                // "Error is NaN... \nfront_env = {:?}| back_env = {:?} \nfront_hc = {} | back_hs = {}. \nError = {}\ntemps={}\nq={}\nsolar_front={}, solar_back={}\nfront_alphas={}\nback_alphas={}\n",
                // front_env,
                // back_env,
                // front_hs,
                // back_hs,
                // err / ((fin - ini) as Float),
                // temps,                    
                // q,
                // solar_front,
                // solar_back,
                // self.front_alphas,
                // self.back_alphas,
            );

            // if count > 10000 {
            //     eprintln!("Err is {}", err / ((fin - ini) as Float))
            // }
            assert!(
                count < 99199000,
                "Excessive number of iterations... \n====\t\tfront_env = {:?}\n\tback_env = {:?}\n\tfront_hc = {}\n\tback_hs = {}.\n\tError = {}\n====\n",
                front_env,
                back_env,
                front_hs,
                back_hs,
                err / ((fin - ini) as Float),
            );
            for (local_i, i) in (ini..fin).into_iter().enumerate() {
                let local_temp = temps.get(local_i, 0).unwrap();
                // temperatures.set(i, 0, local_temp).unwrap();
                global_temperatures.add_to_element(i, 0, local_temp)?;
                global_temperatures.scale_element(i, 0, 0.5)?;
            }

            let max_allowed_error = if count < 100 { 0.01 } else /*if count < 1000*/ { 0.5 }; // else { 1. };

            if err / ((fin - ini) as Float) < max_allowed_error {
                #[cfg(debug_assertions)]
                if count > 100 {
                    dbg!("Breaking after {} iterations... because GOOD!", count);
                }
                break;
            }
            old_err = err;
            count += 1;
        }
        Ok(())
    }
    
    /// Marches one timestep. Returns front and back heat flow    
    pub fn march(
        &self,
        state: &mut SimulationState,
        t_front: Float,
        t_back: Float,
        wind_direction: Float,
        wind_speed: Float,
        dt: Float,
        memory: &mut SurfaceMemory,
    ) -> Result<(Float, Float), String> {
        let tempsssss = self.parent.get_node_temperatures(state);
        memory.temperatures.copy_from(&tempsssss);

        let (rows, ..) = memory.temperatures.size();

        // Calculate and set Front and Back Solar Irradiance
        let mut solar_front = self.parent.front_solar_irradiance(state);
        if solar_front.is_nan() || solar_front < 0.0 {
            solar_front = 0.0;
        }
        let mut solar_back = self.parent.back_solar_irradiance(state);
        if solar_back.is_nan() || solar_front < 0.0 {
            solar_back = 0.0;
        }

        /////////////////////
        // 1st: Calculate the solar absorption in each node
        /////////////////////        
        // memory.q *= 0.0; // clean, just in case
        // self.front_alphas.scale_into(solar_front, &mut memory.q)?;
        let mut solar_radiation = &self.front_alphas * solar_front;
        solar_radiation += &(&self.back_alphas * solar_back);
        // memory.q += &(&self.back_alphas * solar_back);

        /////////////////////
        // 2nd: Calculate the temperature in all no-mass nodes.
        // Also, the heat flow into
        /////////////////////

        let (front_env, back_env, _front_hs, _back_hs) =
            self.calc_border_conditions(state, t_front, t_back, wind_direction, wind_speed);
        let front_rad_hs = 4.
            * self.front_emissivity
            * crate::SIGMA
            * (273.15 + (front_env.rad_temperature + front_env.surface_temperature) / 2.).powi(3);
        let back_rad_hs = 4.
            * self.back_emissivity
            * crate::SIGMA
            * (273.15 + (back_env.rad_temperature + back_env.surface_temperature) / 2.).powi(3);
        
        for (chunk_i,(ini, fin)) in self.nomass_chunks.iter().enumerate() {
            self.march_nomass(
                &mut memory.temperatures, 
                &solar_radiation,// &memory.q,                
                t_front,
                t_back,
                front_rad_hs,
                back_rad_hs,
                wind_direction,
                wind_speed,
                *ini, *fin, 
                &mut memory.nomass_chunks[chunk_i],
                state,
            )?;
        }

        // Calculate final conditions.

        let (front_env, back_env, _front_hs, _back_hs) =
            self.calc_border_conditions(state, t_front, t_back, wind_direction, wind_speed);
        let front_rad_hs = 4.
            * self.front_emissivity
            * crate::SIGMA
            * (273.15 + (front_env.rad_temperature + front_env.surface_temperature) / 2.).powi(3);
        let back_rad_hs = 4.
            * self.back_emissivity
            * crate::SIGMA
            * (273.15 + (back_env.rad_temperature + back_env.surface_temperature) / 2.).powi(3);

        /////////////////////
        // 3rd: Calculate K and C matrices for the massive walls, and march
        /////////////////////
        
        for (chunk_i,(ini, fin)) in self.massive_chunks.iter().enumerate() {            
            self.march_mass(
                &mut memory.temperatures,
                &solar_radiation,// &memory.q,                
                dt,
                t_front,
                t_back,
                front_rad_hs,
                back_rad_hs,
                wind_direction,
                wind_speed,
                *ini, *fin, 
                &mut memory.massive_chunks[chunk_i],
                state,
            )?;
        }

        /////////////////////
        // 4th: Set temperatures, calc heat-flows and return
        /////////////////////
        self.parent.set_node_temperatures(state, &memory.temperatures);

        // Calc heat flow
        let ts_front = memory.temperatures.get(0, 0).unwrap();
        let ts_back = memory.temperatures.get(rows - 1, 0).unwrap();
        let (_front_env, _back_env, front_hs, back_hs) =
            self.calc_border_conditions(state, t_front, t_back, wind_direction, wind_speed);
        self.parent
            .set_front_convection_coefficient(state, front_hs)?;
        self.parent
            .set_back_convection_coefficient(state, back_hs)?;

        let flow_front = (ts_front - t_front) * front_hs;
        let flow_back = (ts_back - t_back) * back_hs;

        Ok((flow_front, flow_back))
    }
}

/// A [`ThermalSurfaceData`] whose parent is a [`Surface`]
pub type ThermalSurface = ThermalSurfaceData<Surface>;

/// A [`ThermalSurfaceData`] whose parent is a [`Fenestration`]
pub type ThermalFenestration = ThermalSurfaceData<Fenestration>;

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
            .set_front_thermal_absorbtance(0.)
            .set_back_thermal_absorbtance(0.)
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
            .set_front_thermal_absorbtance(0.)
            .set_back_thermal_absorbtance(0.)
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
        let mat = Material::new("123123".to_string(), substance.name().clone(), thickness);

        model.add_material(mat)
    }

    #[test]
    fn test_march_massive() {
        let mut model = SimpleModel::default();

        /* SUBSTANCES */
        let brickwork = add_brickwork(&mut model);

        /* MATERIALS */
        let m1 = add_material(&mut model, brickwork, 20. / 1000.);

        /* CONSTRUCTION */
        let mut c = Construction::new("construction".to_string());
        c.materials.push(m1.name().clone());
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
        let s = Surface::new("Surface 1".to_string(), p, c.name().clone());

        let surface = model.add_surface(s);

        // FIRST TEST -- 10 degrees on each side
        let main_dt = 300.0;
        let max_dx = m1.thickness / 2.0;
        let min_dt = 1.0;
        let d = Discretization::new(&c, &model, main_dt, max_dx, min_dt, 1., 0.).unwrap();
        let dt = main_dt / d.tstep_subdivision as Float;
        let normal = geometry3d::Vector3D::new(0., 0., 1.);
        let perimeter = 8. * l;
        let mut state_header = SimulationStateHeader::new();
        let mut ts = ThermalSurface::new(
            &mut state_header,
            &model,
            &None,
            0,
            &surface,
            surface.area(),
            perimeter,
            10.,
            normal,
            &c,
            d,
        )
        .unwrap();

        let mut memory = ts.allocate_memory();

        ts.front_hs = Some(10.);
        ts.back_hs = Some(10.);

        let mut state = state_header.take_values().unwrap();

        // TEST

        // Try marching until q_in and q_out are zero.
        let mut q: Float = 9999000009.0;
        let mut counter: usize = 0;
        let t_environment = 10.;
        let v = crate::SIGMA * (t_environment + 273.15 as Float).powi(4);
        while q.abs() > 0.00015 {            
            ts.parent.set_front_ir_irradiance(&mut state, v).unwrap();
            ts.parent.set_back_ir_irradiance(&mut state, v).unwrap();
            let (q_out, q_in) = ts
                .march(
                    &mut state,
                    t_environment,
                    t_environment,
                    0.0,
                    0.0,
                    dt,
                    &mut memory,
                )
                .unwrap();

            // the same amount of heat needs to leave in each direction
            // println!("q_in = {}, q_out = {} | diff = {}", q_in, q_out, (q_in - q_out).abs());
            assert!(
                (q_in - q_out).abs() < 0.5, //1E-5,
                "diff is {} (count is {counter})",
                (q_in - q_out).abs()
            );

            // q_front is positive
            assert!(q_in >= 0., "q_in = {} | c = {}", q_in, counter);
            assert!(q_out >= 0., "q_out = {} | c = {}", q_out, counter);

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
            assert!(
                (t - 10.0).abs() < 0.002,
                "Error found is {}",
                (t - 10.0).abs()
            );
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
            let (q_front, q_back) = ts
                .march(&mut state, 10.0, 30.0, 0.0, 0.0, dt, &mut memory)
                .unwrap();

            ts.parent
                .set_front_ir_irradiance(&mut state, crate::SIGMA * (10. + 273.15 as Float).powi(4))
                .unwrap();
            ts.parent
                .set_back_ir_irradiance(&mut state, crate::SIGMA * (30. + 273.15 as Float).powi(4))
                .unwrap();

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
    }

    #[test]
    fn test_march_nomass() {
        let mut model = SimpleModel::default();

        /* SUBSTANCE */
        let polyurethane = add_polyurethane(&mut model);

        /* MATERIAL */
        let m1 = add_material(&mut model, polyurethane, 3. / 1000.);

        /* CONSTRUCTION */
        let mut c = Construction::new("Construction");
        c.materials.push(m1.name().clone());
        c.materials.push(m1.name().clone());
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
        let s = Surface::new("WALL".to_string(), p, c.name().clone());
        let surface = model.add_surface(s);
        let mut state_header = SimulationStateHeader::new();

        /* TEST */

        let main_dt = 3.0;
        let max_dx = m1.thickness / 7.0;
        let min_dt = 10.0;
        let d = Discretization::new(&c, &model, main_dt, max_dx, min_dt, 1., 0.).unwrap();
        let dt = main_dt / d.tstep_subdivision as Float;

        let normal = geometry3d::Vector3D::new(0., 0., 1.);
        let perimeter = 8. * l;
        let mut ts = ThermalSurface::new(
            &mut state_header,
            &model,
            &None,
            0,
            &surface,
            surface.area(),
            perimeter,
            10.,
            normal,
            &c,
            d,
        )
        .unwrap();
        ts.front_hs = Some(10.);
        ts.back_hs = Some(10.);
        let mut memory = ts.allocate_memory();
        // assert!(!d.is_massive);

        let mut state = state_header.take_values().unwrap();

        // FIRST TEST -- 10 degrees on each side

        // Try marching until q_in and q_out are zero.

        let (q_in, q_out) = ts
            .march(&mut state, 10.0, 10.0, 0.0, 0.0, dt, &mut memory)
            .unwrap();

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
        let mut model = SimpleModel::default();

        /* SUBSTANCE */
        let polyurethane = add_polyurethane(&mut model);

        /* MATERIAL */
        let m1 = add_material(&mut model, polyurethane, 3. / 1000.);

        /* CONSTRUCTION */
        let mut c = Construction::new("Construction".to_string());
        c.materials.push(m1.name().clone());
        c.materials.push(m1.name().clone());
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
        let s = Surface::new("WALL".to_string(), p, c.name().clone());

        let surface = model.add_surface(s);
        let mut state_header = SimulationStateHeader::new();

        /* TEST */

        let main_dt = 3.0;
        let max_dx = m1.thickness / 7.0;
        let min_dt = 10.0;
        let d = Discretization::new(&c, &model, main_dt, max_dx, min_dt, 1., 0.).unwrap();
        let dt = main_dt / d.tstep_subdivision as Float;

        let normal = geometry3d::Vector3D::new(0., 0., 1.);
        let perimeter = 8. * l;
        let mut ts = ThermalSurface::new(
            &mut state_header,
            &model,
            &None,
            0,
            &surface,
            surface.area(),
            perimeter,
            10.,
            normal,
            &c,
            d,
        )
        .unwrap();
        ts.front_hs = Some(10.);
        ts.back_hs = Some(10.);
        let mut memory = ts.allocate_memory();
        // assert!(!d.is_massive);

        let mut state = state_header.take_values().unwrap();

        // SECOND TEST -- 10 degrees in and 30 out.
        // We expect the final q to be (30-10)/R from
        // outside to inside. I.E. q_in = (30-10)/R,
        // q_out = -(30-10)/R

        let t_front = 10.0;
        let t_back = 30.0;
        let (q_front, q_back) = ts
            .march(&mut state, t_front, t_back, 0.0, 0.0, dt, &mut memory)
            .unwrap();

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
        let mut memory = ChunkMemory {
            k,
            q,
            // These are just to put intermediate data
            aux: Matrix::new(0.0, 2, 1),
            k1: Matrix::new(0.0, 2, 1),
            k2: Matrix::new(0.0, 2, 1),
            k3: Matrix::new(0.0, 2, 1),
            k4: Matrix::new(0.0, 2, 1),
        };
        let dt = 0.01;
        rearrange_k(dt, &c, &mut memory).unwrap();

        let mut time = 0.0;
        loop {
            let temp_a = temperatures.get(0, 0).unwrap();
            let exp_temp_a = temp_a_fn(time);
            let diff_a = (temp_a - exp_temp_a).abs();

            let temp_b = temperatures.get(1, 0).unwrap();
            let exp_temp_b = temp_b_fn(time);
            let diff_b = (temp_b - exp_temp_b).abs();
            const SMOL: Float = 1e-8;
            // println!("{}, {}", diff_a, diff_b);
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

            rk4(&c, &mut memory, &mut temperatures).unwrap();

            time += dt;

            if time > 100. {
                break;
            }
        }
    }
}
