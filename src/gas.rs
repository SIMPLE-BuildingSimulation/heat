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
use polynomial::{poly, Polynomial};

/// A structure containing the data that will describe the thermal
/// behaviour of a gas.
struct Gas {
    /// The thermal conductivity (W/m.K) as a function of the
    /// temperature (in K)
    thermal_conductivity: Polynomial,

    /// The dynamic viscosity (N.s/m2) as a function of the
    /// temperature (in K)
    dynamic_viscosity: Polynomial,

    /// The specific heat capacity (J/kg.K) as a function of the
    /// temperature (in K)
    heat_capacity: Polynomial,

    /// THe Molecular Mass (kg/Mol)
    mass: Float,
}

/// Transforms C into K
fn in_kelvin(t: Float) -> Float {
    t + 273.15
}

impl Gas {
    /// Derives the Thermal Conductivity from a certain Temperature (in C)
    pub fn thermal_conductivity(&self, temp: Float) -> Float {
        self.thermal_conductivity.eval(in_kelvin(temp))
    }

    /// Derives the Dynamic Viscosity from a certain Temperature (in C)
    pub fn dynamic_viscosity(&self, temp: Float) -> Float {
        self.dynamic_viscosity.eval(in_kelvin(temp))
    }

    /// Derives the Soecific Heat Capacity from a certain Temperature (in C)
    pub fn heat_capacity(&self, temp: Float) -> Float {
        self.heat_capacity.eval(in_kelvin(temp))
    }

    /// Retreives the Molecular Mass
    pub fn mass(&self) -> Float {
        self.mass
    }
}

fn air() -> Gas {
    Gas {
        thermal_conductivity: poly![2.873e-3, 7.760e-5],
        dynamic_viscosity: poly![3.723e-6, 4.94e-8],
        heat_capacity: poly![1002.7370, 1.2324e-2],
        mass: 28.97,
    }
}

fn argon() -> Gas {
    Gas {
        thermal_conductivity: poly![2.285e-3, 5.149e-5],
        dynamic_viscosity: poly![3.379e-6, 6.451e-8],
        heat_capacity: poly![521.9285],
        mass: 39.948,
    }
}

fn krypton() -> Gas {
    Gas {
        thermal_conductivity: poly![9.443e-4, 2.826e-5],
        dynamic_viscosity: poly![2.213e-6, 7.777e-8],
        heat_capacity: poly![248.0907],
        mass: 83.8,
    }
}

fn xenon() -> Gas {
    Gas {
        thermal_conductivity: poly![4.538e-4, 1.723e-5],
        dynamic_viscosity: poly![1.069e-6, 7.414e-8],
        heat_capacity: poly![158.3397],
        mass: 131.30,
    }
}

/***********/
/* TESTING */
/***********/

#[cfg(test)]
mod testing {
    use super::*;

    fn check_value(a: Float, b: Float) -> Result<(), String> {
        let err = (a - b).abs() / a.abs();
        if err > 1e-2 {
            //1% error max
            return Err(format!("a = {} | b = {} | err = {}", a, b, err));
        }
        Ok(())
    }

    #[test]
    fn test_thermal_conductivity() {
        check_value(0.0241, crate::gas::air().thermal_conductivity(0.)).unwrap();
        check_value(0.0248, crate::gas::air().thermal_conductivity(10.)).unwrap();

        check_value(0.0163, crate::gas::argon().thermal_conductivity(0.)).unwrap();
        check_value(0.0169, crate::gas::argon().thermal_conductivity(10.)).unwrap();

        check_value(0.0087, crate::gas::krypton().thermal_conductivity(0.)).unwrap();
        check_value(0.0089, crate::gas::krypton().thermal_conductivity(10.)).unwrap();

        check_value(0.0052, crate::gas::xenon().thermal_conductivity(0.)).unwrap();
        check_value(0.0053, crate::gas::xenon().thermal_conductivity(10.)).unwrap();
    }

    #[test]
    fn test_dynamic_viscosity() {
        check_value(1.722e-5, crate::gas::air().dynamic_viscosity(0.)).unwrap();
        check_value(1.771e-5, crate::gas::air().dynamic_viscosity(10.)).unwrap();

        check_value(2.1e-5, crate::gas::argon().dynamic_viscosity(0.)).unwrap();
        check_value(2.165e-5, crate::gas::argon().dynamic_viscosity(10.)).unwrap();

        check_value(2.346e-5, crate::gas::krypton().dynamic_viscosity(0.)).unwrap();
        check_value(2.423e-5, crate::gas::krypton().dynamic_viscosity(10.)).unwrap();

        check_value(2.132e-5, crate::gas::xenon().dynamic_viscosity(0.)).unwrap();
        check_value(2.206e-5, crate::gas::xenon().dynamic_viscosity(10.)).unwrap();
    }

    #[test]
    fn test_heat_capacity() {
        check_value(1006.1034, crate::gas::air().heat_capacity(0.)).unwrap();
        check_value(1006.2265, crate::gas::air().heat_capacity(10.)).unwrap();

        check_value(521.9285, crate::gas::argon().heat_capacity(0.)).unwrap();
        check_value(521.9285, crate::gas::argon().heat_capacity(10.)).unwrap();

        check_value(248.0907, crate::gas::krypton().heat_capacity(0.)).unwrap();
        check_value(248.0907, crate::gas::krypton().heat_capacity(10.)).unwrap();

        check_value(158.3397, crate::gas::xenon().heat_capacity(0.)).unwrap();
        check_value(158.3397, crate::gas::xenon().heat_capacity(10.)).unwrap();
    }

    #[test]
    fn test_mass() {
        check_value(28.97, crate::gas::air().mass()).unwrap();
        check_value(39.948, crate::gas::argon().mass()).unwrap();
        check_value(83.80, crate::gas::krypton().mass()).unwrap();
        check_value(131.3, crate::gas::xenon().mass()).unwrap();
    }
}
