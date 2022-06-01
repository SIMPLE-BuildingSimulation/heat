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

/// Represents a border condition of between a Surface
/// and a Zone or the exterior
#[derive(Debug, Clone, Copy)]
pub struct Environment {
    
    /// The dry bulb temperature of the air, in $`C`$
    pub air_temperature: Float,

    /// The wind speed, in m/2
    pub air_speed: Float,

    /// The incident Infrared Irradiance, in $`W/m^2`$
    pub ir_irrad: Float,

    /// The incident Solar Irradiance, in $`W/m^2`$
    pub solar_radiation: Float,

    /// The environmental emmisivity, used for calculating incident IR
    /// irradiance, if needed
    pub env_emmisivity: Float,
    
}


impl std::default::Default for Environment {
    fn default()->Self{
        const DEFAULT_ENV_EMMISIVITY: Float = 1.;
        const DEFAULT_AIR_TEMP: Float = 22.;
        Self { 
            air_temperature: DEFAULT_AIR_TEMP, 
            air_speed: 0., 
            ir_irrad: crate::SIGMA * DEFAULT_ENV_EMMISIVITY * ( DEFAULT_AIR_TEMP + 273.15 ).powi(4), 
            solar_radiation: 0.,
            env_emmisivity: DEFAULT_ENV_EMMISIVITY 
        }
    }
}

impl Environment {

    pub fn get_hs(&self)->Float{
        // dbg!("Calculate back Rs");
        10.
        // if self.air_speed > 0.{
        //     // 36.34359273
        //     // 32.6
        // }else{
        //     // 5.
        //     // 6.8

        // }
    }

    

}

