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
use simple_model::simulation_state::SimulationState;
use std::rc::Rc;
use simple_model::hvac::cast_hvac;
use simple_model::hvac::{HVACKind, HVAC};
use simple_model::hvac::ideal_heater_cooler::IdealHeaterCooler;
use simple_model::hvac::electric_heater::ElectricHeater;

/// Retrieves a `Vec<(usize, f64)>` containing the amount of heat (the `f64` in W) going into 
/// each space (of index `usize`)
pub fn calc_cooling_heating_power(system: &Rc<dyn HVAC>, state: &SimulationState ) -> Vec<(usize,f64)> {
    
    match system.kind() {
        HVACKind::IdealHeaterCooler => {
            let a = &**system;
            let system = cast_hvac::<IdealHeaterCooler>(a).unwrap();
            let mut ret = Vec::new();
            for space in &system.target_spaces{
                let index = space.index().unwrap();
                let consumption_power = system.heating_cooling_consumption(state).expect("HVAC has not heating/cooling state");
                ret.push((*index,consumption_power))
            }
            ret
            
        },
        HVACKind::ElectricHeater => {
            let a = &**system;
            let system = cast_hvac::<ElectricHeater>(a).unwrap();
            let mut ret = Vec::new();
            if let Ok(space) = system.target_space(){
                let index = space.index().unwrap();
                let consumption_power = system.heating_cooling_consumption(state).expect("HVAC has not heating/cooling state");
                ret.push((*index,consumption_power))
            }
            ret
            
        }
    }
}
