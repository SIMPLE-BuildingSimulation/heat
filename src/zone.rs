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

use std::rc::Rc;

use simple_model::model::SimpleModel;
use simple_model::space::Space;
use gas_properties::air;
use simple_model::simulation_state::SimulationState;
use simple_model::simulation_state_element::SimulationStateElement;



pub struct ThermalZone {
    
    /// The `Space` that this [`Thermal Zone`] represents
    reference_space: Rc<Space>,

    
    /// volume of the zone
    volume: f64,
    // The index containing the temperature of this
    // Zone in the SimulationState
    //temperature_state_index : usize,

    // The index of the state of the heating/cooling in
    // the SimulationState
    //heating_cooling_state_index: Option<usize>,

    // The index of the state of the luminaire in
    // the SimulationState
    //luminaire_state_index: Option<usize>,
}

impl ThermalZone {
    /// This function creates a new ThermalZone from a Space.
    /// It will copy the index of the space, so it should be used
    /// by iterating the spaces in a model (so there is no mismatch).
    pub fn from_space(space: &Rc<Space>, state: &mut SimulationState, space_index: usize) -> Self {
        
        let volume = *space.volume().unwrap();
        // Add Space Temperature state
        let state_index = state.push(
            // start, by default, at 22.0 C
            SimulationStateElement::SpaceDryBulbTemperature(space_index, 22.0),
        );
        space.set_dry_bulb_temperature_index(state_index);

        ThermalZone {            
            reference_space: Rc::clone(space),
            volume,
        }
    }


    pub fn temperature(&self, model: &SimpleModel, state: &SimulationState) -> f64 {
        if let Some(temp) = self.reference_space.dry_bulb_temperature(state) {
            temp
        }else{
            panic!("When retrieving temperature of Thermal Zone: reference space has no temperature field assigned.")
        }
    }

    /// Sets the temperature of a zone.
    pub fn set_temperature(
        &self,
        temperature: f64,
        model: &SimpleModel,
        state: &mut SimulationState,
    ) {
        let space = &model.spaces[self.index];
        let t_index = space.dry_bulb_temperature_index().unwrap();
        state.update_value(
            t_index,
            SimulationStateElement::SpaceDryBulbTemperature(self.index, temperature),
        );
    }

    /// Retrieves the heat capacity of the ThermalZone's air
    pub fn mcp(&self) -> f64 {
        let air_density = air::density(); //kg/m3
        let air_specific_heat = air::specific_heat(); //J/kg.K

        self.volume * air_density * air_specific_heat
    }

    
}
