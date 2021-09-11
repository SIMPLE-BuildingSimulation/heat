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

use building_model::building::Building;
use building_model::space::Space;
use gas_properties::air;
use building_model::simulation_state::SimulationState;
use building_model::simulation_state_element::SimulationStateElement;



pub struct ThermalZone {
    /// The name of the zone
    _name: String,

    /// The position of this zone within
    /// the Thermal Model zones array
    index: usize,

    /// The position of the surfaces with which
    /// this zone is in contact in the Thermal Model
    /// surfaces array
    surface_indexes: Vec<usize>,

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
    /// by iterating the spaces in a building (so there is no mismatch).
    pub fn from_space(space: &Space, state: &mut SimulationState) -> Self {
        let space_index = space.index().unwrap();

        // Add Space Temperature state
        state.push(
            // start, by default, at 22.0 C
            SimulationStateElement::SpaceDryBulbTemperature(space_index, 22.0),
        );

        ThermalZone {
            _name: format!("ThermalZone::{}", space.name),
            index: space_index,
            volume: space.volume().unwrap(),
            surface_indexes: Vec::with_capacity(space.surfaces.len()),
        }
    }


    pub fn push_surface(&mut self, s: usize) {
        self.surface_indexes.push(s);
    }

    pub fn temperature(&self, building: &Building, state: &SimulationState) -> f64 {
        let space = &building.spaces[self.index];
        let t_index = space.dry_bulb_temperature_index().unwrap();

        if let Err(errmsg) = state[t_index].differ_only_in_value(
            SimulationStateElement::SpaceDryBulbTemperature(self.index, 123123.),
        ) {
            panic!("When retrieving Zone Temperature: {}", errmsg);
        }

        state[t_index].get_value()
    }

    /// Sets the temperature of a zone.
    pub fn set_temperature(
        &self,
        temperature: f64,
        building: &Building,
        state: &mut SimulationState,
    ) {
        let space = &building.spaces[self.index];
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
