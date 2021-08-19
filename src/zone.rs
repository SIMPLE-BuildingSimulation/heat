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

    // /// Calculates the amount of heating or cooling being delivered to the
    // /// thermal zone based on the Energy Consumption registered in the SimulationState    
    // pub fn calc_heating_cooling_power(&self, building: &Building, _state: &SimulationState) -> f64 {
    //     let _space = &building.spaces[self.index];

        // 0.0
        // match space.get_heating_cooling_state_index() {
        //     // Has a system... let's do something with it
        //     Some(i) => {
        //         // Check consistency
        //         if let SimulationStateElement::SpaceHeatingCoolingPowerConsumption(space_index, s) =
        //             state[i]
        //         {
        //             debug_assert_eq!(space_index, self.index);

        //             // Get the kind of heater/cooler
        //             let heater_cooler = building.spaces[space_index]
        //                 .get_heating_cooling()
        //                 .unwrap();

        //             // return
        //             calc_cooling_heating_power(heater_cooler, s)
        //         } else {
        //             panic!("Corrupt SimulationState... incorrect SimulationStateElement... found {} at index {}", state[i].to_string(), i)
        //         }
        //     }
        //     // Does not have heating or cooling
        //     None => 0.0,
        // }
    // }

    pub fn calc_lighting_power(&self, building: &Building, _state: &SimulationState) -> f64 {
        let _space = &building.spaces[self.index];
        0.0
        // match space.get_luminaires_state_index() {
        //     // Has a system... let's do something with it
        //     Some(i) => {
        //         // Check consistency
        //         if let SimulationStateElement::SpaceLightingPowerConsumption(space_index, s) =
        //             state[i]
        //         {
        //             if space_index != self.index {
        //                 panic!(
        //                     "Getting Lighting for the wrong Space (expected {}, found {})",
        //                     self.index, space_index
        //                 );
        //             }

        //             s
        //         } else {
        //             panic!("Corrupt BUildingState... incorrect SimulationStateElement found")
        //         }
        //     }
        //     // Does not have heating or cooling
        //     None => 0.0,
        // }
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

    /// Calculates the heat added to a [`ThermalZone`] by lightng, people and appliances
    pub fn get_current_internal_heat_loads(
        &self,
        building: &Building,
        state: &SimulationState,
    ) -> f64 {
        
        // lighting        
        let lighting = self.calc_lighting_power(building, state);

        // people
        let people = 0.;

        // Appliances
        let appliances = 0.;

        // return
        lighting + people + appliances 
    }
}
