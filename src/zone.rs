use building_model::building::Building;
use building_model::object_trait::ObjectTrait;
use building_model::space::Space;
use gas_properties::air;
use simulation_state::simulation_state::SimulationState;
use simulation_state::simulation_state_element::SimulationStateElement;

use crate::heating_cooling::calc_cooling_heating_power;

pub struct ThermalZone {
    /// The name of the zone
    name: String,

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
        // Add Space Temperature state
        state.push(
            // start, by default, at 22.0 C
            SimulationStateElement::SpaceDryBulbTemperature(space.index(), 22.0),
        );

        ThermalZone {
            name: format!("ThermalZone::{}", space.name()),
            index: space.index(),
            volume: space.volume().unwrap(),
            surface_indexes: Vec::with_capacity(space.get_surfaces().len()),
        }
    }

    pub fn calc_heating_cooling_power(&self, building: &Building, state: &SimulationState) -> f64 {
        let space = building.get_space(self.index).unwrap();

        match space.get_heating_cooling_state_index() {
            // Has a system... let's do something with it
            Some(i) => {
                // Check consistency
                if let SimulationStateElement::SpaceHeatingCoolingPowerConsumption(space_index, s) =
                    state[i]
                {
                    debug_assert_eq!(space_index, self.index);

                    // Get the kind of heater/cooler
                    let heater_cooler = building
                        .get_space(space_index)
                        .unwrap()
                        .get_heating_cooling()
                        .unwrap();
                    let p = calc_cooling_heating_power(heater_cooler, s);
                    return p;
                } else {
                    panic!("Corrupt SimulationState... incorrect SimulationStateElement... found {} at index {}", state[i].to_string(), i)
                }
            }
            // Does not have heating or cooling
            None => 0.0,
        }
    }

    pub fn calc_lighting_power(&self, building: &Building, state: &SimulationState) -> f64 {
        let space = building.get_space(self.index).unwrap();

        match space.get_luminaires_state_index() {
            // Has a system... let's do something with it
            Some(i) => {
                // Check consistency
                if let SimulationStateElement::SpaceLightingPowerConsumption(space_index, s) =
                    state[i]
                {
                    if space_index != self.index {
                        panic!(
                            "Getting Lighting for the wrong Space (expected {}, found {})",
                            self.index, space_index
                        );
                    }

                    s
                } else {
                    panic!("Corrupt BUildingState... incorrect SimulationStateElement found")
                }
            }
            // Does not have heating or cooling
            None => 0.0,
        }
    }

    pub fn push_surface(&mut self, s: usize) {
        self.surface_indexes.push(s);
    }

    pub fn temperature(&self, building: &Building, state: &SimulationState) -> f64 {
        let space = building.get_space(self.index).unwrap();
        let t_index = space.get_dry_bulb_temperature_state_index().unwrap();

        if let SimulationStateElement::SpaceDryBulbTemperature(i, v) = state[t_index] {
            if i != self.index {
                panic!(
                    "Incorrect index allocated for Temperature of Space '{}'",
                    self.name
                );
            }

            return v;
        } else {
            panic!(
                "Incorrect StateElement kind allocated for Temperature of Space '{}'... found {}",
                self.name,
                state[t_index].to_string()
            );
        }
    }

    pub fn consume_heat(
        &self,
        accumulated_heat: f64,
        building: &Building,
        state: &mut SimulationState,
    ) {
        let delta_t = accumulated_heat / self.mcp();

        let space = building.get_space(self.index).unwrap();
        let t_index = space.get_dry_bulb_temperature_state_index().unwrap();

        if let SimulationStateElement::SpaceDryBulbTemperature(i, v) = state[t_index] {
            if i != self.index {
                panic!(
                    "Incorrect index allocated for Temperature of Space '{}'",
                    self.name
                );
            }
            state[t_index] = SimulationStateElement::SpaceDryBulbTemperature(i, v + delta_t)
        } else {
            panic!(
                "Incorrect StateElement kind allocated for Temperature of Space '{}'",
                self.name
            );
        }
    }

    pub fn set_temperature(
        &self,
        temperature: f64,
        building: &Building,
        state: &mut SimulationState,
    ) {
        let space = building.get_space(self.index).unwrap();
        let t_index = space.get_dry_bulb_temperature_state_index().unwrap();

        if let SimulationStateElement::SpaceDryBulbTemperature(i, _) = state[t_index] {
            if i != self.index {
                panic!(
                    "Incorrect index allocated for Temperature of Space '{}'",
                    self.name
                );
            }

            state[t_index] = SimulationStateElement::SpaceDryBulbTemperature(i, temperature)
        } else {
            panic!(
                "Incorrect StateElement kind allocated for Temperature of Space '{}'",
                self.name
            );
        }
    }

    /// Retrieves the heat capacity of the ThermalZone's air
    pub fn mcp(&self) -> f64 {
        let air_density = air::density(); //kg/m3
        let air_specific_heat = air::specific_heat(); //J/kg.K

        self.volume * air_density * air_specific_heat
    }
}
