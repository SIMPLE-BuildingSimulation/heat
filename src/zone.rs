use gas_properties::air;
use building_model::building::Building;
use building_model::space::Space;
use building_model::object_trait::ObjectTrait;
use building_model::building_state::BuildingState;
use building_model::building_state_element::BuildingStateElement;


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
    surface_indexes: Vec< usize >,

    /// volume of the zone
    volume: f64,

    /// The index containing the temperature of this
    /// Zone in the BuildingState
    temperature_state_index : usize,

    /// The index of the state of the heating/cooling in 
    /// the BuildingState
    heating_cooling_state_index: Option<usize>,

    /// The index of the state of the luminaire in 
    /// the BuildingState
    luminaire_state_index: Option<usize>,
    
}

impl ThermalZone{
    

    /// This function creates a new ThermalZone from a Space. 
    /// It will copy the index of the space, so it should be used
    /// by iterating the spaces in a building (so there is no mismatch).
    pub fn from_space(space: &Space, state: &mut BuildingState)->Self{

        

        // Add State
        // Add the zone to the State
        let state_index = state.push(
            // Zones start, by default, at 22.0 C
            BuildingStateElement::SpaceDryBulbTemperature(space.index(), 22.0)
        );

        

        ThermalZone {
            name : format!("ThermalZone::{}",space.name()),
            index : space.index(),
            volume : space.volume().unwrap(),
            temperature_state_index: state_index,
            surface_indexes: Vec::with_capacity(space.get_surfaces().len()),
            heating_cooling_state_index: space.get_heating_cooling_state_index(),
            luminaire_state_index: space.get_luminaires_state_index(),
        }
        
        
    }

    pub fn calc_heating_cooling_power(&self, building: &Building, state: &BuildingState)->f64{
        match self.heating_cooling_state_index{
            // Has a system... let's do something with it
            Some(i)=>{
                // Check consistency
                if let BuildingStateElement::SpaceHeatingCoolingPowerConsumption(space_index,s) = state[i]{
                    if space_index == self.index {
                        panic!("Getting Cooling / Heating for the wrong Space");
                    }
                    
                    // Get the kind of heater/cooler
                    let heater_cooler = building.get_space(space_index).unwrap().get_heating_cooling().unwrap();
                    
                    return calc_cooling_heating_power(heater_cooler, s)
                                        
                
                }else{
                    panic!("Corrupt BUildingState... incorrect BuildingStateElement found")
                }
            },
            // Does not have heating or cooling
            None => 0.0
        }
    }

    pub fn calc_lighting_power(&self, state: &BuildingState) -> f64 {
        match self.luminaire_state_index {
            // Has a system... let's do something with it
            Some(i)=>{
                // Check consistency
                if let BuildingStateElement::SpaceLightingPowerConsumption(space_index,s) = state[i]{
                    if space_index == self.index {
                        panic!("Getting Lighting for the wrong Space");
                    }                                        

                    s
                
                }else{
                    panic!("Corrupt BUildingState... incorrect BuildingStateElement found")
                }
            },
            // Does not have heating or cooling
            None => 0.0
        }
    }

    pub fn push_surface(&mut self, s: usize){        
        self.surface_indexes.push(s);        
    }

    pub fn temperature(&self, state: &BuildingState)-> f64{
        if let BuildingStateElement::SpaceDryBulbTemperature(i,v) = state[self.temperature_state_index]{
            if i != self.index {
                panic!("Incorrect index allocated for Temperature of Space '{}'", self.name);
            }
            return v;
        }else{
            panic!("Incorrect StateElement kind allocated for Temperature of Space '{}'", self.name);
        }
    }

    /*
    pub fn accumulate_heat(&mut self, heat: f64){        
        self.accumulated_heat += heat;
    }

    */
    pub fn consume_heat(&self, accumulated_heat: f64, state: &mut BuildingState){

        let delta_t = accumulated_heat/self.mcp();
        
        if let BuildingStateElement::SpaceDryBulbTemperature(i,v) = state[self.temperature_state_index]{
            if i != self.index {
                panic!("Incorrect index allocated for Temperature of Space '{}'", self.name);
            }
            state[self.temperature_state_index] = BuildingStateElement::SpaceDryBulbTemperature(i,v + delta_t)
        }else{
            panic!("Incorrect StateElement kind allocated for Temperature of Space '{}'", self.name);
        }        
    }
    
    pub fn mcp(&self)->f64{

        let air_density = air::density(); //kg/m3
        let air_specific_heat = air::specific_heat();//J/kg.K

        self.volume * air_density * air_specific_heat
    }
}