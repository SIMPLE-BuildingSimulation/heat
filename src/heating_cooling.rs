use building_model::heating_cooling::{HeaterCooler, HeatingCoolingState, HeatingCoolingKind};

pub fn calc_cooling_heating_power(system: &HeaterCooler, state: HeatingCoolingState)->f64{
    
    if let HeatingCoolingState::Off = state {
        return 0.0
    }
    
    match system.get_kind() {
        HeatingCoolingKind::IdealHeaterCooler => {
            match state {
                HeatingCoolingState::Off       => 0.0,
                HeatingCoolingState::Cooling(p)=> -p,
                HeatingCoolingState::Heating(p)=>  p,
            }
        },
        HeatingCoolingKind::ElectricHeating => {
            match state {
                HeatingCoolingState::Off       => 0.0,
                HeatingCoolingState::Cooling(_)=> panic!("Electric Heating cannot be cooling!"),
                HeatingCoolingState::Heating(p)=>  p,
            }
        }
    }
}