use building_model::heating_cooling::{HeaterCooler, HeatingCoolingKind};

pub fn calc_cooling_heating_power(system: &HeaterCooler, consumption_power: f64)->f64{
        
    
    
    match system.get_kind() {
        HeatingCoolingKind::IdealHeaterCooler => consumption_power,
        
        HeatingCoolingKind::ElectricHeating => {
            if consumption_power >= 0.0 {
                return consumption_power
            }else {
                panic!("Electric Heating cannot be cooling!")
            }            
        }
    }
}