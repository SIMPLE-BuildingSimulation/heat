
pub struct ThermalZone {
    
    /// The name of the zone
    name: String,

    /// The position of the surfaces with which
    /// this zone is in contact in the Thermal Model
    /// surfaces array
    surface_indexes: Vec< usize >,

    /// volume of the zone
    volume: f64,

    /// The position of this zone within 
    /// the Thermal Model zones array
    position: usize,

    /// The internal average temperature of 
    /// the zone's air
    temperature: f64,

    /// A value containing the amount of heat
    /// that has been accumulated as the result
    /// of the difference in temperature between the zone
    /// and the adjactent surfaces. It resets to 0.0 every time
    /// the tempearture of the zone is set
    accumulated_heat: f64,
}

impl ThermalZone{
    pub fn new(name: String, volume: f64, position: usize) -> Self{
        ThermalZone{
            name: name,
            surface_indexes: vec![],
            volume: volume,
            position: position,
            temperature: 20.0,
            accumulated_heat: 0.0,
        }
    }

    pub fn push_surface(&mut self, s: usize){        
        self.surface_indexes.push(s);        
    }

    pub fn temperature(&self)-> f64{
        self.temperature
    }

    pub fn accumulate_heat(&mut self, heat: f64){
        self.accumulated_heat += heat;
    }

    pub fn consume_heat(&mut self){
        let delta_t = self.accumulated_heat/self.mcp();
        self.temperature += delta_t;        
        self.accumulated_heat = 0.;
    }

    pub fn mcp(&self)->f64{
        let air_density = 1.225; //kg/m3
        let air_specific_heat = 1003.;//J/kg.K
        self.volume * air_density * air_specific_heat
    }
}