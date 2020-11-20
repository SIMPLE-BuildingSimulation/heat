//use building_model::construction::Construction;
use building_model::building::Building;
use building_model::construction::Construction;

/// Calculates the R-value of the Construction. 
/// # Parameters:
/// * materials: A reference to the vector containing the materials from which this object was built from
pub fn r_value(building: &Building, construction: &Construction)->Result<f64,String>{

    //let construction = building.get_construction(construction_index).unwrap();
    
    //let materials = building.get_materials();
    //let substances = building.get_substances();

    let mut r = 0.0;
    
    for material_index in construction.layers() {
        let material = building.get_material(*material_index).unwrap();
        let substance_index = material.get_substance_index().unwrap();
        let substance = building.get_substance(substance_index).unwrap();        
        let lambda = substance.thermal_conductivity().unwrap();

        r += material.thickness().unwrap()/lambda;
    }

    return Ok(r);
}
