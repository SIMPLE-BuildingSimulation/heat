use communication_protocols::simulation_model::SimulationModel;
use thermal::model::ThermalModel;
use thermal::Float;

use calendar::Date;
use schedule::ScheduleConstant;
use validate::*;
use weather::SyntheticWeather;

use simple_model::{SimulationStateElement, HVAC};
use simple_test_models::{SingleZoneTestBuildingOptions, TestMat, get_single_zone_test_building};

/// A single-zone test model with walls assumed to have
/// no mass. It has a closed solution, which is nice.
///
/// There is no sun.
#[derive(Default)]
struct SingleZoneTestModel {
    /// volume of the zone (m3)
    zone_volume: Float,

    /// Facade area (m2)
    surface_area: Float,

    /// the R-value of the facade
    facade_r: Float,

    /// Infiltration rate (m3/s)
    infiltration_rate: Float,

    /// Heating power (Watts)
    heating_power: Float,

    /// Lighting power (Watts)
    lighting_power: Float,

    /// Temperature outside of the zone
    temp_out: Float,

    /// Temperature at the beginning
    temp_start: Float,
}

impl SingleZoneTestModel {
    fn get_closed_solution(&self) -> Box<impl Fn(Float) -> Float> {
        // heat balance in the form
        // of C*dT/dt = A - B*T
        let air = thermal::gas::Gas::air();
        let rho = air.density(22. + 273.15); //kg/m3
        let cp = air.heat_capacity(22. + 273.15); //J/kg.K
        let u = 1. / self.facade_r;

        let c = self.zone_volume * rho * cp;

        let a = self.heating_power
            + self.lighting_power
            + self.temp_out * u * self.surface_area
            + self.infiltration_rate * rho * cp * self.temp_out;

        let b = u * self.surface_area + rho * self.infiltration_rate * cp;

        let k1 = self.temp_start - a / b;

        let f = move |t: Float| -> Float { a / b + k1 * (-b * t / c).exp() };

        Box::new(f)
    }
}

fn march_with_window() ->(Vec<Float>, Vec<Float>) {
    let surface_area = 4.;
    let window_area = 1.;
    let zone_volume = 40.;

    let (simple_model, mut state_header) = get_single_zone_test_building(
        // &mut state,
        &SingleZoneTestBuildingOptions {
            zone_volume,
            surface_area,
            window_area,
            construction: vec![TestMat::Polyurethane(0.02)],
            emmisivity: 0.0,
            ..Default::default()
        },
    );

    // Finished model the SimpleModel
    let n: usize = 6;
    let main_dt = 60. * 60. / n as Float;
    let thermal_model = ThermalModel::new(&simple_model, &mut state_header, n).unwrap();

    let mut state = state_header.take_values().unwrap();

    // MAP THE STATE
    // model.map_simulation_state(&mut state).unwrap();

    // START TESTING.
    let hs_front = simple_model.surfaces[0]
        .front_convection_coefficient(&state)
        .unwrap();
    let hs_back = simple_model.surfaces[0]
        .back_convection_coefficient(&state)
        .unwrap();
    let r = thermal_model.surfaces[0].discretization.r_value() + 1. / hs_front + 1. / hs_back;

    // Initial T of the zone
    let t_start = thermal_model.zones[0]
        .reference_space
        .dry_bulb_temperature(&state)
        .unwrap();

    let t_out: Float = 30.0; // T of surroundings

    let mut weather = SyntheticWeather::new();
    weather.dry_bulb_temperature = Box::new(ScheduleConstant::new(t_out));

    let dt = main_dt;

    let mut date = Date {
        day: 1,
        hour: 0.0,
        month: 1,
    };

    // test model
    let tester = SingleZoneTestModel {
        zone_volume,
        surface_area, // the window is a hole on the wall... does not add area
        facade_r: r,
        temp_out: t_out,
        temp_start: t_start,
        ..SingleZoneTestModel::default()
    };
    let exp_fn = tester.get_closed_solution();

    // March:
    let n = 80;
    let mut exp = Vec::with_capacity(n);
    let mut found = Vec::with_capacity(n);
    for i in 0..n {
        let time = (i as Float) * dt;
        date.add_seconds(time);

        let found_v = thermal_model.zones[0]
            .reference_space
            .dry_bulb_temperature(&state)
            .unwrap();

        thermal_model
            .march(date, &weather, &simple_model, &mut state)
            .unwrap();

        // Get exact solution.
        let exp_v = exp_fn(time);
        exp.push(exp_v);
        found.push(found_v);
    }
    (exp,found)
}

fn very_simple_march() -> (Vec<Float>, Vec<Float>) {
    let zone_volume = 40.;
    let surface_area = 4.;
    let (simple_model, mut state_header) = get_single_zone_test_building(
        // &mut state,
        &SingleZoneTestBuildingOptions {
            zone_volume,
            surface_area,
            construction: vec![TestMat::Polyurethane(0.02)],
            emmisivity: 0.0,
            ..Default::default()
        },
    );

    let n: usize = 60;
    let main_dt = 60. * 60. / n as Float;
    let thermal_model = ThermalModel::new(&simple_model, &mut state_header, n).unwrap();

    let mut state = state_header.take_values().unwrap();

    let hs_front = simple_model.surfaces[0]
        .front_convection_coefficient(&state)
        .unwrap();
    let hs_back = simple_model.surfaces[0]
        .back_convection_coefficient(&state)
        .unwrap();

    let r = thermal_model.surfaces[0].discretization.r_value() + 1. / hs_front + 1. / hs_back;

    // Initial T of the zone
    let t_start = thermal_model.zones[0]
        .reference_space
        .dry_bulb_temperature(&state)
        .unwrap();

    let t_out: Float = 30.0; // T of surroundings

    // test model
    let tester = SingleZoneTestModel {
        zone_volume,
        surface_area,
        facade_r: r,
        temp_out: t_out,
        temp_start: t_start,
        ..SingleZoneTestModel::default()
    };
    let exp_fn = tester.get_closed_solution();

    let mut weather = SyntheticWeather::new();
    weather.dry_bulb_temperature = Box::new(ScheduleConstant::new(t_out));

    let mut date = Date {
        day: 1,
        hour: 0.0,
        month: 1,
    };

    let n = 1000;
    let mut exp = Vec::with_capacity(n);
    let mut found = Vec::with_capacity(n);
    for i in 0..1000 {
        let time = (i as Float) * main_dt;
        date.add_seconds(time);

        let found_v = thermal_model.zones[0]
            .reference_space
            .dry_bulb_temperature(&state)
            .unwrap();

        thermal_model
            .march(date, &weather, &simple_model, &mut state)
            .unwrap();

        // Get exact solution.
        let exp_v = exp_fn(time);

        exp.push(exp_v);
        found.push(found_v);
    }

    return (exp, found);
}

fn march_with_window_and_luminaire() -> (Vec<Float>, Vec<Float>) {
    let surface_area = 4.;
    let zone_volume = 40.;
    let lighting_power = 100.;

    let (simple_model, mut state_header) = get_single_zone_test_building(
        // &mut state,
        &SingleZoneTestBuildingOptions {
            zone_volume,
            surface_area,
            lighting_power,
            construction: vec![TestMat::Polyurethane(0.02)],
            emmisivity: 0.0,
            ..Default::default()
        },
    );

    // Finished model the SimpleModel

    let n: usize = 20;
    let main_dt = 60. * 60. / n as Float;
    let thermal_model = ThermalModel::new(&simple_model, &mut state_header, n).unwrap();

    let mut state = state_header.take_values().unwrap();

    
    // turn the lights on
    let lum_state_i = simple_model.luminaires[0]
        .power_consumption_index()
        .unwrap();        
    state[lum_state_i] = lighting_power;

    // START TESTING.

    

    let hs_front = simple_model.surfaces[0]
        .front_convection_coefficient(&state)
        .unwrap();
    let hs_back = simple_model.surfaces[0]
        .back_convection_coefficient(&state)
        .unwrap();
    let r = thermal_model.surfaces[0].discretization.r_value() + 1./hs_front + 1./hs_back;

    // Initial T of the zone
    let t_start = 22.;
    
    thermal_model.zones[0]
        .reference_space
        .set_dry_bulb_temperature(&mut state, t_start);
        

    let t_out: Float = 30.0; // T of surroundings
    
    // test model
    let tester = SingleZoneTestModel {
        zone_volume,
        surface_area, // the window is a hole on the wall... does not add area
        lighting_power,
        facade_r: r,
        temp_out: t_out,
        temp_start: t_start,
        ..SingleZoneTestModel::default()
    };
    let exp_fn = tester.get_closed_solution();

    let mut weather = SyntheticWeather::new();
    weather.dry_bulb_temperature = Box::new(ScheduleConstant::new(t_out));

    let dt = main_dt; // / model.dt_subdivisions() as Float;

    let mut date = Date {
        day: 1,
        hour: 0.0,
        month: 1,
    };

    // March:
    let n = 800;
    let mut exp = Vec::with_capacity(n);
    let mut found = Vec::with_capacity(n);
    for i in 0..n {
        
        
        let time = (i as Float) * dt;
        date.add_seconds(time);

        let found_v = thermal_model.zones[0]
            .reference_space
            .dry_bulb_temperature(&state)
            .unwrap();

        thermal_model
            .march(date, &weather, &simple_model, &mut state)
            .unwrap();

        // Get exact solution.
        let exp_v = exp_fn(time);

        exp.push(exp_v);
        found.push(found_v);

        
    }

    (exp, found)
}

fn march_with_window_and_heater() -> (Vec<Float>, Vec<Float>) {
    let surface_area = 4.;
    let zone_volume = 40.;
    let heating_power = 100.;

    let (simple_model, mut state_header) = get_single_zone_test_building(
        // &mut state,
        &SingleZoneTestBuildingOptions {
            zone_volume,
            surface_area,
            heating_power,
            construction: vec![TestMat::Polyurethane(0.02)],
            emmisivity: 0.0,
            ..Default::default()
        },
    );

    // Finished model the SimpleModel

    let n: usize = 20;
    let main_dt = 60. * 60. / n as Float;
    let thermal_model = ThermalModel::new(&simple_model, &mut state_header, n).unwrap();

    let mut state = state_header.take_values().unwrap();
    // MAP THE STATE
    // model.map_simulation_state(&mut state).unwrap();

    // turn the heater on
    if let HVAC::ElectricHeater(heater) = &simple_model.hvacs[0] {
        let hvac_state_i = heater.heating_cooling_consumption_index().unwrap();
        state[hvac_state_i] = heating_power;
    }

    // START TESTING.
    // assert!(!model.surfaces[0].is_massive());

    let hs_front = simple_model.surfaces[0]
        .front_convection_coefficient(&state)
        .unwrap();
    let hs_back = simple_model.surfaces[0]
        .back_convection_coefficient(&state)
        .unwrap();
    let r = thermal_model.surfaces[0].discretization.r_value() + 1./hs_front + 1./hs_back;

    // Initial T of the zone
    let t_start = thermal_model.zones[0]
        .reference_space
        .dry_bulb_temperature(&state)
        .unwrap();
    let t_out: Float = 30.0; // T of surroundings

    // test model
    let tester = SingleZoneTestModel {
        zone_volume,
        surface_area, // the window is a hole on the wall... does not add area
        heating_power,
        facade_r: r,
        temp_out: t_out,
        temp_start: t_start,
        ..SingleZoneTestModel::default()
    };
    let exp_fn = tester.get_closed_solution();

    let mut weather = SyntheticWeather::new();
    weather.dry_bulb_temperature = Box::new(ScheduleConstant::new(t_out));

    let dt = main_dt; // / model.dt_subdivisions() as Float;

    let mut date = Date {
        day: 1,
        hour: 0.0,
        month: 1,
    };

    // March:
    let n = 800;
    let mut exp = Vec::with_capacity(n);
    let mut found = Vec::with_capacity(n);
    for i in 0..n {
        let time = (i as Float) * dt;
        date.add_seconds(time);

        let found_v = thermal_model.zones[0]
            .reference_space
            .dry_bulb_temperature(&state)
            .unwrap();

        thermal_model
            .march(date, &weather, &simple_model, &mut state)
            .unwrap();

        // Get exact solution.
        let exp_v = exp_fn(time);

        exp.push(exp_v);
        found.push(found_v);
    }
    (exp, found)
}

fn march_with_window_heater_and_infiltration() -> (Vec<Float>, Vec<Float>)  {
    let surface_area = 4.;
    let zone_volume = 40.;
    let heating_power = 10.;
    let infiltration_rate = 0.1;
    let t_out: Float = 30.0; // T of surroundings

    let (simple_model, mut state_header) = get_single_zone_test_building(
        // &mut state,
        &SingleZoneTestBuildingOptions {
            zone_volume,
            surface_area,
            heating_power,
            infiltration_rate,
            emmisivity: 0.0,
            construction: vec![TestMat::Polyurethane(0.02)],
            ..Default::default()
        },
    );

    // Finished model the SimpleModel

    let n: usize = 20;
    let main_dt = 60. * 60. / n as Float;
    let thermal_model = ThermalModel::new(&simple_model, &mut state_header, n).unwrap();

    // Set infiltration
    let inf_vol_index = state_header.push(
        SimulationStateElement::SpaceInfiltrationVolume(0),
        infiltration_rate,
    );
    simple_model.spaces[0].set_infiltration_volume_index(inf_vol_index);
    let inf_temp_index = state_header.push(
        SimulationStateElement::SpaceInfiltrationTemperature(0),
        t_out,
    );
    simple_model.spaces[0].set_infiltration_temperature_index(inf_temp_index);

    // MAP THE STATE

    let mut state = state_header.take_values().unwrap();

    // turn the heater on
    if let HVAC::ElectricHeater(heater) = &simple_model.hvacs[0] {
        let hvac_state_i = heater.heating_cooling_consumption_index().unwrap();
        state[hvac_state_i] = heating_power;
    }

    // START TESTING.

    let hs_front = simple_model.surfaces[0]
        .front_convection_coefficient(&state)
        .unwrap();
    let hs_back = simple_model.surfaces[0]
        .back_convection_coefficient(&state)
        .unwrap();
    let r = thermal_model.surfaces[0].discretization.r_value() + 1./hs_front + 1./hs_back;

    // Initial T of the zone
    let t_start = thermal_model.zones[0]
        .reference_space
        .dry_bulb_temperature(&state)
        .unwrap();

    // test model
    let tester = SingleZoneTestModel {
        zone_volume,
        surface_area, // the window is a hole on the wall... does not add area
        heating_power,
        facade_r: r,
        temp_out: t_out,
        temp_start: t_start,
        infiltration_rate,
        ..SingleZoneTestModel::default()
    };
    let exp_fn = tester.get_closed_solution();

    let mut weather = SyntheticWeather::new();
    weather.dry_bulb_temperature = Box::new(ScheduleConstant::new(t_out));

    let dt = main_dt; // / model.dt_subdivisions() as Float;

    let mut date = Date {
        day: 1,
        hour: 0.0,
        month: 1,
    };

    // March:
    let n = 22;
    let mut exp = Vec::with_capacity(n);
    let mut found = Vec::with_capacity(n);
    for i in 0..n {
        let time = (i as Float) * dt;
        date.add_seconds(time);

        let found_v = thermal_model.zones[0]
            .reference_space
            .dry_bulb_temperature(&state)
            .unwrap();

        thermal_model
            .march(date, &weather, &simple_model, &mut state)
            .unwrap();

        // Get exact solution.
        let exp_v = exp_fn(time);

        exp.push(exp_v);
        found.push(found_v);
    }
    (exp, found)
}


fn march_solar_radiation_massive() -> (Vec<Float>, Vec<Float>)  {
    let surface_area = 20. * 3.;
    let zone_volume = 600.;

    let (simple_model, mut state_header) = get_single_zone_test_building(
        // &mut state,
        &SingleZoneTestBuildingOptions {
            zone_volume,
            surface_area,
            construction: vec![TestMat::Concrete(0.2)],
            emmisivity: 0.9,
            ..Default::default()
        },
    );

    // Finished model the SimpleModel

    let n: usize = 20;
    // let main_dt = 60. * 60. / n as Float;
    let thermal_model = ThermalModel::new(&simple_model, &mut state_header, n).unwrap();
    
    let mut state = state_header.take_values().unwrap();

    let cols = validate::from_csv("./tests/solar_radiation_massive/eplusout.csv", &[3, 9,10,11, 13]);
    let incident_solar_radiation = &cols[0]; //3
    let _indoor_thermal_heat_gain = &cols[1]; //9
    let outdoor_temp = &cols[2]; //10
    let outdoor_thermal_heat_gain = &cols[3]; //11
    let exp_zone_air_temp = &cols[4]; //13

    // Set initial temperature
    simple_model.spaces[0].set_dry_bulb_temperature(&mut state, exp_zone_air_temp[0]);

    let mut date = Date {
        month: 1,
        day: 1,
        hour: 0.0,
    };
    let n  = outdoor_temp.len();
    let mut exp = Vec::with_capacity(n);
    let mut found = Vec::with_capacity(n);
    for i in 0..n {        
        // Get zone's temp
        let found_temp = simple_model.spaces[0].dry_bulb_temperature(&state).unwrap();
        let exp_temp = exp_zone_air_temp[i];
        if i > 300 { // skip warmup
            exp.push(exp_temp);
            found.push(found_temp);                       
        }

        // Set outdoor temp
        let mut weather = SyntheticWeather::new();
        weather.dry_bulb_temperature = Box::new(ScheduleConstant::new(outdoor_temp[i]));

        let surface = &simple_model.surfaces[0];

        // Set Solar Radiation
        if incident_solar_radiation[i] > 1e-5{
            surface.set_back_incident_solar_irradiance(&mut state, incident_solar_radiation[i]);
        }
        surface.set_back_incident_solar_irradiance(&mut state, incident_solar_radiation[i]);

        // Set IR radiation
        // thermal_heat_gain[i] / surface_area/  0.9  = crate::SIGMA  ( tout.powi(4) - ts.powi(4) );
        let ts = surface.first_node_temperature(&state).unwrap();
        surface.set_back_ir_irradiance(&mut state, outdoor_thermal_heat_gain[i] / surface_area / 0.9 + thermal::SIGMA * (ts + 273.15).powi(4));
        // let ts = surface.last_node_temperature(&state).unwrap();
        // let v = indoor_thermal_heat_gain[i] / surface_area / 0.9 + thermal::SIGMA * (ts + 273.15).powi(4);
        // surface.set_front_ir_irradiance(&mut state, v);
        surface.set_front_ir_irradiance(&mut state, thermal::SIGMA * (exp_temp + 273.15 as Float).powi(4));

        // March
        thermal_model
            .march(date, &weather, &simple_model, &mut state)
            .unwrap();

        // Advance
        date.add_hours(1. / n as Float);
        

    }
    (exp, found)
}

fn march_solar_radiation_mixed_mass() -> (Vec<Float>, Vec<Float>)  {
    let surface_area = 20. * 3.;
    let zone_volume = 600.;

    let (simple_model, mut state_header) = get_single_zone_test_building(
        // &mut state,
        &SingleZoneTestBuildingOptions {
            zone_volume,
            surface_area,
            construction: vec![TestMat::Polyurethane(0.02), TestMat::Concrete(0.2), TestMat::Polyurethane(0.02)],
            emmisivity: 0.9,
            ..Default::default()
        },
    );

    // Finished model the SimpleModel

    let n: usize = 20;
    // let main_dt = 60. * 60. / n as Float;
    let thermal_model = ThermalModel::new(&simple_model, &mut state_header, n).unwrap();
    
    let mut state = state_header.take_values().unwrap();

    let cols = validate::from_csv("./tests/solar_radiation_mixed/eplusout.csv", &[3, 9,10,11, 13]);
    let incident_solar_radiation = &cols[0];
    let _indoor_thermal_heat_gain = &cols[1];    
    let outdoor_temp = &cols[2];
    let outdoor_thermal_heat_gain = &cols[3];
    let exp_zone_air_temp = &cols[4];

    // Set initial temperature
    simple_model.spaces[0].set_dry_bulb_temperature(&mut state, exp_zone_air_temp[0]);

    let mut date = Date {
        month: 1,
        day: 1,
        hour: 0.0,
    };
    let n  = outdoor_temp.len();
    let mut exp = Vec::with_capacity(n);
    let mut found = Vec::with_capacity(n);
    for i in 0..n {        
        // Get zone's temp
        let found_temp = simple_model.spaces[0].dry_bulb_temperature(&state).unwrap();
        let exp_temp = exp_zone_air_temp[i];
        if i > 300 { // skip warmup
            exp.push(exp_temp);
            found.push(found_temp);                       
        }

        // Set outdoor temp
        let mut weather = SyntheticWeather::new();
        weather.dry_bulb_temperature = Box::new(ScheduleConstant::new(outdoor_temp[i]));

        let surface = &simple_model.surfaces[0];

        // Set Solar Radiation
        surface.set_back_incident_solar_irradiance(&mut state, incident_solar_radiation[i]);

        // Set IR radiation
        // thermal_heat_gain[i] / surface_area/  0.9  = crate::SIGMA  ( tout.powi(4) - ts.powi(4) );
        let ts = surface.first_node_temperature(&state).unwrap();
        surface.set_back_ir_irradiance(&mut state, outdoor_thermal_heat_gain[i] / surface_area / 0.9 + thermal::SIGMA * (ts + 273.15).powi(4));
        let _ts = surface.last_node_temperature(&state).unwrap();                
        // let v = indoor_thermal_heat_gain[i] / surface_area / 0.9 + thermal::SIGMA * (ts + 273.15).powi(4);
        // surface.set_front_ir_irradiance(&mut state, v);
        surface.set_front_ir_irradiance(&mut state, thermal::SIGMA * (exp_temp + 273.15 as Float).powi(4));


        // March
        thermal_model
            .march(date, &weather, &simple_model, &mut state)
            .unwrap();

        // Advance
        date.add_hours(1. / n as Float);
        

    }
    (exp, found)
}



#[test]
fn validate() {
    
    let mut validations = Validations::new("report.md", ".");
    
    let (expected, found) = very_simple_march();
    let v = validate::SeriesValidator {
        x_label: Some("time step"),
        y_label: Some("Zone Temperature"),
        y_units: Some("C"),

        expected,
        found,

        ..validate::SeriesValidator::default()
    };

    validations.push(Box::new(v));
    validations.validate().unwrap();





    let (expected, found) = march_with_window();
    let v = validate::SeriesValidator {
        x_label: Some("time step"),
        y_label: Some("Zone Temperature"),
        y_units: Some("C"),

        expected,
        found,

        ..validate::SeriesValidator::default()
    };
    validations.push(Box::new(v));
    validations.validate().unwrap();


    let (expected, found) = march_with_window_and_luminaire();
    let v = validate::SeriesValidator {
        x_label: Some("time step"),
        y_label: Some("Zone Temperature"),
        y_units: Some("C"),

        expected,
        found,

        ..validate::SeriesValidator::default()
    };
    validations.push(Box::new(v));
    validations.validate().unwrap();



    let (expected, found) = march_with_window_and_heater();
    let v = validate::SeriesValidator {
        x_label: Some("time step"),
        y_label: Some("Zone Temperature"),
        y_units: Some("C"),

        expected,
        found,

        ..validate::SeriesValidator::default()
    };
    validations.push(Box::new(v));
    validations.validate().unwrap();


    let (expected, found) = march_with_window_heater_and_infiltration();
    let v = validate::SeriesValidator {
        x_label: Some("time step"),
        y_label: Some("Zone Temperature"),
        y_units: Some("C"),

        expected,
        found,

        ..validate::SeriesValidator::default()
    };
    validations.push(Box::new(v));
    validations.validate().unwrap();


    
    let (expected, found) = march_solar_radiation_massive();
    let v = validate::SeriesValidator {
        x_label: Some("time step"),
        y_label: Some("Zone Temperature"),
        y_units: Some("C"),
        found_name: "Simple",
        expected_name: "EnergyPlus",

        expected,
        found,

        ..validate::SeriesValidator::default()
    };
    validations.push(Box::new(v));
    validations.validate().unwrap();


    let (expected, found) = march_solar_radiation_mixed_mass();
    let v = validate::SeriesValidator {
        x_label: Some("time step"),
        y_label: Some("Zone Temperature"),
        y_units: Some("C"),
        found_name: "Simple",
        expected_name: "EnergyPlus",

        expected,
        found,

        ..validate::SeriesValidator::default()
    };
    validations.push(Box::new(v));
    validations.validate().unwrap();
}