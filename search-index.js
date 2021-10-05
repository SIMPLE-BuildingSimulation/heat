var searchIndex = JSON.parse('{\
"thermal":{"doc":"A Finite Difference-based Thermal simulation module","t":[6,0,0,0,0,0,5,5,5,5,5,5,5,5,5,3,11,11,11,12,11,12,11,11,12,11,11,11,11,11,11,11,11,11,12,11,11,11,12,13,13,4,3,11,12,11,12,11,11,11,11,11,11,11,11,11,11,12,11,12,12,11,11,11,11,12,11,12,11,12,11,11,11,11,12,11,12,11,11,11,12,11,11,11,11,11,11,3,11,11,11,11,11,11,12,11,11,11,12],"n":["Float","construction","heating_cooling","model","surface","zone","build_thermal_network","calc_c_matrix","calc_full_rs_back","calc_full_rs_front","calc_k_matrix","calc_n_total_nodes","discretize_construction","get_first_and_last_massive_elements","calc_cooling_heating_power","ThermalModel","borrow","borrow_mut","calculate_zones_abc","dt","dt_subdivisions","dt_subdivisions","estimate_zones_future_temperatures","estimate_zones_mean_future_temperatures","fenestrations","from","get_current_zones_temperatures","get_thermal_fenestration","get_thermal_surface","get_thermal_zone","into","march","module_name","new","surfaces","try_from","try_into","type_id","zones","Fenestration","Surface","ThermalSurface","ThermalSurfaceData","area","area","back_boundary","back_boundary","back_temperature","borrow","borrow","borrow_mut","borrow_mut","calc_heat_flow","data","from","from","front_boundary","front_boundary","front_temperature","full_rs_back","full_rs_front","get_node_temperatures","into","into","is_massive","kt4_func","march","massive","mut_data","n_nodes","new","new_fenestration","new_surface","rs_back","rs_back","rs_front","rs_front","set_back_boundary","set_front_boundary","set_node_temperatures","total_r","try_from","try_from","try_into","try_into","type_id","type_id","ThermalZone","borrow","borrow_mut","from","from_space","into","mcp","reference_space","try_from","try_into","type_id","volume"],"q":["thermal","","","","","","thermal::construction","","","","","","","","thermal::heating_cooling","thermal::model","","","","","","","","","","","","","","","","","","","","","","","","thermal::surface","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","thermal::zone","","","","","","","","","","",""],"d":["","","","","","","Builds the necessary data for marching forward through …","","","","Calculates the <code>K</code> matrix (i.e., the thermal network) for …","This function calculates the number of nodes that result …","Given a Maximum element thickness ($<code>\\\\Delta x_{max}</code>$) and …","In a discretization scheme, this function finds the first …","Retrieves a <code>Vec<(usize, Float)></code> containing the amount of …","","","","This estimation assumes nothing changes during this time. …","The model’s dt (i.e., main_dt / self.dt_subdivisions)","Retrieves the dt_subdivisions (i.e. the number of …","The number of steps that this model needs to take in …","Uses an analytical solution to estimate the future Zones …","Uses an analytical solution to estimate an average …","All the Fenestrations in the model","","Retrieves a vector of the current temperatures of all the …","Retrieves a THermalFenestration","Retrieves a ThermalSurface","Retrieves a ThermalZone","","Advances one main_timestep through time. That is, it …","","Creates a new ThermalModel from a SimpleModel.","All the surfaces in the model","","","","All the Thermal Zones in the model","","","","This is a Surface from the point of view of our thermal …","Gets the <code>area</code> of the surface","The area of the Surface","","The location of the back boundary zone in the Zones array …","","","","","","Calculates the heat flow out of the layer, based on the …","Borrows the <code>ThermalSurfaceData</code> of this surface","","","","The location of the front boundary zone in the Zones …","","The exterior (i.e. back side) resistance after the last …","The interior (i.e. front side) resistance before any …","Retrieves the state of the Surface as a Matrix object.","","","Checks whether a wall has thermal mass","","Marches one timestep. Returns front and back heat flow    ","Has thermal mass at all?","Borrows a mutable version of <code>ThermalSurfaceData</code> of this …","The number of nodes after discretizing the construction","Constructs a new ThermalSurface object.        ","","","Gets the <code>rs_back</code>","The back side convection coefficient","Gets the <code>rs_front</code>","The front side convection coefficient","","","","","","","","","","","","","","","This function creates a new ThermalZone from a Space. It …","","Retrieves the heat capacity of the ThermalZone’s air","The <code>Space</code> that this [<code>Thermal Zone</code>] represents","","","","volume of the zone"],"i":[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,0,0,2,3,2,3,2,2,3,2,3,2,2,2,3,2,3,2,3,3,2,2,3,2,3,2,3,2,3,3,2,2,2,3,2,3,2,2,2,3,2,3,2,3,2,3,0,4,4,4,4,4,4,4,4,4,4,4],"f":[null,null,null,null,null,null,[[["usize",15],["rc",3],["f64",15]],[["string",3],["result",4,["string"]]]],[[["usize",15],["rc",3]],[["vec",3,["f64"]],["f64",15]]],[[["usize",15],["rc",3],["f64",15]],["f64",15]],[[["usize",15],["rc",3],["f64",15]],["f64",15]],[[["usize",15],["rc",3],["f64",15]],["matrix",3]],[[],[["usize",15],["string",3],["result",4,["usize","string"]]]],[[["rc",3],["f64",15]]],[[],[["result",4,["string"]],["string",3]]],[[["rc",3],["simulationstate",6]],["vec",3]],null,[[]],[[]],[[["simplemodel",3],["simulationstate",6]]],null,[[],["usize",15]],null,[[["f64",15]],[["vec",3,["f64"]],["f64",15]]],[[["f64",15]],[["vec",3,["f64"]],["f64",15]]],null,[[]],[[["simulationstate",6]],[["vec",3,["f64"]],["f64",15]]],[[["usize",15]],[["string",3],["result",4,["thermalsurface","string"]],["thermalsurface",4]]],[[["usize",15]],[["string",3],["result",4,["thermalsurface","string"]],["thermalsurface",4]]],[[["usize",15]],[["result",4,["thermalzone","string"]],["thermalzone",3],["string",3]]],[[]],[[["simulationstate",6],["date",3],["weather",8],["simplemodel",3]],[["result",4,["string"]],["string",3]]],[[],["str",15]],[[["usize",15],["simplemodel",3],["simulationstateheader",3]],[["result",4,["string"]],["string",3]]],null,[[],["result",4]],[[],["result",4]],[[],["typeid",3]],null,null,null,null,null,[[],["f64",15]],null,[[],["option",4]],null,[[["simulationstate",6]],["f64",15]],[[]],[[]],[[]],[[]],[[["simulationstate",6],["f64",15]]],[[],["thermalsurfacedata",3]],[[]],[[]],[[],["option",4]],null,[[["simulationstate",6]],["f64",15]],null,null,[[["simulationstate",6]],["matrix",3]],[[]],[[]],[[],["bool",15]],null,[[["f64",15],["simulationstate",6]]],null,[[],["thermalsurfacedata",3]],null,[[["usize",15],["simulationstateheader",3],["bool",15],["rc",3],["f64",15]],[["result",4,["string"]],["string",3]]],[[["rc",3],["simulationstateheader",3],["f64",15]],[["result",4,["string"]],["string",3]]],[[["simulationstateheader",3],["rc",3],["f64",15]],[["result",4,["string"]],["string",3]]],[[],["f64",15]],null,[[],["f64",15]],null,[[["boundary",4]]],[[["boundary",4]]],[[["matrix",3],["simulationstate",6]]],null,[[],["result",4]],[[],["result",4]],[[],["result",4]],[[],["result",4]],[[],["typeid",3]],[[],["typeid",3]],null,[[]],[[]],[[]],[[["usize",15],["rc",3],["simulationstateheader",3]]],[[]],[[],["f64",15]],null,[[],["result",4]],[[],["result",4]],[[],["typeid",3]],null],"p":[[3,"ThermalModel"],[4,"ThermalSurface"],[3,"ThermalSurfaceData"],[3,"ThermalZone"]]}\
}');
if (window.initSearch) {window.initSearch(searchIndex)};