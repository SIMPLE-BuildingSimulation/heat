var searchIndex = JSON.parse('{\
"thermal":{"doc":"A Finite Difference-based Thermal simulation module.","t":[6,17,17,0,0,0,0,0,0,0,0,0,3,12,11,11,11,11,12,12,11,11,12,12,11,12,11,11,11,11,11,13,13,3,13,13,4,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,12,11,11,12,11,11,11,11,11,12,11,11,11,12,12,3,12,12,11,11,11,11,11,12,11,11,11,11,12,12,11,11,11,11,3,11,11,11,11,11,11,11,11,11,12,11,11,11,12,5,11,11,11,12,5,5,5,5,5,5,11,11,12,11,11,11,11,11,3,11,12,11,12,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,12,11,12,11,12,11,11,11,11,5,3,11,11,11,12,11,12,11,11,12,11,11,11,11,11,11,11,12,11,11,11,12,8,6,6,3,10,10,10,10,10,10,10,10,10,12,12,12,10,12,10,12,10,11,11,11,12,10,11,12,12,10,12,10,12,10,11,11,11,10,11,12,11,12,12,5,11,10,11,10,11,11,11,11,3,11,11,11,11,11,11,12,11,11,11,12],"n":["Float","PI","SIGMA","cavity","construction","environment","gas","glazing","heating_cooling","model","surface","zone","Cavity","angle","borrow","borrow_mut","clone","clone_into","ein","eout","fmt","from","gas","height","into","thickness","to_owned","try_from","try_into","type_id","u_value","Back","Cavity","Discretization","None","Solid","UValue","borrow","borrow","borrow_mut","borrow_mut","build","chunk_segments","clone","clone_into","default","discretize_construction","fmt","from","from","get_chunks","get_k_q","into","into","n_elements","new","r_value","segments","to_owned","try_from","try_from","try_into","try_into","tstep_subdivision","type_id","type_id","u_value","0","0","Environment","air_speed","air_temperature","borrow","borrow_mut","clone","clone_into","default","env_emmisivity","fmt","from","get_hs","into","ir_irrad","solar_radiation","to_owned","try_from","try_into","type_id","Gas","air","argon","borrow","borrow_mut","cavity_convection","clone","clone_into","density","dynamic_viscosity","dynamic_viscosity","fmt","from","heat_capacity","heat_capacity","in_kelvin","into","krypton","mass","mass","nu_0_60","nu_60","nu_60_90","nu_90","nu_90_180","nusselt","raleigh","thermal_conductivity","thermal_conductivity","to_owned","try_from","try_into","type_id","xenon","Glazing","alpha_back","alpha_back","alpha_front","alpha_front","alphas","borrow","borrow_mut","clone","clone_into","combine","combine_layers","combined_alphas","combined_rho_back","combined_rho_front","combined_tau","fmt","from","get_back_glazing_system","get_front_glazing_system","get_glazing_from_iter","into","new","rho_back","rho_back","rho_front","rho_front","tau","tau","to_owned","try_from","try_into","type_id","calc_cooling_heating_power","ThermalModel","borrow","borrow_mut","calculate_zones_abc","dt","dt_subdivisions","dt_subdivisions","estimate_zones_future_temperatures","estimate_zones_mean_future_temperatures","fenestrations","from","get_current_zones_temperatures","get_thermal_zone","into","march","module_name","new","surfaces","try_from","try_into","type_id","zones","SurfaceTrait","ThermalFenestration","ThermalSurface","ThermalSurfaceData","add_back_convection_state","add_back_convective_heatflow_state","add_back_ir_irradiance_state","add_back_solar_irradiance_state","add_front_convection_state","add_front_convective_heatflow_state","add_front_ir_irradiance_state","add_front_solar_irradiance_state","add_node_temperature_states","area","back_alphas","back_boundary","back_convection_coefficient","back_emmisivity","back_infrared_irradiance","back_solar_absorbtance","back_solar_irradiance","back_temperature","borrow","borrow_mut","discretization","first_node_temperature_index","from","front_alphas","front_boundary","front_convection_coefficient","front_emmisivity","front_infrared_irradiance","front_solar_absorbtance","front_solar_irradiance","front_temperature","get_node_temperatures","into","last_node_temperature_index","march","massive_chunks","new","nomass_chunks","parent","rk4","set_back_boundary","set_back_convection_coefficient","set_front_boundary","set_front_convection_coefficient","set_node_temperatures","try_from","try_into","type_id","ThermalZone","borrow","borrow_mut","from","from_space","into","mcp","reference_space","try_from","try_into","type_id","volume"],"q":["thermal","","","","","","","","","","","","thermal::cavity","","","","","","","","","","","","","","","","","","","thermal::construction","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","thermal::construction::UValue","","thermal::environment","","","","","","","","","","","","","","","","","","","thermal::gas","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","thermal::glazing","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","thermal::heating_cooling","thermal::model","","","","","","","","","","","","","","","","","","","","","","thermal::surface","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","thermal::zone","","","","","","","","","","",""],"d":["The kind of Floating point number used in the library… …","Well, Pi","The Stefan–Boltzmann constant (in $<code>W m^{-2} K^4</code>$), …","","","","","","","","","","Represents some gas enclosed by two solid materials","The angle of the cavity in radians. $<code>0</code>$ is horizontal; $…","","","","","The thermal emissivity of the material at the inner side …","The thermal emissivity of the material at the outer side …","","Returns the argument unchanged.","The gas contained","Height of the <code>Cavity</code>. Defined by ISO15099/2003 as “the …","Calls <code>U::from(self)</code>.","The distance between the two materials, in $<code>m</code>$","","","","","Calculates the <code>U-value</code>—including convective and …","The resistance is a surface coefficient.","A cavity, comprised of a gas","Represents the discretization of a <code>Construction</code> for heat …","Undefined yet","A normal (i.e., $<code>\\\\lambda/\\\\Delta x</code>$) U-value","Represents a thermal connection in the thermal network. It …","","","","","Creates the <code>segments</code> of the <code>Discretization</code>.","Auxiliary function for <code>get_chunks()</code> function","","","","Given a Maximum element thickness ($<code>\\\\Delta x_{max}</code>$) and a …","","Returns the argument unchanged.","Returns the argument unchanged.","Gets a a the segments that correspond to Massive and …","Produces $<code>\\\\overline{K}</code>$ and $<code>\\\\vec{q}</code>$ (as in the equation $…","Calls <code>U::from(self)</code>.","Calls <code>U::from(self)</code>.","The number of elements on each layer","Creates a new <code>Discretization</code>.","Calculates the R value of the whole system","Contains the node’s mass and the <code>UValue</code> of each segment","","","","","","Contains the minimum number of timesteps per model …","","","Gets the U-value of a <code>UValue</code> object","","","Represents a border condition of between a Surface and a …","The wind speed, in m/2","The dry bulb temperature of the air, in $<code>C</code>$","","","","","","The environmental emmisivity, used for calculating …","","Returns the argument unchanged.","","Calls <code>U::from(self)</code>.","The incident Infrared Irradiance, in $<code>W/m^2</code>$","The incident Solar Irradiance, in $<code>W/m^2</code>$","","","","","A structure containing the data that will describe the …","Returns a gas with the properties of Air","Returns a gas with the properties of argon","","","Calculates the convective heat transfer coefficient within …","","","Derives the density based on the temperature (in $<code>K</code>$)","Derives the Dynamic Viscosity at a certain Temperature (in …","The dynamic viscosity ( $<code>{N.s}/{m^2}</code>$) as a function of the","","Returns the argument unchanged.","Derives the Soecific Heat Capacity at a certain …","The specific heat capacity ($<code>{J}/{kg.K}</code>$) as a function of …","Transforms C into K","Calls <code>U::from(self)</code>.","Returns a gas with the properties of krypton","Retreives the Molecular Mass","THe Molecular Mass ($<code>{kg}/{Mol}</code>$)","Calculates the Nusselt number for cavities tilted between $…","Calculates the Nusselt number for cavities tilted $<code>60^o</code>$","Calculates the Nusselt number for cavities tilted between $…","Calculates the Nusselt number for cavities tilted $<code>90^o</code>$","Calculates the Nusselt number for cavities tilted between $…","Calculates the Nusselt of a cavity number based on the …","Calculates the Raleigh number of a <code>Gas</code> cavity based on its …","Derives the Thermal Conductivity at a certain Temperature …","The thermal conductivity ($<code>{W}/{m.K}</code>$) as a function of the","","","","","Returns a gas with the properties of xenon","An abstraction of a glazing layer for optical purposes.","Gets the back absorbtance","Absorbtance $<code>\\\\alpha_b</code>$ on the back side","Gets the front absorbtance","Absorbtance $<code>\\\\alpha_f</code>$ on the front side","Calculates the absorbtances of each <code>Glazing</code> of the system, …","","","","","Combines two <code>Glazing</code> into a new <code>Glazing</code>","Combines several <code>Glazing</code> into a new <code>Glazing</code>","Calculates the front solar absorbtance of two <code>Glazing</code> …","Calculates the overall back reflectance of a system of two …","Calculates the overall front reflectance of a system of …","Calculates the overall transmittance of a system of two …","","Returns the argument unchanged.","","","","Calls <code>U::from(self)</code>.","Creates a new <code>Glazing</code>","Gets the back reflectance","Reflectance $<code>\\\\rho_b</code>$ on the back side","Gets the front reflectance","Reflectance $<code>\\\\rho_f</code>$ on the front side.","Gets the transmittance","Transmittance $<code>\\\\tau</code>$","","","","","Retrieves a <code>Vec&lt;(usize, Float)&gt;</code> containing the amount of …","","","","This estimation assumes nothing changes during this time. …","The model’s dt (i.e., main_dt / self.dt_subdivisions)","Retrieves the dt_subdivisions (i.e. the number of …","The number of steps that this model needs to take in order …","Uses an analytical solution to estimate the future Zones …","Uses an analytical solution to estimate an average …","All the Fenestrations in the model","Returns the argument unchanged.","Retrieves a vector of the current temperatures of all the …","Retrieves a ThermalZone","Calls <code>U::from(self)</code>.","Advances one main_timestep through time. That is, it …","","Creates a new ThermalModel from a SimpleModel.","All the surfaces in the model","","","","All the Thermal Zones in the model","","","","This is a Surface from the point of view of our thermal …","","","","","","","","","","The area of the Surface","The absorbtances of each node in the system, proportional …","The location of the back boundary zone in the Zones array …","","The thermal absorbtance on the back side (from 0 to 1)","","The solar absorbtance on the back side (from 0 to 1)","","","","","The <code>Discretization</code> that represents this <code>ThermalSurfaceData</code>","","Returns the argument unchanged.","The absorbtances of each node in the system, proportional …","The location of the front boundary zone in the Zones array …","","The thermal absorbtance on the front side (from 0 to 1)","","The solar absorbtance on the front side (from 0 to 1)","","","","Calls <code>U::from(self)</code>.","","Marches one timestep. Returns front and back heat flow    ","The chunks of nodes that have mass","","The chunks of nodes that have nomass","","Marches forward through time, solving the Ordinary …","","","","","","","","","","","","Returns the argument unchanged.","This function creates a new ThermalZone from a Space. It …","Calls <code>U::from(self)</code>.","Retrieves the heat capacity of the ThermalZone’s air","The <code>Space</code> that this [<code>Thermal Zone</code>] represents","","","","volume of the zone"],"i":[0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,12,12,0,12,12,0,10,12,10,12,10,10,12,12,12,10,12,10,12,10,10,10,12,10,10,10,10,12,10,12,10,12,10,10,12,12,32,33,0,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,0,15,15,15,15,15,15,15,15,15,15,15,15,15,15,0,15,15,15,15,0,0,0,0,0,0,15,15,15,15,15,15,15,15,0,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,0,0,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,0,0,0,0,29,29,29,29,29,29,29,29,29,30,30,30,29,30,29,30,29,29,30,30,30,29,30,30,30,29,30,29,30,29,29,29,30,29,30,30,30,30,30,0,30,29,30,29,29,30,30,30,0,22,22,22,22,22,22,22,22,22,22,22],"f":[0,0,0,0,0,0,0,0,0,0,0,0,0,0,[[]],[[]],[1,1],[[]],0,0,[[1,2],3],[[]],0,0,[[]],0,[[]],[[],4],[[],4],[[],5],[[1,6,6],6],0,0,0,0,0,0,[[]],[[]],[[]],[[]],[[7,8,[9,[8]],6,6],[[4,[10,11]]]],[10,9],[12,12],[[]],[[],12],[[7,6,6,6]],[[12,2],3],[[]],[[]],[10],[[10,8,8,13,14,6,6,14,6,6]],[[]],[[]],0,[[7,6,6,6,6,6],[[4,[10,11]]]],[10,6],0,[[]],[[],4],[[],4],[[],4],[[],4],0,[[],5],[[],5],[[12,6,6],6],0,0,0,0,0,[[]],[[]],[14,14],[[]],[[],14],0,[[14,2],3],[[]],[14,6],[[]],0,0,[[]],[[],4],[[],4],[[],5],0,[[],15],[[],15],[[]],[[]],[[15,6,6,6,6,6],6],[15,15],[[]],[[15,6],6],[[15,6],6],0,[[15,2],3],[[]],[[15,6],6],0,[6,6],[[]],[[],15],[15,6],0,[[6,6,6],6],[[6,6],6],[[6,6,6],6],[[6,6],6],[[6,6,6],6],[[6,6,6],6],[[15,6,6,6],6],[[15,6],6],0,[[]],[[],4],[[],4],[[],5],[[],15],0,[16,6],0,[16,6],0,[[],[[9,[6]]]],[[]],[[]],[16,16],[[]],[[16,16],16],[[],16],[[16,16]],[[16,16],6],[[16,16],6],[[16,16],6],[[16,2],3],[[]],[17,[[4,[[9,[16]],11]]]],[17,[[4,[[9,[16]],11]]]],[8,[[4,[[9,[16]],11]]]],[[]],[[6,6,6],16],[16,6],0,[16,6],0,[16,6],0,[[]],[[],4],[[],4],[[],5],[[18,19],9],0,[[]],[[]],[[20,21,19]],0,[20,8],0,[[20,6],[[9,[6]]]],[[20,6],[[9,[6]]]],0,[[]],[[20,19],[[9,[6]]]],[[20,8],[[4,[22,11]]]],[[]],[[20,23,24,21,19],[[4,[11]]]],[[],25],[[26,21,27,8],[[4,[20,11]]]],0,[[],4],[[],4],[[],5],0,0,0,0,0,[[27,8]],[[27,8]],[[27,8]],[[27,8]],[[27,8]],[[27,8]],[[27,8]],[[27,8]],[[27,8,8]],0,0,0,[19,[[28,[6]]]],0,[19,6],0,[19,6],[19,6],[[]],[[]],0,[[],8],[[]],0,0,[19,[[28,[6]]]],0,[19,6],0,[19,6],[19,6],[19,13],[[]],[[],8],[[[30,[29]],19,6,6,6]],0,[[27,8,7,6,7,10],[[4,[[30,[29]],11]]]],0,0,[[6,13,13,13,13]],[[[30,[29]],31]],[[19,6]],[[[30,[29]],31]],[[19,6]],[[19,13]],[[],4],[[],4],[[],5],0,[[]],[[]],[[]],[[7,27,8],22],[[]],[[22,6],6],0,[[],4],[[],4],[[],5],0],"p":[[3,"Cavity"],[3,"Formatter"],[6,"Result"],[4,"Result"],[3,"TypeId"],[6,"Float"],[3,"Rc"],[15,"usize"],[3,"Vec"],[3,"Discretization"],[3,"String"],[4,"UValue"],[6,"Matrix"],[3,"Environment"],[3,"Gas"],[3,"Glazing"],[3,"Construction"],[4,"HVAC"],[6,"SimulationState"],[3,"ThermalModel"],[3,"SimpleModel"],[3,"ThermalZone"],[3,"Date"],[8,"Weather"],[15,"str"],[3,"MetaOptions"],[3,"SimulationStateHeader"],[4,"Option"],[8,"SurfaceTrait"],[3,"ThermalSurfaceData"],[4,"Boundary"],[13,"Solid"],[13,"Cavity"]]}\
}');
if (typeof window !== 'undefined' && window.initSearch) {window.initSearch(searchIndex)};
if (typeof exports !== 'undefined') {exports.searchIndex = searchIndex};