# Welcome to SIMPLE's Thermal Simulation Module

![build badge](https://github.com/SIMPLE-BuildingSimulation/thermal/actions/workflows/build.yaml/badge.svg)
![docs badge](https://github.com/SIMPLE-BuildingSimulation/thermal/actions/workflows/docs.yaml/badge.svg)
![tests badge](https://github.com/SIMPLE-BuildingSimulation/thermal/actions/workflows/tests.yaml/badge.svg)


This module reads a [SIMPLE model](https://github.com/SIMPLE-BuildingSimulation/simple_model) and estimates the temperatures. **It is still under development**. Check the documentation [Here](https://simple-buildingsimulation.github.io/thermal/)

Some features:

* Walls are modelled through Finite Difference method, and the solution is found through a [Runge-Kutta method](https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods)
* Materials are described as they are, the module will choose which ones should be considered no-mass
* The temperature of thermal zones are updated through [an analytical equation](https://simple-buildingsimulation.github.io/thermal/thermal/model/struct.ThermalModel.html#method.calculate_zones_abc)

Some To Do's (we would love your help)

* Calculate Solar Heat Gains. I have not developed this because I am thinking on how to calculate these things through [Ray-Tracing](https://github.com/SIMPLE-BuildingSimulation/rendering). Open for discussion
* Proper handling of fenestration (for now, we are treating them as walls, i.e., as surfaces with several layers)
* ... others



