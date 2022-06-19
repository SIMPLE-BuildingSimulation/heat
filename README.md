# Welcome to SIMPLE's Thermal Simulation Module

![build badge](https://github.com/SIMPLE-BuildingSimulation/thermal/actions/workflows/build.yaml/badge.svg)
![docs badge](https://github.com/SIMPLE-BuildingSimulation/thermal/actions/workflows/docs.yaml/badge.svg)
![tests badge](https://github.com/SIMPLE-BuildingSimulation/thermal/actions/workflows/tests.yaml/badge.svg)
[![codecov](https://codecov.io/gh/SIMPLE-BuildingSimulation/thermal/branch/main/graph/badge.svg?token=X6RV5WE0UL)](https://codecov.io/gh/SIMPLE-BuildingSimulation/thermal)
[![Clippy check](https://github.com/SIMPLE-BuildingSimulation/thermal/actions/workflows/style.yaml/badge.svg)](https://github.com/SIMPLE-BuildingSimulation/thermal/actions/workflows/style.yaml)


This module reads a [SIMPLE model](https://github.com/SIMPLE-BuildingSimulation/simple_model) and estimates thermodynamic properties. **It is still under development**. Check the documentation [Here](https://simple-buildingsimulation.github.io/thermal/rustdoc/doc/thermal/index.html)

# Is it accurate?

Check the automatic Validation report [HERE](https://simple-buildingsimulation.github.io/thermal/validation/walls.html)

## Some features

* Walls are modelled through Finite Difference method, and the solution is found through a [Runge-Kutta method](https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods)
* Materials are described as they are, the module will choose which ones should be considered no-mass.
* Fenestration and Walls are all treated equally. (Photons don't care whether they are reaching a window or a wall or a door; they just bounce and heat stuff up)
* The temperature of thermal zones are updated through [an analytical equation](https://simple-buildingsimulation.github.io/thermal/thermal/model/struct.ThermalModel.html#method.calculate_zones_abc)

Some To Do's (we would love your help)

* We have some issue in materials with Mass and No-Mass materials (check [HERE](https://simple-buildingsimulation.github.io/thermal/validation/walls.html))



