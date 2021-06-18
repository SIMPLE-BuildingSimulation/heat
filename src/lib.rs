//! A Finite Difference-based Thermal simulation module
//!
//! It uses finite differences for marching forward in time and also
//! for calculating the heat transfer through walls.

pub mod construction;
pub mod heating_cooling;
pub mod model;
pub mod surface;
pub mod zone;
