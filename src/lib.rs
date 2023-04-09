#![doc = include_str!("../README.md")]

pub mod common;
mod hashgrid;
pub mod utils;

pub use hashgrid::{HashGrid, HashCell, HashGridError, PeriodicImage};