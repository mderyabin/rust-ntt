pub mod ntmath;
pub mod polymath;
pub use ntmath::*;
pub use polymath::*;

#[cfg(test)]
mod ntmath_tests;

#[cfg(test)]
mod polymath_tests;