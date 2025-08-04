use rand::Rng;
use std::ops::{AddAssign, MulAssign, Neg};

// Core polynomial ring trait with miniumum set of operations for implementing CKKS.
pub trait PolyRing<const DEGREE: usize>:
    Clone + for<'a> AddAssign<&'a Self> + for<'a> MulAssign<&'a Self> + Neg<Output = Self>
{
    type Context;

    fn zero(context: &Self::Context) -> Self;
    fn from_coeffs(coeffs: &[i64], context: &Self::Context) -> Self;
    fn to_coeffs(&self) -> [i64; DEGREE];

    fn context(&self) -> &Self::Context;
}

// Sampling trait - provides common sampling operations for polynomials
pub trait PolySampler<const DEGREE: usize>: PolyRing<DEGREE> {
    fn sample_uniform<R: Rng>(rng: &mut R, context: &Self::Context) -> Self;
    fn sample_gaussian<R: Rng>(rng: &mut R, std_dev: f64, context: &Self::Context) -> Self;
    fn sample_tribits<R: Rng>(rng: &mut R, hamming_weight: usize, context: &Self::Context) -> Self;
    fn sample_noise<R: Rng>(rng: &mut R, variance: f64, context: &Self::Context) -> Self;
}
