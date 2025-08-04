use crate::context::NttContext;
use std::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};
use std::sync::Arc;

/// Polynomial in NTT-friendly form with shared context
#[derive(Debug, Clone)]
pub struct NttPolynomial<const DEGREE: usize> {
    coeffs: [u64; DEGREE],
    context: Arc<NttContext<DEGREE>>,
}

impl<const DEGREE: usize> NttPolynomial<DEGREE> {
    /// Create polynomial from coefficients
    pub fn from_coeffs(
        coeffs: [u64; DEGREE],
        context: Arc<NttContext<DEGREE>>,
    ) -> Self {
        Self { coeffs, context }
    }

    /// Create zero polynomial
    pub fn zero(context: Arc<NttContext<DEGREE>>) -> Self {
        Self {
            coeffs: [0u64; DEGREE],
            context,
        }
    }

    /// Get coefficients
    pub fn coeffs(&self) -> &[u64; DEGREE] {
        &self.coeffs
    }

    /// Get mutable coefficients
    pub fn coeffs_mut(&mut self) -> &mut [u64; DEGREE] {
        &mut self.coeffs
    }

    /// Get context
    pub fn context(&self) -> &Arc<NttContext<DEGREE>> {
        &self.context
    }

    // NTT operations
    pub fn ntt_forward(&mut self) {
        // Cooley-Tukey forward negacyclic NTT
        // using algorithm from https://eprint.iacr.org/2016/504.pdf

        let mut t = DEGREE >> 1;
        let mut n = 1;

        while n < DEGREE {
            for i in 0..n {
                let j1 = 2 * i * t;
                let j2 = j1 + t - 1;
                let s = self.context.tf[n + i];

                for j in j1..=j2 {
                    let u = self.coeffs[j];
                    let v = self.context.class.modmul(self.coeffs[j + t], s);

                    self.coeffs[j] = self.context.class.modadd(u, v);
                    self.coeffs[j + t] = self.context.class.modsub(u, v);
                }
            }

            n <<= 1;
            t >>= 1;
        }
    }

    pub fn ntt_inverse(&mut self) {
        // Gentleman-Sande inverse negacyclic NTT
        let mut t = 1;
        let mut h = DEGREE >> 1;

        while h > 0 {
            let mut j1 = 0;

            for i in 0..h {
                let j2 = j1 + t - 1;
                let s = self.context.itf[h + i];

                for j in j1..=j2 {
                    let u = self.coeffs[j];
                    let v = self.coeffs[j + t];

                    self.coeffs[j] = self.context.class.modadd(u, v);
                    self.coeffs[j + t] = self.context.class.modsub(u, v);
                    self.coeffs[j + t] =
                        self.context.class.modmul(self.coeffs[j + t], s);
                }

                j1 += t << 1;
            }

            h >>= 1;
            t <<= 1;
        }

        // Final normalization
        for coeff in &mut self.coeffs {
            *coeff = self.context.class.modmul(*coeff, self.context.inv_n);
        }
    }

    pub fn ntt_forward_shoup(&mut self) {
        // Cooley-Tukey forward negacyclic NTT with Shoup multiplication
        let mut t = DEGREE >> 1;
        let mut n = 1;

        while n < DEGREE {
            for i in 0..n {
                let j1 = 2 * i * t;
                let j2 = j1 + t - 1;
                let s = self.context.tf[n + i];
                let s_shoup = self.context.tf_shoup[n + i];

                for j in j1..=j2 {
                    let v = self.context.class.modmul_shoup(
                        self.coeffs[j + t],
                        s,
                        s_shoup,
                    );

                    self.coeffs[j + t] =
                        self.context.class.modsub(self.coeffs[j], v);
                    self.context.class.modadd_eq(&mut self.coeffs[j], v);
                }
            }

            n <<= 1;
            t >>= 1;
        }
    }

    pub fn ntt_inverse_shoup(&mut self) {
        // Gentleman-Sande inverse negacyclic NTT with Shoup multiplication
        let mut t = 1;
        let mut h = DEGREE >> 1;

        while h > 0 {
            let mut j1 = 0;

            for i in 0..h {
                let j2 = j1 + t - 1;
                let s = self.context.itf[h + i];
                let s_shoup = self.context.itf_shoup[h + i];

                for j in j1..=j2 {
                    let u = self.coeffs[j];
                    let v = self.coeffs[j + t];

                    self.coeffs[j] = self.context.class.modadd(u, v);
                    self.coeffs[j + t] = self.context.class.modsub(u, v);
                    self.context.class.modmul_shoup_eq(
                        &mut self.coeffs[j + t],
                        s,
                        s_shoup,
                    );
                }

                j1 += t << 1;
            }

            h >>= 1;
            t <<= 1;
        }

        // Final normalization with Shoup
        for coeff in &mut self.coeffs {
            self.context.class.modmul_shoup_eq(
                coeff,
                self.context.inv_n,
                self.context.inv_n_shoup,
            );
        }
    }

    // Convolution methods
    pub fn negacyclic_convolution(&self, other: &Self) -> Self {
        debug_assert_eq!(
            self.context.modulus(),
            other.context.modulus(),
            "Cannot convolve polynomials with different moduli"
        );

        let mut result = self.clone();
        let mut other_copy = other.clone();

        result.ntt_forward();
        other_copy.ntt_forward();

        // Pointwise multiplication in NTT domain
        for i in 0..DEGREE {
            result.coeffs[i] = self
                .context
                .class
                .modmul(result.coeffs[i], other_copy.coeffs[i]);
        }

        result.ntt_inverse();
        result
    }

    pub fn negacyclic_convolution_shoup(&self, other: &Self) -> Self {
        debug_assert_eq!(
            self.context.modulus(),
            other.context.modulus(),
            "Cannot convolve polynomials with different moduli"
        );

        let mut result = self.clone();
        let mut other_copy = other.clone();

        result.ntt_forward_shoup();
        other_copy.ntt_forward_shoup();

        // Pointwise multiplication in NTT domain
        for i in 0..DEGREE {
            result.coeffs[i] = self
                .context
                .class
                .modmul(result.coeffs[i], other_copy.coeffs[i]);
        }

        result.ntt_inverse_shoup();
        result
    }

    // Sampling utility
    pub fn sample_random(context: Arc<NttContext<DEGREE>>) -> Self {
        use rand::{Rng, rng};

        let mut generator = rng();
        let mut coeffs = [0u64; DEGREE];

        for coeff in &mut coeffs {
            *coeff = generator.random_range(1..context.modulus());
        }

        Self { coeffs, context }
    }
}

// Trait implementations - this is where the math logic lives

impl<const DEGREE: usize> Add for &NttPolynomial<DEGREE> {
    type Output = NttPolynomial<DEGREE>;

    fn add(self, rhs: Self) -> Self::Output {
        debug_assert_eq!(
            self.context.modulus(),
            rhs.context.modulus(),
            "Cannot add polynomials with different moduli"
        );

        let mut result_coeffs = [0u64; DEGREE];
        for i in 0..DEGREE {
            result_coeffs[i] =
                self.context.class.modadd(self.coeffs[i], rhs.coeffs[i]);
        }

        NttPolynomial {
            coeffs: result_coeffs,
            context: Arc::clone(&self.context),
        }
    }
}

impl<const DEGREE: usize> Add<&NttPolynomial<DEGREE>> for NttPolynomial<DEGREE> {
    type Output = NttPolynomial<DEGREE>;

    fn add(self, rhs: &NttPolynomial<DEGREE>) -> Self::Output {
        &self + rhs
    }
}

impl<const DEGREE: usize> AddAssign<&NttPolynomial<DEGREE>>
    for NttPolynomial<DEGREE>
{
    fn add_assign(&mut self, rhs: &NttPolynomial<DEGREE>) {
        debug_assert_eq!(
            self.context.modulus(),
            rhs.context.modulus(),
            "Cannot add polynomials with different moduli"
        );

        for i in 0..DEGREE {
            self.context
                .class
                .modadd_eq(&mut self.coeffs[i], rhs.coeffs[i]);
        }
    }
}

impl<const DEGREE: usize> Sub for &NttPolynomial<DEGREE> {
    type Output = NttPolynomial<DEGREE>;

    fn sub(self, rhs: Self) -> Self::Output {
        debug_assert_eq!(
            self.context.modulus(),
            rhs.context.modulus(),
            "Cannot subtract polynomials with different moduli"
        );

        let mut result_coeffs = [0u64; DEGREE];
        for i in 0..DEGREE {
            result_coeffs[i] =
                self.context.class.modsub(self.coeffs[i], rhs.coeffs[i]);
        }

        NttPolynomial {
            coeffs: result_coeffs,
            context: Arc::clone(&self.context),
        }
    }
}

impl<const DEGREE: usize> Sub<&NttPolynomial<DEGREE>> for NttPolynomial<DEGREE> {
    type Output = NttPolynomial<DEGREE>;

    fn sub(self, rhs: &NttPolynomial<DEGREE>) -> Self::Output {
        &self - rhs
    }
}

impl<const DEGREE: usize> SubAssign<&NttPolynomial<DEGREE>>
    for NttPolynomial<DEGREE>
{
    fn sub_assign(&mut self, rhs: &NttPolynomial<DEGREE>) {
        debug_assert_eq!(
            self.context.modulus(),
            rhs.context.modulus(),
            "Cannot subtract polynomials with different moduli"
        );

        for i in 0..DEGREE {
            self.context
                .class
                .modsub_eq(&mut self.coeffs[i], rhs.coeffs[i]);
        }
    }
}

impl<const DEGREE: usize> Mul for &NttPolynomial<DEGREE> {
    type Output = NttPolynomial<DEGREE>;

    fn mul(self, rhs: Self) -> Self::Output {
        self.negacyclic_convolution(rhs)
    }
}

impl<const DEGREE: usize> Mul<&NttPolynomial<DEGREE>> for NttPolynomial<DEGREE> {
    type Output = NttPolynomial<DEGREE>;

    fn mul(self, rhs: &NttPolynomial<DEGREE>) -> Self::Output {
        &self * rhs
    }
}

impl<const DEGREE: usize> MulAssign<&NttPolynomial<DEGREE>>
    for NttPolynomial<DEGREE>
{
    fn mul_assign(&mut self, rhs: &NttPolynomial<DEGREE>) {
        *self = &*self * rhs;
    }
}

impl<const DEGREE: usize> Neg for NttPolynomial<DEGREE> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        let mut result_coeffs = [0u64; DEGREE];
        for i in 0..DEGREE {
            result_coeffs[i] = self.context.class.modneg(self.coeffs[i]);
        }

        NttPolynomial {
            coeffs: result_coeffs,
            context: self.context,
        }
    }
}

impl<const DEGREE: usize> Neg for &NttPolynomial<DEGREE> {
    type Output = NttPolynomial<DEGREE>;

    fn neg(self) -> Self::Output {
        let mut result_coeffs = [0u64; DEGREE];
        for i in 0..DEGREE {
            result_coeffs[i] = self.context.class.modneg(self.coeffs[i]);
        }

        NttPolynomial {
            coeffs: result_coeffs,
            context: Arc::clone(&self.context),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::context::NttContext;
    use crate::ntmath::find_first_prime_up;

    #[test]
    fn test_polynomial_creation() {
        const N: usize = 4;
        let q = find_first_prime_up(10, N);
        let ctx = NttContext::<N>::new(q);

        let zero = NttPolynomial::zero(Arc::clone(&ctx));
        assert_eq!(zero.coeffs(), &[0u64; N]);

        let coeffs = [1, 2, 3, 4];
        let poly = NttPolynomial::from_coeffs(coeffs, ctx);
        assert_eq!(poly.coeffs(), &coeffs);
    }

    #[test]
    fn test_basic_arithmetic() {
        const N: usize = 4;
        let q = find_first_prime_up(10, N);
        let ctx = NttContext::<N>::new(q);

        let a = NttPolynomial::from_coeffs([1, 2, 3, 4], Arc::clone(&ctx));
        let b = NttPolynomial::from_coeffs([2, 3, 4, 5], Arc::clone(&ctx));

        // Test addition
        let c = &a + &b;
        let expected_add = [3, 5, 7, 9];
        assert_eq!(c.coeffs(), &expected_add);

        // Test subtraction
        let d = &b - &a;
        let expected_sub = [1, 1, 1, 1];
        assert_eq!(d.coeffs(), &expected_sub);

        // Test negation
        let neg_a = -&a;
        for i in 0..N {
            let expected = if a.coeffs()[i] == 0 {
                0
            } else {
                q - a.coeffs()[i]
            };
            assert_eq!(neg_a.coeffs()[i], expected);
        }
    }

    #[test]
    fn test_assign_operations() {
        const N: usize = 4;
        let q = find_first_prime_up(10, N);
        let ctx = NttContext::<N>::new(q);

        let mut a = NttPolynomial::from_coeffs([1, 2, 3, 4], Arc::clone(&ctx));
        let b = NttPolynomial::from_coeffs([2, 3, 4, 5], Arc::clone(&ctx));

        // Test AddAssign
        let original_a = a.clone();
        a += &b;
        let expected = [3, 5, 7, 9];
        assert_eq!(a.coeffs(), &expected);

        // Test SubAssign
        a -= &b;
        assert_eq!(a.coeffs(), original_a.coeffs());

        // Test MulAssign with simple case
        let mut simple = NttPolynomial::from_coeffs([1, 0, 0, 0], Arc::clone(&ctx));
        let identity = simple.clone();
        simple *= &identity;
        assert_eq!(simple.coeffs(), &[1, 0, 0, 0]); // 1 * 1 = 1
    }

    #[test]
    fn test_ntt_forward_inverse() {
        const N: usize = 4;
        let q = find_first_prime_up(10, N);
        let ctx = NttContext::<N>::new(q);

        let original = NttPolynomial::sample_random(Arc::clone(&ctx));
        let mut test_poly = original.clone();

        // Forward then inverse should give back original
        test_poly.ntt_forward();
        test_poly.ntt_inverse();

        assert_eq!(test_poly.coeffs(), original.coeffs());
    }

    #[test]
    fn test_ntt_shoup_forward_inverse() {
        const N: usize = 4;
        let q = find_first_prime_up(10, N);
        let ctx = NttContext::<N>::new(q);

        let original = NttPolynomial::sample_random(Arc::clone(&ctx));
        let mut test_poly = original.clone();

        // Forward then inverse should give back original (Shoup version)
        test_poly.ntt_forward_shoup();
        test_poly.ntt_inverse_shoup();

        assert_eq!(test_poly.coeffs(), original.coeffs());
    }

    #[test]
    fn test_convolution_consistency() {
        const N: usize = 4;
        let q = find_first_prime_up(10, N);
        let ctx = NttContext::<N>::new(q);

        let a = NttPolynomial::sample_random(Arc::clone(&ctx));
        let b = NttPolynomial::sample_random(Arc::clone(&ctx));

        // Both convolution methods should give same result
        let result1 = a.negacyclic_convolution(&b);
        let result2 = a.negacyclic_convolution_shoup(&b);

        assert_eq!(result1.coeffs(), result2.coeffs());
    }

    #[test]
    fn test_multiplication_via_convolution() {
        const N: usize = 4;
        let q = find_first_prime_up(10, N);
        let ctx = NttContext::<N>::new(q);

        // Test simple multiplication: (1 + x) * (1 + x) = 1 + 2x + x^2
        let poly1 = NttPolynomial::from_coeffs([1, 1, 0, 0], Arc::clone(&ctx));
        let poly2 = poly1.clone();

        let result = &poly1 * &poly2;

        // In negacyclic convolution: x^4 = -1, so x^2 stays as x^2
        let expected = [1, 2, 1, 0]; // 1 + 2x + x^2
        assert_eq!(result.coeffs(), &expected);
    }

    #[test]
    fn test_negacyclic_property() {
        const N: usize = 4;
        let q = find_first_prime_up(10, N);
        let ctx = NttContext::<N>::new(q);

        // Test x^N = -1 property
        // Multiply x^(N-1) by x should give -1 (i.e., q-1)
        let x_to_n_minus_1 =
            NttPolynomial::from_coeffs([0, 0, 0, 1], Arc::clone(&ctx));
        let x = NttPolynomial::from_coeffs([0, 1, 0, 0], Arc::clone(&ctx));

        let result = &x_to_n_minus_1 * &x;

        // x^3 * x = x^4 = -1 ≡ q-1 (mod q)
        let expected = [q - 1, 0, 0, 0];
        assert_eq!(result.coeffs(), &expected);
    }

    #[test]
    fn test_zero_properties() {
        const N: usize = 4;
        let q = find_first_prime_up(10, N);
        let ctx = NttContext::<N>::new(q);

        let zero = NttPolynomial::zero(Arc::clone(&ctx));
        let a = NttPolynomial::sample_random(Arc::clone(&ctx));

        // a + 0 = a
        let a_plus_zero = &a + &zero;
        assert_eq!(a_plus_zero.coeffs(), a.coeffs());

        // a * 0 = 0
        let a_times_zero = &a * &zero;
        assert_eq!(a_times_zero.coeffs(), &[0u64; N]);

        // -0 = 0
        let neg_zero = -&zero;
        assert_eq!(neg_zero.coeffs(), &[0u64; N]);
    }

    #[test]
    fn test_sample_random() {
        const N: usize = 8;
        let q = find_first_prime_up(10, N);
        let ctx = NttContext::<N>::new(q);

        let poly1 = NttPolynomial::sample_random(Arc::clone(&ctx));
        let poly2 = NttPolynomial::sample_random(Arc::clone(&ctx));

        // Very unlikely to be identical
        assert_ne!(poly1.coeffs(), poly2.coeffs());

        // All coefficients should be in valid range [1, q)
        for &coeff in poly1.coeffs() {
            assert!(coeff > 0 && coeff < q);
        }
    }
}
