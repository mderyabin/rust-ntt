use rand::{SeedableRng, rngs::StdRng};
use rust_ntt::*;
use std::sync::Arc;

const N: usize = 1usize << 2; // 4

#[test]
fn test_polynomial_addition() {
    let q = find_first_prime_up(20, N);
    let ctx = NttContext::<N>::new(q);
    let mut rng = StdRng::seed_from_u64(42);

    let ax = NttPolynomial::sample_random(Arc::clone(&ctx), &mut rng);
    let bx = NttPolynomial::sample_random(Arc::clone(&ctx), &mut rng);

    // Test addition operator
    let result = &ax + &bx;

    // Verify coefficient-wise with naive modular addition
    let mut expected = [0u64; N];
    for i in 0..N {
        expected[i] = modadd_naive(ax.coeffs()[i], bx.coeffs()[i], q);
    }

    assert_eq!(result.coeffs(), &expected);
}

#[test]
fn test_polynomial_addition_assign() {
    let q = find_first_prime_up(20, N);
    let ctx = NttContext::<N>::new(q);
    let mut rng = StdRng::seed_from_u64(42);

    let mut ax = NttPolynomial::sample_random(Arc::clone(&ctx), &mut rng);
    let bx = NttPolynomial::sample_random(Arc::clone(&ctx), &mut rng);

    // Calculate expected result before mutation
    let mut expected = [0u64; N];
    for i in 0..N {
        expected[i] = modadd_naive(ax.coeffs()[i], bx.coeffs()[i], q);
    }

    // Test += operator
    ax += &bx;

    assert_eq!(ax.coeffs(), &expected);
}

#[test]
fn test_polynomial_pointwise_multiplication() {
    let q = find_first_prime_up(20, N);
    let ctx = NttContext::<N>::new(q);
    let mut rng = StdRng::seed_from_u64(42);

    // Create polynomials and transform to NTT domain for pointwise ops
    let mut ax = NttPolynomial::sample_random(Arc::clone(&ctx), &mut rng);
    let mut bx = NttPolynomial::sample_random(Arc::clone(&ctx), &mut rng);

    ax.ntt_forward();
    bx.ntt_forward();

    // Pointwise multiplication in NTT domain
    let mut result = ax.clone();
    for i in 0..N {
        result.coeffs_mut()[i] = ctx.class().modmul(ax.coeffs()[i], bx.coeffs()[i]);
    }

    // Verify against naive pointwise multiplication
    let mut expected = [0u64; N];
    for i in 0..N {
        expected[i] = modmul_naive(ax.coeffs()[i], bx.coeffs()[i], q);
    }

    assert_eq!(result.coeffs(), &expected);
}

#[test]
fn test_ntt_inverse_is_identity() {
    let q = find_first_prime_down(58, N);
    let ctx = NttContext::<N>::new(q);
    let mut rng = StdRng::seed_from_u64(42);

    for _ in 0..10 {
        let original = NttPolynomial::sample_random(Arc::clone(&ctx), &mut rng);
        let mut test_poly = original.clone();

        // Forward then inverse should give back original
        test_poly.ntt_forward();
        test_poly.ntt_inverse();

        assert_eq!(original.coeffs(), test_poly.coeffs());
    }
}

#[test]
fn test_negacyclic_convolution() {
    let q = find_first_prime_down(58, N);
    let ctx = NttContext::<N>::new(q);
    let mut rng = StdRng::seed_from_u64(42);

    for _ in 0..10 {
        let ax = NttPolynomial::sample_random(Arc::clone(&ctx), &mut rng);
        let bx = NttPolynomial::sample_random(Arc::clone(&ctx), &mut rng);

        // Fast NTT convolution
        let cx_ntt = ax.negacyclic_convolution(&bx);

        // Reference naive convolution
        let cx_naive = ax.naive_negacyclic_convolution(&bx);

        assert_eq!(cx_ntt.coeffs(), cx_naive.coeffs());
    }
}

#[test]
fn test_ntt_shoup_inverse_is_identity() {
    let q = find_first_prime_down(58, N);
    let ctx = NttContext::<N>::new(q);
    let mut rng = StdRng::seed_from_u64(42);

    for _ in 0..10 {
        let original = NttPolynomial::sample_random(Arc::clone(&ctx), &mut rng);
        let mut test_poly = original.clone();

        // Shoup version: forward then inverse should give back original
        test_poly.ntt_forward_shoup();
        test_poly.ntt_inverse_shoup();

        assert_eq!(original.coeffs(), test_poly.coeffs());
    }
}

#[test]
fn test_negacyclic_convolution_shoup() {
    let q = find_first_prime_down(58, N);
    let ctx = NttContext::<N>::new(q);
    let mut rng = StdRng::seed_from_u64(42);

    for _ in 0..10 {
        let ax = NttPolynomial::sample_random(Arc::clone(&ctx), &mut rng);
        let bx = NttPolynomial::sample_random(Arc::clone(&ctx), &mut rng);

        // Shoup-optimized NTT convolution
        let cx_shoup = ax.negacyclic_convolution_shoup(&bx);

        // Reference naive convolution
        let cx_naive = ax.naive_negacyclic_convolution(&bx);

        assert_eq!(cx_shoup.coeffs(), cx_naive.coeffs());
    }
}

#[test]
fn test_convolution_methods_equivalent() {
    let q = find_first_prime_down(58, N);
    let ctx = NttContext::<N>::new(q);
    let mut rng = StdRng::seed_from_u64(42);

    for _ in 0..5 {
        let ax = NttPolynomial::sample_random(Arc::clone(&ctx), &mut rng);
        let bx = NttPolynomial::sample_random(Arc::clone(&ctx), &mut rng);

        // All three methods should give identical results
        let cx_naive = ax.naive_negacyclic_convolution(&bx);
        let cx_ntt = ax.negacyclic_convolution(&bx);
        let cx_shoup = ax.negacyclic_convolution_shoup(&bx);
        let cx_mul_op = &ax * &bx; // Using * operator

        assert_eq!(cx_ntt.coeffs(), cx_naive.coeffs());
        assert_eq!(cx_shoup.coeffs(), cx_naive.coeffs());
        assert_eq!(cx_mul_op.coeffs(), cx_naive.coeffs());
    }
}
