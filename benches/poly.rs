use concrete_ntt::prime64::Plan;
use criterion::{Criterion, criterion_group, criterion_main};
use rand::{SeedableRng, rngs::StdRng};
use rust_ntt::*;
use std::hint::black_box;
use std::sync::Arc;

// const Q: u64 = 741507920154517877;
const N: usize = 1usize << 12;

fn bench_naive_neg_conv(c: &mut Criterion) {
    let q: u64 = find_first_prime_down(58, N);
    let mut rng = StdRng::seed_from_u64(42); // Deterministic seed

    let ctx = NttContext::<N>::new(q);

    let ax = NttPolynomial::sample_random(Arc::clone(&ctx), &mut rng);
    let bx = NttPolynomial::sample_random(Arc::clone(&ctx), &mut rng);

    c.bench_function("naive convolution ", |b| {
        b.iter(|| {
            ax.negacyclic_convolution(black_box(&bx));
        })
    });
}

fn bench_ntt_neg_conv(c: &mut Criterion) {
    let q: u64 = find_first_prime_down(58, N);
    let ctx = NttContext::<N>::new(q);
    let mut rng = StdRng::seed_from_u64(42); // Deterministic seed

    let ax = NttPolynomial::sample_random(Arc::clone(&ctx), &mut rng);
    let bx = NttPolynomial::sample_random(Arc::clone(&ctx), &mut rng);

    c.bench_function("ntt convolution ", |b| {
        b.iter(|| {
            ax.negacyclic_convolution(black_box(&bx));
        })
    });
}

fn bench_ntt_neg_conv_shoup(c: &mut Criterion) {
    let q: u64 = find_first_prime_down(58, N);
    let mut rng = StdRng::seed_from_u64(42); // Deterministic seed
    let ctx = NttContext::<N>::new(q);

    let ax = NttPolynomial::sample_random(Arc::clone(&ctx), &mut rng);
    let bx = NttPolynomial::sample_random(Arc::clone(&ctx), &mut rng);

    c.bench_function("ntt convolution shoup", |b| {
        b.iter(|| {
            ax.negacyclic_convolution_shoup(black_box(&bx));
        })
    });
}

fn bench_ntt_forward(c: &mut Criterion) {
    let q: u64 = find_first_prime_down(58, N);
    let mut rng = StdRng::seed_from_u64(42); // Deterministic seed
    let ctx = NttContext::<N>::new(q);

    let ax = NttPolynomial::sample_random(Arc::clone(&ctx), &mut rng);

    c.bench_function("ntt forward barrett", |b| {
        b.iter(|| {
            let mut poly = ax.clone();
            poly.ntt_forward();
            black_box(poly);
        })
    });
}

fn bench_ntt_inverse(c: &mut Criterion) {
    let q: u64 = find_first_prime_down(58, N);
    let mut rng = StdRng::seed_from_u64(42); // Deterministic seed
    let ctx = NttContext::<N>::new(q);

    let mut ax = NttPolynomial::sample_random(Arc::clone(&ctx), &mut rng);
    ax.ntt_forward(); // Start with NTT-transformed data

    c.bench_function("ntt inverse barrett", |b| {
        b.iter(|| {
            let mut poly = ax.clone();
            poly.ntt_inverse();
            black_box(poly);
        })
    });
}

fn bench_ntt_forward_shoup(c: &mut Criterion) {
    let q: u64 = find_first_prime_down(58, N);
    let mut rng = StdRng::seed_from_u64(42); // Deterministic seed
    let ctx = NttContext::<N>::new(q);

    let mut ax = NttPolynomial::sample_random(Arc::clone(&ctx), &mut rng);
    ax.ntt_forward(); // Start with NTT-transformed data

    c.bench_function("ntt forward shoup", |b| {
        b.iter(|| {
            let mut poly = ax.clone();
            poly.ntt_forward_shoup();
            black_box(poly);
        })
    });
}

fn bench_ntt_inverse_shoup(c: &mut Criterion) {
    let q: u64 = find_first_prime_down(58, N);
    let mut rng = StdRng::seed_from_u64(42); // Deterministic seed
    let ctx = NttContext::<N>::new(q);

    let mut ax = NttPolynomial::sample_random(Arc::clone(&ctx), &mut rng);
    ax.ntt_forward(); // Start with NTT-transformed data

    c.bench_function("ntt inverse shoup", |b| {
        b.iter(|| {
            let mut poly = ax.clone();
            poly.ntt_inverse_shoup();
            black_box(poly);
        })
    });
}

fn bench_concrete_forward(c: &mut Criterion) {
    let q: u64 = find_first_prime_down(58, N);
    let mut rng = StdRng::seed_from_u64(42); // Deterministic seed
    let ctx = NttContext::<N>::new(q);

    let plan = Plan::try_new(N, q).unwrap();
    let sample_poly = NttPolynomial::sample_random(Arc::clone(&ctx), &mut rng);
    let mut ax = sample_poly.coeffs().to_vec();

    c.bench_function("concrete forward ", |b| {
        b.iter(|| {
            plan.fwd(&mut ax);
        })
    });
}

fn bench_concrete_inverse(c: &mut Criterion) {
    let q: u64 = find_first_prime_down(58, N);
    let mut rng = StdRng::seed_from_u64(42); // Deterministic seed
    let ctx = NttContext::<N>::new(q);

    let plan = Plan::try_new(N, q).unwrap();
    let sample_poly = NttPolynomial::sample_random(Arc::clone(&ctx), &mut rng);
    let mut ax = sample_poly.coeffs().to_vec();

    c.bench_function("concrete inverse ", |b| {
        b.iter(|| {
            plan.inv(&mut ax);
        })
    });
}

criterion_group!(
    poly,
    bench_naive_neg_conv,
    bench_ntt_neg_conv,
    bench_ntt_neg_conv_shoup,
    bench_ntt_forward,
    bench_ntt_inverse,
    bench_ntt_forward_shoup,
    bench_ntt_inverse_shoup,
    bench_concrete_forward,
    bench_concrete_inverse,
);
criterion_main!(poly);
