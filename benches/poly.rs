use criterion::{Criterion, criterion_group, criterion_main};

use rust_ntt::*;

use std::hint::black_box;

use concrete_ntt::prime64::Plan;

// const Q: u64 = 741507920154517877;
const N: usize = 1usize<<12;


fn bench_naive_neg_conv(c: &mut Criterion) {
    let q: u64 = find_first_prime_down(58, N);

    let ring = PolyRing::<N>::new(q);

    let ax = ring.sample_random();
    let bx = ring.sample_random();

    c.bench_function("naive convolution ", |b| {
        b.iter(|| { ring.naive_negacyclic_convolution(black_box(&ax), black_box(&bx)); })
    });
}

fn bench_ntt_neg_conv(c: &mut Criterion) {
    let q: u64 = find_first_prime_down(58, N);

    let ring = PolyRing::<N>::new(q);

    let ax = ring.sample_random();
    let bx = ring.sample_random();

    c.bench_function("ntt convolution ", |b| {
        b.iter(|| { ring.ntt_negacyclic_convolution(black_box(&ax), black_box(&bx)); })
    });
}

fn bench_ntt_forward(c: &mut Criterion) {
    let q: u64 = find_first_prime_down(58, N);

    let ring = PolyRing::<N>::new(q);

    let ax = ring.sample_random();

    c.bench_function("ntt forward ", |b| {
        b.iter(|| { ring.ntt_forward(&mut black_box(ax)); })
    });
}

fn bench_ntt_inverse(c: &mut Criterion) {
    let q: u64 = find_first_prime_down(58, N);

    let ring = PolyRing::<N>::new(q);

    let mut ax = ring.sample_random();

    c.bench_function("ntt inverse ", |b| {
        b.iter(|| { ring.ntt_inverse(&mut ax); })
    });
}

fn bench_concrete_forward(c: &mut Criterion) {
    let q: u64 = find_first_prime_down(58, N);

    let plan = Plan::try_new(N, q).unwrap();
    let ring = PolyRing::<N>::new(q);

    let mut ax = ring.sample_random().to_vec();

    c.bench_function("concrete forward ", |b| {
        b.iter(|| { plan.fwd(&mut ax); })
    });
}

fn bench_concrete_inverse(c: &mut Criterion) {
    let q: u64 = find_first_prime_down(58, N);

    let plan = Plan::try_new(N, q).unwrap();
    let ring = PolyRing::<N>::new(q);

    let mut ax = ring.sample_random().to_vec();

    c.bench_function("concrete inverse ", |b| {
        b.iter(|| { plan.inv(&mut ax); })
    });
}


criterion_group!(poly, bench_naive_neg_conv, 
                       bench_ntt_neg_conv, 
                       bench_ntt_forward,
                       bench_ntt_inverse,
                       bench_concrete_forward,
                       bench_concrete_inverse
                );
criterion_main!(poly);