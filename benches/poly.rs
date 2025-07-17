use criterion::{Criterion, criterion_group, criterion_main};

use rust_ntt::*;

use std::hint::black_box;

const Q: u64 = 741507920154517877;


fn bench_naive_neg_conv(c: &mut Criterion) {
    const N: usize = 1usize<<12;
    let q: u64 = Q.clone();

    let ring = PolyRing::<N>::new(q);

    let ax = ring.sample_random();
    let bx = ring.sample_random();

    c.bench_function("naive convolution ", |b| {
        b.iter(|| { ring.naive_negacyclic_convolution(black_box(&ax), black_box(&bx)); })
    });
}

fn bench_ntt_neg_conv(c: &mut Criterion) {
    const N: usize = 1usize<<12;
    let q: u64 = Q.clone();

    let ring = PolyRing::<N>::new(q);

    let ax = ring.sample_random();
    let bx = ring.sample_random();

    c.bench_function("ntt convolution ", |b| {
        b.iter(|| { ring.ntt_negacyclic_convolution(black_box(&ax), black_box(&bx)); })
    });
}

fn bench_ntt_forward(c: &mut Criterion) {
    const N: usize = 1usize<<12;
    let q: u64 = Q.clone();

    let ring = PolyRing::<N>::new(q);

    let ax = ring.sample_random();

    c.bench_function("ntt forward ", |b| {
        b.iter(|| { ring.ntt_forward(&mut black_box(ax)); })
    });
}

fn bench_ntt_inverse(c: &mut Criterion) {
    const N: usize = 1usize<<12;
    let q: u64 = Q.clone();

    let ring = PolyRing::<N>::new(q);

    let ax = ring.sample_random();

    c.bench_function("ntt inverse ", |b| {
        b.iter(|| { ring.ntt_inverse(&mut black_box(ax)); })
    });
}

criterion_group!(poly, bench_naive_neg_conv, 
                       bench_ntt_neg_conv, 
                       bench_ntt_forward,
                       bench_ntt_inverse
                );
criterion_main!(poly);