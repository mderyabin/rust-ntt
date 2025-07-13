use criterion::{Criterion, criterion_group, criterion_main};

use rust_ntt::*;

use std::hint::black_box;

const Q: u64 = 741507920154517877;


fn bench_naive_neg_conv(c: &mut Criterion) {
    let q: u64 = Q.clone();

    let ring = PolyRing::new(q, 1usize << 12);

    let ax = ring.sample_random();
    let bx = ring.sample_random();

    c.bench_function("modmul barrett struct eq", |b| {
        b.iter(|| ring.naive_negacyclic_convolution(&black_box(&ax), &black_box(&bx)))
    });
}

criterion_group!(poly, bench_naive_neg_conv);
criterion_main!(poly);