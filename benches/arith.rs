use criterion::{Criterion, criterion_group, criterion_main};
use rand::{Rng, rng};

use rust_ntt::*;

use std::hint::black_box;

const Q: u64 = 741507920154517877;

fn benchmark_modadd_naive(c: &mut Criterion) {
    let mut generator = rng();

    let in1: u64 = generator.random_range(1..Q);
    let in2: u64 = generator.random_range(1..Q);

    c.bench_function("modadd naive", |b| {
        b.iter(|| modadd_naive(black_box(in1), black_box(in2), black_box(Q)))
    });
}

fn benchmark_modadd(c: &mut Criterion) {
    let mut generator = rng();

    let in1: u64 = generator.random_range(1..Q);
    let in2: u64 = generator.random_range(1..Q);

    c.bench_function("modadd", |b| {
        b.iter(|| modadd(black_box(in1), black_box(in2), black_box(Q)))
    });
}

fn benchmark_modadd_struct(c: &mut Criterion) {
    let mut generator = rng();

    let in1: u64 = generator.random_range(1..Q);
    let in2: u64 = generator.random_range(1..Q);

    let barrett64 = CongruenceClass::new(Q);

    c.bench_function("modadd struct", |b| {
        b.iter(|| barrett64.modadd(black_box(in1), black_box(in2)))
    });
}

fn benchmark_modadd_eq_struct(c: &mut Criterion) {
    let mut generator = rng();

    let mut in1: u64 = generator.random_range(1..Q);
    let in2: u64 = generator.random_range(1..Q);

    let barrett64 = CongruenceClass::new(Q);

    c.bench_function("modadd eq struct", |b| {
        b.iter(|| barrett64.modadd_eq(&mut in1, black_box(in2)))
    });
}

fn benchmark_modsub_struct(c: &mut Criterion) {
    let mut generator = rng();

    let in1: u64 = generator.random_range(1..Q);
    let in2: u64 = generator.random_range(1..Q);

    let barrett64 = CongruenceClass::new(Q);

    c.bench_function("modsub struct", |b| {
        b.iter(|| barrett64.modsub(black_box(in1), black_box(in2)))
    });
}

fn benchmark_modsub_eq_struct(c: &mut Criterion) {
    let mut generator = rng();

    let mut in1: u64 = generator.random_range(1..Q);
    let in2: u64 = generator.random_range(1..Q);

    let barrett64 = CongruenceClass::new(Q);

    c.bench_function("modsub eq struct", |b| {
        b.iter(|| barrett64.modsub_eq(&mut in1, black_box(in2)))
    });
}

fn benchmark_modsub(c: &mut Criterion) {
    let mut generator = rng();

    let in1: u64 = generator.random_range(1..Q);
    let in2: u64 = generator.random_range(1..Q);

    c.bench_function("modsub", |b| {
        b.iter(|| modsub(black_box(in1), black_box(in2), black_box(Q)))
    });
}

fn benchmark_modmul_naive(c: &mut Criterion) {
    let mut generator = rng();

    let in1: u64 = generator.random_range(1..Q);
    let in2: u64 = generator.random_range(1..Q);

    c.bench_function("modmul naive", |b| {
        b.iter(|| modmul_naive(black_box(in1), black_box(in2), black_box(Q)))
    });
}

fn benchmark_modmul_barrett_old(c: &mut Criterion) {
    let mut generator = rng();

    let in1: u64 = generator.random_range(1..Q);
    let in2: u64 = generator.random_range(1..Q);

    let (mu, logq) = barrett_precompute_old(Q);

    c.bench_function("modmul barrett old", |b| {
        b.iter(|| {
            modmul_barrett_old(
                black_box(in1),
                black_box(in2),
                black_box(Q),
                black_box(mu),
                black_box(logq),
            )
        })
    });
}

fn benchmark_modmul_barrett_old_eq(c: &mut Criterion) {
    let mut generator = rng();

    let mut in1: u64 = generator.random_range(1..Q);
    let in2: u64 = generator.random_range(1..Q);

    let (mu, logq) = barrett_precompute_old(Q);

    c.bench_function("modmul barrett old eq", |b| {
        b.iter(|| {
            modmul_barrett_old_eq(
                &mut in1,
                black_box(in2),
                black_box(Q),
                black_box(mu),
                black_box(logq),
            )
        })
    });
}

fn benchmark_modmul_barrett(c: &mut Criterion) {
    let mut generator: rand::prelude::ThreadRng = rng();

    let in1: u64 = generator.random_range(1..Q);
    let in2: u64 = generator.random_range(1..Q);

    let mu = barrett_precompute(Q);

    c.bench_function("modmul barrett", |b| {
        b.iter(|| modmul_barrett(black_box(in1), black_box(in2), black_box(Q), black_box(mu)))
    });
}

fn benchmark_modmul_barrett_eq(c: &mut Criterion) {
    let mut generator = rng();

    let mut in1: u64 = generator.random_range(1..Q);
    let in2: u64 = generator.random_range(1..Q);

    let mu = barrett_precompute(Q);

    c.bench_function("modmul barrett eq", |b| {
        b.iter(|| modmul_barrett_eq(&mut in1, black_box(in2), black_box(Q), black_box(mu)))
    });
}

fn benchmark_modmul_barrett_struct(c: &mut Criterion) {
    let mut generator: rand::prelude::ThreadRng = rng();

    let in1: u64 = generator.random_range(1..Q);
    let in2: u64 = generator.random_range(1..Q);

    let barrett64 = CongruenceClass::new(Q);

    c.bench_function("modmul barrett struct", |b| {
        b.iter(|| barrett64.modmul(black_box(in1), black_box(in2)))
    });
}

fn benchmark_modmul_barrett_eq_struct(c: &mut Criterion) {
    let mut generator = rng();

    let mut in1: u64 = generator.random_range(1..Q);
    let in2: u64 = generator.random_range(1..Q);

    let barrett64 = CongruenceClass::new(Q);

    c.bench_function("modmul barrett struct eq", |b| {
        b.iter(|| barrett64.modmul_eq(&mut in1, black_box(in2)))
    });
}

fn benchmark_modmul_shoup_struct(c: &mut Criterion) {
    let mut generator: rand::prelude::ThreadRng = rng();

    let in1: u64 = generator.random_range(1..Q);
    let in2: u64 = generator.random_range(1..Q);

    let class = CongruenceClass::new(Q);

    let prec = class.precompute_shoup(in2);


    c.bench_function("modmul shoup struct", |b| {
        b.iter(|| class.modmul_shoup(black_box(in1), black_box(in2), black_box(prec)))
    });
}

fn benchmark_modmul_shoup_as64_struct(c: &mut Criterion) {
    let mut generator: rand::prelude::ThreadRng = rng();

    let in1: u64 = generator.random_range(1..Q);
    let in2: u64 = generator.random_range(1..Q);

    let class = CongruenceClass::new(Q);

    let prec = class.precompute_shoup(in2);


    c.bench_function("modmul shoup struct as64", |b| {
        b.iter(|| class.modmul_shoup_as64(black_box(in1), black_box(in2), black_box(prec)))
    });
}

fn benchmark_modmul_shoup_eq_struct(c: &mut Criterion) {
    let mut generator: rand::prelude::ThreadRng = rng();

    let mut in1: u64 = generator.random_range(1..Q);
    let in2: u64 = generator.random_range(1..Q);

    let class = CongruenceClass::new(Q);

    let prec = class.precompute_shoup(in2);


    c.bench_function("modmul shoup eq struct", |b| {
        b.iter(|| class.modmul_shoup_eq(&mut in1, black_box(in2), black_box(prec)))
    });
}

criterion_group!(
    arith,
    benchmark_modadd_naive,
    benchmark_modadd,
    benchmark_modsub,
    benchmark_modadd_struct,
    benchmark_modadd_eq_struct,
    benchmark_modsub_struct,
    benchmark_modsub_eq_struct,
    benchmark_modmul_naive,
    benchmark_modmul_barrett_old,
    benchmark_modmul_barrett_old_eq,
    benchmark_modmul_barrett,
    benchmark_modmul_barrett_eq,
    benchmark_modmul_barrett_struct,
    benchmark_modmul_barrett_eq_struct,
    benchmark_modmul_shoup_struct,
    benchmark_modmul_shoup_as64_struct,
    benchmark_modmul_shoup_eq_struct
);
criterion_main!(arith);
