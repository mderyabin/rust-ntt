use rand::{rng, Rng};

use ntt::*;

fn main() {
    let mut generator = rng();
    let q : u64 = 741507920154517877;


    let ring = PolyRing::new(q, 4);

    let ax = vec![1u64, 2, 3, 4];
    let bx = vec![5u64, 6, 7, 8];

    let cx = ring.naive_negacyclic_convolution(&ax, &bx);

    dbg!(&cx);

    // let a : u64 = generator.random_range(1..q);
    // let b : u64 = generator.random_range(1..q);

    // println!("a = {}", a);
    // println!("b = {}", b);


    // println!("checking addition");

    // let c1 = modadd(a, b, q);
    // let c2 = modadd_naive(a, b, q);

    // println!("c1 = {}", c1);
    // println!("c2 = {}", c2);


    // println!("checking subtraction");

    // let c1 = modsub(a, b, q);
    // let c2 = (q+a-b) % q;

    // println!("c1 = {}", c1);
    // println!("c2 = {}", c2);

    
    // println!("checking barrett multiplication");

    // let mu = barrett_precompute(q);

    // let c1 = modmul_barrett(a, b, q, mu);
    // let c2 = modmul_naive(a, b, q);

    // println!("c1 = {}", c1);
    // println!("c2 = {}", c2);

}
