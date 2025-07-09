fn main() {
    let q = 741507920154517877u64;
    let a = 280429249880250689u64;
    let b = 530127388764165774u64;
    let mu = precompute_mu(q);

    let naive = modmul_naive(a, b, q);
    let fast = modmul_barrett(a, b, q, mu);

    println!("a = {a}, b = {b}");
    println!("naive  = {naive}");
    println!("barrett= {fast}");
    assert_eq!(fast, naive);
}

fn modmul_naive(a: u64, b: u64, q: u64) -> u64 {
    ((a as u128 * b as u128) % (q as u128)) as u64
}

/// μ = floor(2^126 / q), shift = 62
fn precompute_mu(q: u64) -> u128 {
    (1u128 << 126) / q as u128
}

/// (a * b) mod q using Barrett reduction with μ = 2^126 / q
fn modmul_barrett(a: u64, b: u64, q: u64, mu: u128) -> u64 {
    let c = (a as u128) * (b as u128);

    let tmp1 = c >> 62;
    let tmp2 = (tmp1 * mu) >> 64;

    let mut r = c.wrapping_sub(tmp2 * q as u128);

    if r >= q as u128 {
        r -= q as u128;
    }
    if r >= q as u128 {
        r -= q as u128;
    }
    if r >= q as u128 {
        r -= q as u128;
    }

    r as u64
}
