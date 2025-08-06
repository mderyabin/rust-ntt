use rand::{Rng, SeedableRng, rngs::StdRng};
use rust_ntt::*;

fn main() {
    println!("🔢 Barrett Reduction vs Naive Modular Multiplication");
    println!("═══════════════════════════════════════════════════");

    // Test with the original specific case
    let q = 741507920154517877u64;
    let a = 280429249880250689u64;
    let b = 530127388764165774u64;

    println!("📊 Testing specific large values:");
    println!("q = {} (~2^59.7 bits)", q);
    println!("a = {}", a);
    println!("b = {}", b);

    test_modmul_methods(a, b, q);

    println!("\n🎲 Testing with random values:");
    println!("{}", "─".repeat(40));

    // Test with multiple random cases using deterministic seed
    let mut rng = StdRng::seed_from_u64(42);
    let class = CongruenceClass::new(q);

    for i in 1..=5 {
        let a = rng.random_range(1..q);
        let b = rng.random_range(1..q);

        println!("\n🧪 Test case {}:", i);
        println!("a = {}", a);
        println!("b = {}", b);

        let naive = modmul_naive(a, b, q);
        let barrett_lib = class.modmul(a, b);
        let barrett_manual = modmul_barrett_manual(a, b, q);

        println!("naive     = {}", naive);
        println!("barrett   = {} ✓", barrett_lib);
        println!("manual    = {}", barrett_manual);

        // Verify all methods agree
        assert_eq!(naive, barrett_lib, "Library Barrett failed!");
        assert_eq!(naive, barrett_manual, "Manual Barrett failed!");
        println!("✅ All methods match");
    }

    println!("\n🚀 Performance note: Barrett reduction avoids expensive division");
    println!(
        "   by precomputing μ = ⌊2^126/q⌋ and using bit shifts + multiplication"
    );
    println!("   Useful for repeated modular multiplications with same modulus.");
}

fn test_modmul_methods(a: u64, b: u64, q: u64) {
    let class = CongruenceClass::new(q);

    let naive = modmul_naive(a, b, q);
    let barrett_lib = class.modmul(a, b); // Uses library's Barrett implementation
    let barrett_manual = modmul_barrett_manual(a, b, q);

    println!("\n🧮 Results:");
    println!("naive     = {}", naive);
    println!("barrett   = {}", barrett_lib);
    println!("manual    = {}", barrett_manual);

    // Verify correctness
    assert_eq!(barrett_lib, naive, "Library Barrett mismatch!");
    assert_eq!(barrett_manual, naive, "Manual Barrett mismatch!");
    println!("✅ All implementations match!");
}

/// Reference implementation using 128-bit arithmetic
fn modmul_naive(a: u64, b: u64, q: u64) -> u64 {
    ((a as u128 * b as u128) % (q as u128)) as u64
}

/// Manual Barrett reduction implementation for demonstration
/// μ = ⌊2^126 / q⌋, using shift = 62
fn modmul_barrett_manual(a: u64, b: u64, q: u64) -> u64 {
    let mu = precompute_mu(q);
    let c = (a as u128) * (b as u128);

    // Barrett reduction algorithm:
    // 1. Estimate quotient: q_est ≈ (ab >> 62) * μ >> 64
    let tmp1 = c >> 62;
    let tmp2 = (tmp1 * mu) >> 64;

    // 2. Compute remainder: r = ab - q_est * q
    let mut r = c.wrapping_sub(tmp2 * q as u128);

    // 3. Final corrections (at most 2-3 needed)
    while r >= q as u128 {
        r -= q as u128;
    }

    r as u64
}

/// Precompute Barrett parameter: μ = ⌊2^126 / q⌋
fn precompute_mu(q: u64) -> u128 {
    (1u128 << 126) / q as u128
}
