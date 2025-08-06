/*
Mathematical Background:

For NTT to work, we need:
1. Prime q ≡ 1 (mod 2N) - ensures primitive 2N-th root of unity exists
2. Generator g such that:
   - g^N ≡ -1 (mod q)  [primitive 2N-th root property]
   - g^(2N) ≡ 1 (mod q) [order divides 2N]

The constraint q ≡ 1 (mod 2N) means q = k×2N + 1 for some integer k.
This ensures that 2N divides (q-1), so primitive 2N-th roots exist in Z_q*.

SageMath verification for small example:
```sage
q = 17  # Example: 17 = 2×4×2 + 1, so q ≡ 1 (mod 8) for N=4
N = 4
F = GF(q)
g = F.multiplicative_generator()^((q-1)/(2*N))  # Find primitive 2N-th root
print(f"g = {g}")
print(f"g^N = {g^N} (should be -1) (g^4 = 16 = -1 mod 17)")
print(f"g^(2N) = {g^(2*N)} (should be 1)")
```
*/
use rust_ntt::*;

fn main() {
    println!("🔍 NTT-Friendly Prime Search & Generator Finding");
    println!("══════════════════════════════════════════════════");

    const N: usize = 1usize << 10; // 1024
    let target_bits = 52;

    println!("🎯 Target: Find primes q ≡ 1 (mod 2N) for NTT");
    println!("📐 N (degree) = {} = 2^{}", N, N.trailing_zeros());
    println!("🔢 Target bits ≈ {}", target_bits);
    println!("📏 Constraint: q ≡ 1 (mod 2N = {})\n", 2 * N);

    // Find several NTT-friendly primes
    println!("🔍 Searching for NTT-friendly primes:");
    println!(
        "{:<15} {:<20} {:<10} {:<15}",
        "Type", "Prime q", "Bits", "q mod 2N"
    );
    println!("{}", "─".repeat(65));

    let mut q = find_first_prime_up(target_bits, N);
    println!(
        "{:<15} {:<20} {:<10} {:<15}",
        "first_up",
        q,
        bit_length(q),
        q % (2 * N as u64)
    );

    // Find next few primes up
    for i in 1..=3 {
        q = find_next_prime_up(q, N);
        println!(
            "{:<15} {:<20} {:<10} {:<15}",
            format!("next_up_{}", i),
            q,
            bit_length(q),
            q % (2 * N as u64)
        );
    }

    // Also try finding primes downward from a higher starting point
    let q_down = find_first_prime_down(target_bits + 2, N);
    println!(
        "{:<15} {:<20} {:<10} {:<15}",
        "first_down",
        q_down,
        bit_length(q_down),
        q_down % (2 * N as u64)
    );

    // Pick one prime for detailed analysis
    let chosen_q = q;
    println!("\n🧪 Detailed analysis of q = {}:", chosen_q);
    analyze_prime_and_generator(chosen_q, N);

    // Test with smaller degree for easier verification
    println!("\n📚 Educational example with N = 4:");
    const SMALL_N: usize = 4;
    let small_q = find_first_prime_up(10, SMALL_N);
    analyze_prime_and_generator(small_q, SMALL_N);
}

fn analyze_prime_and_generator(q: u64, n: usize) {
    let class = CongruenceClass::new(q);

    println!("  📊 Prime properties:");
    println!("    q = {} ({} bits)", q, bit_length(q));
    println!("    q - 1 = {} = {} × 2N", q - 1, (q - 1) / (2 * n as u64));
    println!("    Verification: q ≡ {} (mod 2N) ✓", q % (2 * n as u64));

    // Find primitive root and NTT generator
    let g0 = find_primitive_root(q);
    let g = find_generator(q, n);

    println!("  🌱 Generator analysis:");
    println!("    Primitive root g₀ = {}", g0);
    println!("    NTT generator g = g₀^((q-1)/2N) = {}", g);

    // Verify generator properties
    let g_to_n = class.modexp(g, n as u64);
    let g_to_2n = class.modexp(g, (2 * n) as u64);

    println!("  ✅ Generator verification:");
    println!(
        "    g^N = {} = {} (should be q-1 = -1)",
        g_to_n,
        if g_to_n == q - 1 { "q-1 ✓" } else { "❌" }
    );
    println!("    g^(2N) = {} (should be 1)", g_to_2n);

    if g_to_n == q - 1 && g_to_2n == 1 {
        println!("    🎉 Perfect! g is a primitive 2N-th root of unity");
    } else {
        println!("    ❌ Generator verification failed!");
    }

    // Show some powers of g for educational purposes
    if n <= 16 {
        println!("  📈 Powers of g (first {} values):", n.min(8));
        for i in 0..n.min(8) {
            let power = class.modexp(g, i as u64);
            println!("    g^{} = {}", i, power);
        }
    }
}

fn bit_length(n: u64) -> u32 {
    if n == 0 { 0 } else { 64 - n.leading_zeros() }
}
