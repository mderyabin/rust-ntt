use rand::{SeedableRng, rngs::StdRng};
use rust_ntt::*;
use std::sync::Arc;

fn main() {
    const N: usize = 4;
    let q = find_first_prime_up(10, N); // Find prime q == 1 (mod 2N)
    let ctx = NttContext::<N>::new(q);

    println!("🌟 Negacyclic Convolution Demo");
    println!("═══════════════════════════════");
    println!("📐 Degree: N = {}", N);
    println!("🔢 Modulus: q = {} (prime, q ≡ 1 mod {})", q, 2 * N);
    println!("🎯 Ring: Z_q[x]/(x^N + 1)\n");

    let mut rng = StdRng::seed_from_u64(42); // Deterministic seed

    let ax = NttPolynomial::sample_random(Arc::clone(&ctx), &mut rng);
    let bx = NttPolynomial::sample_random(Arc::clone(&ctx), &mut rng);

    println!("📊 Input polynomials:");
    println!("a(x) = {:?}", ax.coeffs());
    println!("b(x) = {:?}", bx.coeffs());

    // Naive O(N^2) convolution for reference
    let cx_naive = ax.naive_negacyclic_convolution(&bx);
    println!("\n🧮 Naive convolution result:");
    println!("c(x) = a(x) * b(x) = {:?}", cx_naive.coeffs());

    // ⚡ Fast NTT-based convolution
    let cx_ntt = ax.negacyclic_convolution(&bx);
    println!("\n🚀 NTT convolution result:");
    println!("c(x) = a(x) * b(x) = {:?}", cx_ntt.coeffs());

    // Verify both methods give same result
    let results_match = cx_naive.coeffs() == cx_ntt.coeffs();
    println!(
        "\n✨ Results match: {}",
        if results_match { "✅" } else { "❌" }
    );
}

/*
Sage math example, can copy paste!

q = 1033
N = 4
R.<x> = PolynomialRing(GF(q))
S = R.quotient(x^N + 1)

# Define polynomials from Rust output
a_coeffs = [544, 561, 657, 419]
b_coeffs = [36, 429, 762, 877]

# Create polynomials in quotient ring
a = S(a_coeffs)
b = S(b_coeffs)

# Compute negacyclic convolution
c = a * b
print("SageMath result:", list(c))

# Let's also verify the negacyclic property manually
print("\nManual verification:")
print("a(x) =", a)
print("b(x) =", b)
print("c(x) = a(x) * b(x) =", c)

# Show that x^4 = -1 in this ring
print(f"\nRing property: x^{N} =", S(x^N), "= -1 mod", q)

# Results;

SageMath result: [30, 631, 453, 128]

Manual verification:
a(x) = 419*xbar^3 + 657*xbar^2 + 561*xbar + 544
b(x) = 877*xbar^3 + 762*xbar^2 + 429*xbar + 36
c(x) = a(x) * b(x) = 128*xbar^3 + 453*xbar^2 + 631*xbar + 30

Ring property: x^4 = 1032 = -1 mod 1033
*/
