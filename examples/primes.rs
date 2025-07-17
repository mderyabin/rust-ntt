use rust_ntt::*;

fn main() {
    let n: usize = 1usize << 10;
    let mut q = find_first_prime_up(52, n);

    let class = CongruenceClass::new(q);

    println!("n = {}", n);
    println!("first prime = {}", q);
    println!("check prime = {}", q % (n as u64));

    // q = find_first_prime_down(10, n);
    // println!("first prime = {}", q);
    // println!("check prime = {}", q % (n as u64));

    q = find_next_prime_up(q, n);
    println!("first prime = {}", q);
    println!("check prime = {}", q % (n as u64));

    // q = find_next_prime_down(q, n);
    // println!("first prime = {}", q);
    // println!("check prime = {}", q % (n as u64));

    // q = find_next_prime_down(q, n);
    // println!("first prime = {}", q);
    // println!("check prime = {}", q % (n as u64));

    // q = find_next_prime_down(q, n);
    // println!("first prime = {}", q);
    // println!("check prime = {}", q % (n as u64));

    let g0 = find_primitive_root(q);
    println!("root = {}", g0);

    let g = find_generator(q, n);
    println!("generator = {}", g);

    println!("check 1 = {}", class.modexp(g, n as u64));
    println!("check 2 = {}", class.modexp(g, (n << 1) as u64));
}
