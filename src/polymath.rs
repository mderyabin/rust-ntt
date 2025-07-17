// Number Theoretic operations over polynomials and vectors
use crate::ntmath::*;
use rand::{Rng, rng};

#[derive(Debug, Clone, Copy)]
pub struct PolyRing<const N: usize>  {
    pub class: CongruenceClass,

    g: u64,
    invN: u64,
    tf:  [u64; N],
    itf: [u64; N],
}

impl<const N: usize>  PolyRing<N> {
    pub fn new(q: u64) -> Self {
        let class = CongruenceClass::new(q);

        let g = find_generator(q, N);

        let tf = compute_twiddle_factors::<N>(&class, g, false);
        let itf= compute_twiddle_factors::<N>(&class, g, true);
        
        let invN = class.modinv(N as u64);

        Self { class, g, invN, tf, itf }
    }

    pub fn sample_random(&self) -> [u64; N] {
        let mut generator = rng();
        let mut res = [0u64; N];
        for i in 0..N {
            res[i] = generator.random_range(1..self.class.q);
        }
        res
    }

    pub fn add(&self, ax: &[u64], bx: &[u64]) -> [u64; N] {
        assert_eq!(ax.len(), N, "input must have len = {}", N);
        assert_eq!(bx.len(), N, "input must have len = {}", N);

        let mut cx = [0u64; N];

        for i in 0..N {
            cx[i] = self.class.modadd(ax[i], bx[i])
        }

        cx
    }

    pub fn add_eq(&self, ax: &mut [u64], bx: &[u64]) {
        assert_eq!(ax.len(), N, "input must have len = {}", N);
        assert_eq!(bx.len(), N, "input must have len = {}", N);

        for (a, &b) in ax.iter_mut().zip(bx) {
            self.class.modadd_eq(&mut *a, b);
        }
    }

    pub fn sub(&self, ax: &[u64], bx: &[u64]) -> [u64; N] {
        assert_eq!(ax.len(), N, "input must have len = {}", N);
        assert_eq!(bx.len(), N, "input must have len = {}", N);

        let mut cx = [0u64; N];

        for i in 0..N {
            cx[i] = self.class.modsub(ax[i], bx[i])
        }

        cx
    }

    pub fn sub_eq(&self, ax: &mut [u64], bx: &[u64]) {
        assert_eq!(ax.len(), N, "input must have len = {}", N);
        assert_eq!(bx.len(), N, "input must have len = {}", N);

        for (a, &b) in ax.iter_mut().zip(bx) {
            self.class.modsub_eq(&mut *a, b);
        }
    }

    pub fn neg(&self, ax: &[u64]) -> Vec<u64> {
        assert_eq!(ax.len(), N, "input must have len = {}", N);

        ax.iter().map(|&a| self.class.modneg(a)).collect()
    }

    pub fn neg_eq(&self, ax: &mut [u64]) {
        assert_eq!(ax.len(), N, "input must have len = {}", N);

        for a in ax.iter_mut() {
            self.class.modneg_eq(&mut *a);
        }
    }

    pub fn mul(&self, ax: &[u64], bx: &[u64]) -> [u64; N] {
        assert_eq!(ax.len(), N, "input must have len = {}", N);
        assert_eq!(bx.len(), N, "input must have len = {}", N);

        let mut cx = [0u64; N];

        for i in 0..N {
            cx[i] = self.class.modmul(ax[i], bx[i])
        }

        cx
    }

    pub fn mul_eq(&self, ax: &mut [u64], bx: &[u64]) {
        assert_eq!(ax.len(), N, "input must have len = {}", N);
        assert_eq!(bx.len(), N, "input must have len = {}", N);

        for (a, &b) in ax.iter_mut().zip(bx) {
            self.class.modmul_eq(&mut *a, b);
        }
    }

    pub fn naive_negacyclic_convolution(&self, ax: &[u64], bx: &[u64]) -> [u64; N] {
        assert_eq!(ax.len(), N, "input must have len = {}", N);
        assert_eq!(bx.len(), N, "input must have len = {}", N);

        let mut cx = [0u64; N];

        for i in 0..N {
            for j in 0..=i {
                self.class
                    .modadd_eq(&mut cx[i], self.class.modmul(ax[j], bx[i - j]));
            }
            for j in (i + 1)..N {
                self.class
                    .modsub_eq(&mut cx[i], self.class.modmul(ax[j], bx[N + i - j]));
            }
        }

        cx
    }

    pub fn ntt_forward(&self, ax: & mut [u64; N])  {
        // Cooley-Tukey forward negacyclic NTT
        // using algorithm from https://eprint.iacr.org/2016/504.pdf

        let mut t = N >> 1;
        let mut n = 1;
        while n < N {
            for i in 0..n {
                let j1 = 2 * i * t;
                let j2 = j1 + t - 1;
                let s = self.tf[n + i]; 
                for j in j1..=j2 {
                    let u = ax[j];
                    let v = self.class.modmul(ax[j + t], s);
                    ax[j] = self.class.modadd(u, v);
                    ax[j + t] = self.class.modsub(u, v);
                }
            }
            n <<= 1;
            t >>= 1;
        }

    }

    pub fn ntt_inverse(&self, ax: &mut [u64; N]) {
        // Gentelman-Sande inverse negacyclic NTT
        // using algorithm from https://eprint.iacr.org/2016/504.pdf

        let mut t = 1;
        let mut h = N>>1;
        while h > 0 {
            let mut j1 = 0;
            for i in 0..h {
                let j2 = j1 + t - 1;
                let s = self.itf[h + i]; 
                for j in j1..=j2 {
                    let u = ax[j];
                    let v = ax[j + t];
                    ax[j] = self.class.modadd(u, v);
                    ax[j + t] = self.class.modsub(u, v);
                    ax[j + t] = self.class.modmul(ax[j + t], s);
                    // self.class.modmul_eq(&mut ax[j + t], s);
                }
                j1 += t<<1;
            }
            h >>= 1;
            t <<= 1;
        }

        for j in 0..N {
            ax[j] = self.class.modmul(ax[j], self.invN);
            // self.class.modmul_eq(&mut ax[j], self.invN);
        }

    }

    pub fn ntt_negacyclic_convolution(&self, ax: &[u64; N], bx: &[u64; N]) -> [u64; N] {
        let mut ax_ntt= ax.clone();
        let mut bx_ntt = bx.clone();

        self.ntt_forward(&mut ax_ntt);
        self.ntt_forward(&mut bx_ntt);

        let mut cx_ntt = self.mul(&ax_ntt, &bx_ntt);

        self.ntt_inverse(&mut cx_ntt);

        cx_ntt
    }

}

fn bit_reverse(number: usize, length: usize) -> usize {
    let mut reversed = 0;
    for i in 0..length {
        if (number >> i) & 1 != 0 {
            reversed |= 1 << (length - 1 - i);
        }
    }

    reversed
}

fn compute_twiddle_factors<const N: usize>(class: &CongruenceClass, g: u64, is_inverse: bool) -> [u64; N] {
    let mut tf=  [0u64; N];
    let mut tf_direct=  [0u64; N];
    let degree_len = usize::BITS - N.leading_zeros() - 1;

    let gg = if is_inverse { class.modinv(g) } else { g };

    tf_direct[0] = 1;
    for i in 1..N {
        tf_direct[i] = class.modmul(tf_direct[i-1], gg);
    }

    for i in 0..N {

        tf[i] = tf_direct[bit_reverse(i, degree_len as usize)];
    }

    tf
}