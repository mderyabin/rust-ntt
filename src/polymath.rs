// Number Theoretic operations over polynomials and vectors
use crate::ntmath::*;
use rand::{Rng, rng};

#[derive(Debug, Clone, Copy)]
pub struct PolyRing<const N: usize> {
    pub class: CongruenceClass,

    inv_n: u64,
    inv_n_shoup: u64,
    tf: [u64; N],
    tf_shoup: [u64; N],
    itf: [u64; N],
    itf_shoup: [u64; N],
}

impl<const N: usize> PolyRing<N> {
    pub fn new(q: u64) -> Self {
        let class = CongruenceClass::new(q);

        let g = find_generator(q, N);

        let tf = compute_twiddle_factors::<N>(&class, g, false);
        let itf = compute_twiddle_factors::<N>(&class, g, true);

        let mut tf_shoup = [0u64; N];
        tf_shoup
            .iter_mut()
            .enumerate()
            .for_each(|(i, prec)| *prec = class.precompute_shoup(tf[i]));

        let mut itf_shoup = [0u64; N];
        itf_shoup
            .iter_mut()
            .enumerate()
            .for_each(|(i, prec)| *prec = class.precompute_shoup(itf[i]));

        let inv_n = class.modinv(N as u64);
        let inv_n_shoup = class.precompute_shoup(inv_n);

        Self {
            class,
            inv_n,
            inv_n_shoup,
            tf,
            tf_shoup,
            itf,
            itf_shoup,
        }
    }

    pub fn sample_random(&self) -> [u64; N] {
        let mut generator = rng();
        let mut res = [0u64; N];
        res.iter_mut()
            .for_each(|r| *r = generator.random_range(1..self.class.q));
        res
    }

    pub fn add(&self, ax: &[u64], bx: &[u64]) -> [u64; N] {
        assert_eq!(ax.len(), N, "input must have len = {}", N);
        assert_eq!(bx.len(), N, "input must have len = {}", N);

        let mut cx = [0u64; N];
        cx.iter_mut()
            .enumerate()
            .for_each(|(i, c)| *c = self.class.modadd(ax[i], bx[i]));

        cx
    }

    pub fn add_eq(&self, ax: &mut [u64], bx: &[u64]) {
        assert_eq!(ax.len(), N, "input must have len = {}", N);
        assert_eq!(bx.len(), N, "input must have len = {}", N);

        ax.iter_mut().zip(bx).for_each(|(a, &b)| {
            self.class.modadd_eq(a, b);
        });
    }

    pub fn sub(&self, ax: &[u64], bx: &[u64]) -> [u64; N] {
        assert_eq!(ax.len(), N, "input must have len = {}", N);
        assert_eq!(bx.len(), N, "input must have len = {}", N);

        let mut cx = [0u64; N];
        cx.iter_mut()
            .enumerate()
            .for_each(|(i, c)| *c = self.class.modsub(ax[i], bx[i]));

        cx
    }

    pub fn sub_eq(&self, ax: &mut [u64], bx: &[u64]) {
        assert_eq!(ax.len(), N, "input must have len = {}", N);
        assert_eq!(bx.len(), N, "input must have len = {}", N);

        ax.iter_mut().zip(bx).for_each(|(a, &b)| {
            self.class.modsub_eq(a, b);
        });
    }

    pub fn neg(&self, ax: &[u64]) -> Vec<u64> {
        assert_eq!(ax.len(), N, "input must have len = {}", N);

        ax.iter().map(|&a| self.class.modneg(a)).collect()
    }

    pub fn neg_eq(&self, ax: &mut [u64]) {
        assert_eq!(ax.len(), N, "input must have len = {}", N);

        ax.iter_mut().for_each(|a| self.class.modneg_eq(a));
    }

    pub fn mul(&self, ax: &[u64], bx: &[u64]) -> [u64; N] {
        assert_eq!(ax.len(), N, "input must have len = {}", N);
        assert_eq!(bx.len(), N, "input must have len = {}", N);

        let mut cx = [0u64; N];
        cx.iter_mut()
            .enumerate()
            .for_each(|(i, c)| *c = self.class.modmul(ax[i], bx[i]));

        cx
    }

    pub fn mul_eq(&self, ax: &mut [u64], bx: &[u64]) {
        assert_eq!(ax.len(), N, "input must have len = {}", N);
        assert_eq!(bx.len(), N, "input must have len = {}", N);

        ax.iter_mut().zip(bx).for_each(|(a, &b)| {
            self.class.modmul_eq(a, b);
        });
    }

    pub fn naive_negacyclic_convolution(&self, ax: &[u64], bx: &[u64]) -> [u64; N] {
        assert_eq!(ax.len(), N, "input must have len = {}", N);
        assert_eq!(bx.len(), N, "input must have len = {}", N);

        let mut cx = [0u64; N];

        // Вычисляем свертку для каждого выходного элемента
        (0..N).for_each(|i| {
            // Суммируем произведения для j <= i
            (0..=i).for_each(|j| {
                self.class
                    .modadd_eq(&mut cx[i], self.class.modmul(ax[j], bx[i - j]));
            });

            // Вычитаем произведения для j > i (отрицательные индексы)
            ((i + 1)..N).for_each(|j| {
                self.class
                    .modsub_eq(&mut cx[i], self.class.modmul(ax[j], bx[N + i - j]));
            });
        });

        cx
    }

    #[inline]
    pub fn ntt_forward(&self, ax: &mut [u64; N]) {
        // Cooley-Tukey forward negacyclic NTT
        // using algorithm from https://eprint.iacr.org/2016/504.pdf

        let mut t = N >> 1;
        let mut n = 1;

        while n < N {
            // Обрабатываем все группы бабочек на текущем уровне
            (0..n).for_each(|i| {
                let j1 = 2 * i * t;
                let j2 = j1 + t - 1;
                let s = self.tf[n + i];

                // Обрабатываем пары элементов в текущей группе
                (j1..=j2).for_each(|j| {
                    let u = ax[j];
                    let v = self.class.modmul(ax[j + t], s);

                    // Операция бабочки
                    ax[j] = self.class.modadd(u, v);
                    ax[j + t] = self.class.modsub(u, v);
                });
            });

            n <<= 1;
            t >>= 1;
        }
    }

    #[inline]
    pub fn ntt_forward_shoup(&self, ax: &mut [u64; N]) {
        // Cooley-Tukey forward negacyclic NTT
        // using algorithm from https://eprint.iacr.org/2016/504.pdf

        let mut t = N >> 1;
        let mut n = 1;

        while n < N {
            // Обрабатываем все группы бабочек на текущем уровне
            (0..n).for_each(|i| {
                let j1 = 2 * i * t;
                let j2 = j1 + t - 1;
                let s = self.tf[n + i];
                let s_shoup = self.tf_shoup[n + i];

                // Обрабатываем пары элементов в текущей группе
                (j1..=j2).for_each(|j| {
                    // let u = ax[j];
                    let v = self.class.modmul_shoup(ax[j + t], s, s_shoup);

                    // Операция бабочки
                    // ax[j] = self.class.modadd(u, v);
                    // ax[j + t] = self.class.modsub(u, v);
                    ax[j + t] = self.class.modsub(ax[j], v);
                    self.class.modadd_eq(&mut ax[j], v);
                });
            });

            n <<= 1;
            t >>= 1;
        }
    }

    #[inline]
    pub fn ntt_inverse(&self, ax: &mut [u64; N]) {
        // Gentelman-Sande inverse negacyclic NTT
        // using algorithm from https://eprint.iacr.org/2016/504.pdf

        let mut t = 1;
        let mut h = N >> 1;

        while h > 0 {
            let mut j1 = 0;

            // Обрабатываем все группы бабочек на текущем уровне
            (0..h).for_each(|i| {
                let j2 = j1 + t - 1;
                let s = self.itf[h + i];

                // Обрабатываем пары элементов в текущей группе
                (j1..=j2).for_each(|j| {
                    let u = ax[j];
                    let v = ax[j + t];

                    // Операция бабочки
                    ax[j] = self.class.modadd(u, v);
                    ax[j + t] = self.class.modsub(u, v);
                    ax[j + t] = self.class.modmul(ax[j + t], s);
                });

                j1 += t << 1;
            });

            h >>= 1;
            t <<= 1;
        }

        // Финальная нормализация
        ax.iter_mut().for_each(|x| {
            *x = self.class.modmul(*x, self.inv_n);
        });
    }

    #[inline]
    pub fn ntt_inverse_shoup(&self, ax: &mut [u64; N]) {
        // Gentelman-Sande inverse negacyclic NTT
        // using algorithm from https://eprint.iacr.org/2016/504.pdf

        let mut t = 1;
        let mut h = N >> 1;

        while h > 0 {
            let mut j1 = 0;

            // Обрабатываем все группы бабочек на текущем уровне
            (0..h).for_each(|i| {
                let j2 = j1 + t - 1;
                let s = self.itf[h + i];
                let s_shoup = self.itf_shoup[h + i];

                // Обрабатываем пары элементов в текущей группе
                (j1..=j2).for_each(|j| {
                    let u = ax[j];
                    let v = ax[j + t];

                    // Операция бабочки
                    ax[j] = self.class.modadd(u, v);
                    ax[j + t] = self.class.modsub(u, v);
                    // ax[j + t] = self.class.modmul_shoup(ax[j + t], s, s_shoup);
                    self.class.modmul_shoup_eq(&mut ax[j + t], s, s_shoup);
                });

                j1 += t << 1;
            });

            h >>= 1;
            t <<= 1;
        }

        // Финальная нормализация
        ax.iter_mut().for_each(|x| {
            // *x = self.class.modmul_shoup(*x, self.inv_n, self.inv_n_shoup);
            self.class
                .modmul_shoup_eq(&mut *x, self.inv_n, self.inv_n_shoup);
        });
    }

    pub fn ntt_negacyclic_convolution(&self, ax: &[u64; N], bx: &[u64; N]) -> [u64; N] {
        let mut ax_ntt = *ax;
        let mut bx_ntt = *bx;

        self.ntt_forward(&mut ax_ntt);
        self.ntt_forward(&mut bx_ntt);

        let mut cx_ntt = self.mul(&ax_ntt, &bx_ntt);

        self.ntt_inverse(&mut cx_ntt);

        cx_ntt
    }

    pub fn ntt_negacyclic_convolution_shoup(&self, ax: &[u64; N], bx: &[u64; N]) -> [u64; N] {
        let mut ax_ntt = *ax;
        let mut bx_ntt = *bx;

        self.ntt_forward_shoup(&mut ax_ntt);
        self.ntt_forward_shoup(&mut bx_ntt);

        let mut cx_ntt = self.mul(&ax_ntt, &bx_ntt);

        self.ntt_inverse_shoup(&mut cx_ntt);

        cx_ntt
    }
}

fn bit_reverse(number: usize, length: usize) -> usize {
    let mut reversed = 0;
    (0..length).for_each(|i| {
        if (number >> i) & 1 != 0 {
            reversed |= 1 << (length - 1 - i);
        }
    });
    reversed
}

fn compute_twiddle_factors<const N: usize>(
    class: &CongruenceClass,
    g: u64,
    is_inverse: bool,
) -> [u64; N] {
    let mut tf = [0u64; N];
    let mut tf_direct = [0u64; N];
    let degree_len = usize::BITS - N.leading_zeros() - 1;

    let gg = if is_inverse { class.modinv(g) } else { g };

    tf_direct[0] = 1;
    (1..N).for_each(|i| {
        tf_direct[i] = class.modmul(tf_direct[i - 1], gg);
    });

    (0..N).for_each(|i| {
        tf[i] = tf_direct[bit_reverse(i, degree_len as usize)];
    });

    tf
}
