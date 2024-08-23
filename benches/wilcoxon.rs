//! Criterion benchmark to measure the performance of the Wilcoxon signed-rank test.

use criterion::black_box;
use criterion::{criterion_group, criterion_main, Criterion};

use core::ops::{Add, Sub};
use rand::prelude::SliceRandom;
use stattest::test::WilcoxonWTest;
use stattest::traits::Bounded;

trait WrappingAdd<Rhs = Self> {
    type Output;

    fn wrapping_add(self, rhs: Rhs) -> Self::Output;
}

impl<T: Bounded + Sub<T, Output = T> + Add<T, Output = T> + PartialOrd + PartialEq + Copy>
    WrappingAdd for T
{
    type Output = T;

    fn wrapping_add(self, rhs: T) -> Self::Output {
        if self >= T::UPPER_BOUND - rhs {
            T::LOWER_BOUND
        } else {
            self + rhs
        }
    }
}

fn generate_test_cases<const N: usize, F: Default + Copy + WrappingAdd<F, Output = F>>(
    step: F,
) -> ([F; N], [F; N]) {
    let mut x = [F::default(); N];
    let mut y = [F::default(); N];
    let mut start = F::default();
    for i in 0..N {
        start = start.wrapping_add(step);
        x[i] = start;
        start = start.wrapping_add(step);
        y[i] = start;
    }

    // We shuffle the arrays to make the test more realistic
    x.shuffle(&mut rand::thread_rng());
    y.shuffle(&mut rand::thread_rng());

    (x, y)
}

macro_rules! bench_float_wilcoxon {
    ($group:ident, $float:ty, $($quantizer:ty),*) => {
        let test_cases = (0..10)
            .map(|_| generate_test_cases::<200_000, $float>(0.1454829354839453473))
            .collect::<Vec<_>>();

        $group.bench_function(&format!(
            "sort_unstable_{}",
            stringify!($float)
        ), |b| {
            b.iter(|| {
                for (x, y) in test_cases.iter() {
                    WilcoxonWTest::paired(black_box(x), black_box(y)).unwrap();
                }
            })
        });

        #[cfg(feature="voracious_radix_sort")]
        $group.bench_function(&format!(
            "voracious_{}",
            stringify!($float)
        ), |b| {
            b.iter(|| {
                for (x, y) in test_cases.iter() {
                    WilcoxonWTest::voracious_paired(black_box(x), black_box(y)).unwrap();
                }
            })
        });

        $(
            $group.bench_function(
                &format!(
                    "quantized_sort_unstable_{}_to_{}",
                    stringify!($float),
                    stringify!($quantizer))
                    , |b| {
                b.iter(|| {
                    for (x, y) in test_cases.iter() {
                        WilcoxonWTest::quantized_paired::<_, _, $quantizer>(black_box(x), black_box(y)).unwrap();
                    }
                })
            });

            #[cfg(feature="voracious_radix_sort")]
            $group.bench_function(&format!(
                "quantized_voracious_{}_to_{}",
                stringify!($float),
                stringify!($quantizer)
            ), |b| {
                b.iter(|| {
                    for (x, y) in test_cases.iter() {
                        WilcoxonWTest::voracious_quantized_paired::<_, _, $quantizer>(black_box(x), black_box(y)).unwrap();
                    }
                })
            });
        )*
    };
}

fn bench_wilcoxon(c: &mut Criterion) {
    let mut group = c.benchmark_group("Wilcoxon signed-rank test");

    bench_float_wilcoxon!(group, f32, i8, i16);
    bench_float_wilcoxon!(group, f64, i8, i16, i32);

    let test_cases = (0..10)
        .map(|_| generate_test_cases::<200_000, i64>(1))
        .collect::<Vec<_>>();

    group.bench_function("sort_unstable_i64", |b| {
        b.iter(|| {
            for (x, y) in test_cases.iter() {
                WilcoxonWTest::paired(black_box(x), black_box(y)).unwrap();
            }
        })
    });

    #[cfg(feature = "voracious_radix_sort")]
    group.bench_function("voracious_i64", |b| {
        b.iter(|| {
            for (x, y) in test_cases.iter() {
                WilcoxonWTest::voracious_paired(black_box(x), black_box(y)).unwrap();
            }
        })
    });

    let test_cases = (0..10)
        .map(|_| generate_test_cases::<200_000, i32>(1))
        .collect::<Vec<_>>();

    group.bench_function("sort_unstable_i32", |b| {
        b.iter(|| {
            for (x, y) in test_cases.iter() {
                WilcoxonWTest::paired(black_box(x), black_box(y)).unwrap();
            }
        })
    });

    #[cfg(feature = "voracious_radix_sort")]
    group.bench_function("voracious_i32", |b| {
        b.iter(|| {
            for (x, y) in test_cases.iter() {
                WilcoxonWTest::voracious_paired(black_box(x), black_box(y)).unwrap();
            }
        })
    });

    let test_cases = (0..10)
        .map(|_| generate_test_cases::<200_000, i16>(1))
        .collect::<Vec<_>>();

    group.bench_function("sort_unstable_i16", |b| {
        b.iter(|| {
            for (x, y) in test_cases.iter() {
                WilcoxonWTest::paired(black_box(x), black_box(y)).unwrap();
            }
        })
    });

    #[cfg(feature = "voracious_radix_sort")]
    group.bench_function("voracious_i16", |b| {
        b.iter(|| {
            for (x, y) in test_cases.iter() {
                WilcoxonWTest::voracious_paired(black_box(x), black_box(y)).unwrap();
            }
        })
    });

    let test_cases = (0..10)
        .map(|_| generate_test_cases::<200_000, i8>(1))
        .collect::<Vec<_>>();

    group.bench_function("sort_unstable_i8", |b| {
        b.iter(|| {
            for (x, y) in test_cases.iter() {
                WilcoxonWTest::paired(black_box(x), black_box(y)).unwrap();
            }
        })
    });

    #[cfg(feature = "voracious_radix_sort")]
    group.bench_function("voracious_i8", |b| {
        b.iter(|| {
            for (x, y) in test_cases.iter() {
                WilcoxonWTest::voracious_paired(black_box(x), black_box(y)).unwrap();
            }
        })
    });

    group.finish();
}

criterion_group!(benches, bench_wilcoxon);

criterion_main!(benches);
