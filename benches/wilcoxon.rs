//! Criterion benchmark to measure the performance of the Wilcoxon signed-rank test.

use criterion::black_box;
use criterion::{criterion_group, criterion_main, Criterion};

use core::ops::{Add, Sub};
use rand::prelude::SliceRandom;
use stattest::test::WilcoxonWTest;

trait Bounded {
    const UPPER_BOUND: Self;
    const LOWER_BOUND: Self;
}

impl Bounded for i8 {
    const UPPER_BOUND: i8 = i8::MAX;
    const LOWER_BOUND: i8 = i8::MIN;
}

impl Bounded for i16 {
    const UPPER_BOUND: i16 = i16::MAX;
    const LOWER_BOUND: i16 = i16::MIN;
}

impl Bounded for i32 {
    const UPPER_BOUND: i32 = i32::MAX;
    const LOWER_BOUND: i32 = i32::MIN;
}

impl Bounded for i64 {
    const UPPER_BOUND: i64 = i64::MAX;
    const LOWER_BOUND: i64 = i64::MIN;
}

impl Bounded for i128 {
    const UPPER_BOUND: i128 = i128::MAX;
    const LOWER_BOUND: i128 = i128::MIN;
}

impl Bounded for f32 {
    const UPPER_BOUND: f32 = f32::MAX;
    const LOWER_BOUND: f32 = f32::MIN;
}

impl Bounded for f64 {
    const UPPER_BOUND: f64 = f64::MAX;
    const LOWER_BOUND: f64 = f64::MIN;
}

trait WrappingAdd<Rhs = Self> {
    type Output;

    fn wrapping_add(self, rhs: Rhs) -> Self::Output;
}

impl<T: Bounded + Sub<T, Output=T> + Add<T, Output = T> + PartialOrd + PartialEq + Copy> WrappingAdd for T {
    type Output = T;

    fn wrapping_add(self, rhs: T) -> Self::Output {
        if self >= T::UPPER_BOUND - rhs{
            T::LOWER_BOUND
        } else {
            self + rhs
        }
    }
}


fn generate_test_cases<const N: usize, F: Default + Copy + WrappingAdd<F, Output=F>>(
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

fn bench_wilcoxon(c: &mut Criterion) {
    let mut group = c.benchmark_group("Wilcoxon signed-rank test");

    let test_cases = (0..10)
        .map(|_| generate_test_cases::<200_000, f64>(0.1454829354839453473))
        .collect::<Vec<_>>();

    group.bench_function("sort_unstable_f64", |b| {
        b.iter(|| {
            for (x, y) in test_cases.iter() {
                WilcoxonWTest::paired(black_box(x), black_box(y)).unwrap();
            }
        })
    });

    #[cfg(feature="voracious_radix_sort")]
    group.bench_function("voracious_f64", |b| {
        b.iter(|| {
            for (x, y) in test_cases.iter() {
                WilcoxonWTest::voracious_paired(black_box(x), black_box(y)).unwrap();
            }
        })
    });

    let test_cases = (0..10)
        .map(|_| generate_test_cases::<200_000, f32>(0.1454829354839453473))
        .collect::<Vec<_>>();

    group.bench_function("sort_unstable_f32", |b| {
        b.iter(|| {
            for (x, y) in test_cases.iter() {
                WilcoxonWTest::paired(black_box(x), black_box(y)).unwrap();
            }
        })
    });

    #[cfg(feature="voracious_radix_sort")]
    group.bench_function("voracious_f32", |b| {
        b.iter(|| {
            for (x, y) in test_cases.iter() {
                WilcoxonWTest::voracious_paired(black_box(x), black_box(y)).unwrap();
            }
        })
    });

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

    #[cfg(feature="voracious_radix_sort")]
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

    #[cfg(feature="voracious_radix_sort")]
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

    #[cfg(feature="voracious_radix_sort")]
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

    #[cfg(feature="voracious_radix_sort")]
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
