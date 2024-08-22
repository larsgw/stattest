//! Criterion benchmark to measure the performance of the Wilcoxon signed-rank test.

use criterion::black_box;
use criterion::{criterion_group, criterion_main, Criterion};

use core::ops::Add;
use rand::prelude::SliceRandom;
use stattest::test::WilcoxonWTest;

fn generate_test_cases<const N: usize, F: Default + Copy + Add<F, Output = F>>(
    step: F,
) -> ([F; N], [F; N]) {
    let mut x = [F::default(); N];
    let mut y = [F::default(); N];
    let mut start = F::default();
    for i in 0..N {
        start = start + step;
        x[i] = start;
        start = start + step;        
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
        .map(|_| generate_test_cases::<200_000, f64>(0.1))
        .collect::<Vec<_>>();

    group.bench_function("f64", |b| {
        b.iter(|| {
            for (x, y) in test_cases.iter() {
                WilcoxonWTest::paired(black_box(x), black_box(y)).unwrap();
            }
        })
    });

    let test_cases = (0..10)
        .map(|_| generate_test_cases::<200_000, f32>(0.1))
        .collect::<Vec<_>>();

    group.bench_function("f32", |b| {
        b.iter(|| {
            for (x, y) in test_cases.iter() {
                WilcoxonWTest::paired(black_box(x), black_box(y)).unwrap();
            }
        })
    });

    let test_cases = (0..10)
        .map(|_| generate_test_cases::<200_000, i64>(1))
        .collect::<Vec<_>>();

    group.bench_function("i64", |b| {
        b.iter(|| {
            for (x, y) in test_cases.iter() {
                WilcoxonWTest::paired(black_box(x), black_box(y)).unwrap();
            }
        })
    });
    
    let test_cases = (0..10)
        .map(|_| generate_test_cases::<200_000, i32>(1))
        .collect::<Vec<_>>();

    group.bench_function("i32", |b| {
        b.iter(|| {
            for (x, y) in test_cases.iter() {
                WilcoxonWTest::paired(black_box(x), black_box(y)).unwrap();
            }
        })
    });

    group.finish();
}

criterion_group!(benches, bench_wilcoxon);

criterion_main!(benches);
