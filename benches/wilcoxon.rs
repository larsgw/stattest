//! Criterion benchmark to measure the performance of the Wilcoxon signed-rank test.

use criterion::black_box;
use criterion::{criterion_group, criterion_main, Criterion};

use rand::prelude::SliceRandom;
use stattest::test::WilcoxonWTest;

fn generate_test_cases<const N: usize>() -> ([f64; N], [f64; N]) {
    let mut x = [0.0; N];
    let mut y = [0.0; N];
    for i in 0..N {
        x[i] = i as f64;
        y[i] = (i as f64) + 0.1;
    }

    // We shuffle the arrays to make the test more realistic
    x.shuffle(&mut rand::thread_rng());
    y.shuffle(&mut rand::thread_rng());

    (x, y)
}

fn bench_wilcoxon(c: &mut Criterion) {
    let test_cases = (0..100)
        .map(|_| generate_test_cases::<1000>())
        .collect::<Vec<_>>();

    c.bench_function("Wilcoxon signed-rank test", |b| {
        b.iter(|| {
            for (x, y) in test_cases.iter() {
                WilcoxonWTest::paired(black_box(x), black_box(y)).unwrap();
            }
        })
    });
}

criterion_group!(benches, bench_wilcoxon);

criterion_main!(benches);
