//! Submodule providing the Bounded trait.

pub trait Bounded {
    const UPPER_BOUND: Self;
    const LOWER_BOUND: Self;
}

impl Bounded for i8 {
    const UPPER_BOUND: i8 = i8::MAX;
    const LOWER_BOUND: i8 = i8::MIN + 1;
}

impl Bounded for i16 {
    const UPPER_BOUND: i16 = i16::MAX;
    const LOWER_BOUND: i16 = i16::MIN + 1;
}

impl Bounded for i32 {
    const UPPER_BOUND: i32 = i32::MAX;
    const LOWER_BOUND: i32 = i32::MIN + 1;
}

impl Bounded for i64 {
    const UPPER_BOUND: i64 = i64::MAX;
    const LOWER_BOUND: i64 = i64::MIN + 1;
}

impl Bounded for i128 {
    const UPPER_BOUND: i128 = i128::MAX;
    const LOWER_BOUND: i128 = i128::MIN + 1;
}

impl Bounded for f32 {
    const UPPER_BOUND: f32 = f32::MAX;
    const LOWER_BOUND: f32 = f32::MIN;
}

impl Bounded for f64 {
    const UPPER_BOUND: f64 = f64::MAX;
    const LOWER_BOUND: f64 = f64::MIN;
}
