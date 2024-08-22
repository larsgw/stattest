//! Submodule implementing the `Zero` trait for all numeric types.

/// Trait for zero value.
pub trait Zero {
    /// Zero value.
    const ZERO: Self;
}

macro_rules! impl_zero {
    ($($t:ty),*) => {
        $(
            impl Zero for $t {
                const ZERO: Self = 0 as $t;
            }
        )*
    };
}

impl_zero!(u8, u16, u32, u64, u128, usize, i8, i16, i32, i64, i128, isize, f32, f64);
