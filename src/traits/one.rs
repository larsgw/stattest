//! Submodule implementing the `One` trait for all numeric types.

/// Trait for One value.
pub trait One {
    /// One value.
    const ONE: Self;
}

macro_rules! impl_one {
    ($($t:ty),*) => {
        $(
            impl One for $t {
                const ONE: Self = 1 as $t;
            }
        )*
    };
}

impl_one!(u8, u16, u32, u64, u128, usize, i8, i16, i32, i64, i128, isize, f32, f64);
