//! Submodule providing the trait Abs and a blanket implementation for all types that implement Zero.

/// Trait for absolute value.
pub trait Abs {
    type Output;

    /// Absolute value.
    fn abs(self) -> Self::Output;
}

impl Abs for f32 {
    type Output = f32;

    #[inline]
    fn abs(self) -> Self::Output {
        f32::abs(self)
    }
}

impl Abs for f64 {
    type Output = f64;

    #[inline]
    fn abs(self) -> Self::Output {
        f64::abs(self)
    }
}

impl Abs for i8 {
    type Output = i8;

    #[inline]
    fn abs(self) -> Self::Output {
        i8::abs(self)
    }
}

impl Abs for i16 {
    type Output = i16;

    #[inline]
    fn abs(self) -> Self::Output {
        i16::abs(self)
    }
}

impl Abs for i32 {
    type Output = i32;

    #[inline]
    fn abs(self) -> Self::Output {
        i32::abs(self)
    }
}

impl Abs for i64 {
    type Output = i64;

    #[inline]
    fn abs(self) -> Self::Output {
        i64::abs(self)
    }
}
