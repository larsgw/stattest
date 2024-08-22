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
        self.abs()
    }
}

impl Abs for f64 {
    type Output = f64;

    #[inline]
    fn abs(self) -> Self::Output {
        self.abs()
    }
}

impl Abs for &f32 {
    type Output = f32;

    #[inline]
    fn abs(self) -> Self::Output {
        self.abs()
    }
}

impl Abs for &f64 {
    type Output = f64;

    #[inline]
    fn abs(self) -> Self::Output {
        self.abs()
    }
}

impl Abs for i8 {
    type Output = i8;

    #[inline]
    fn abs(self) -> Self::Output {
        self.abs()
    }
}

impl Abs for i16 {
    type Output = i16;

    #[inline]
    fn abs(self) -> Self::Output {
        self.abs()
    }
}

impl Abs for i32 {
    type Output = i32;

    #[inline]
    fn abs(self) -> Self::Output {
        self.abs()
    }
}

impl Abs for i64 {
    type Output = i64;

    #[inline]
    fn abs(self) -> Self::Output {
        self.abs()
    }
}

impl Abs for i128 {
    type Output = i128;

    #[inline]
    fn abs(self) -> Self::Output {
        self.abs()
    }
}

impl Abs for isize {
    type Output = isize;

    #[inline]
    fn abs(self) -> Self::Output {
        self.abs()
    }
}

impl Abs for &i8 {
    type Output = i8;

    #[inline]
    fn abs(self) -> Self::Output {
        self.abs()
    }
}

impl Abs for &i16 {
    type Output = i16;

    #[inline]
    fn abs(self) -> Self::Output {
        self.abs()
    }
}

impl Abs for &i32 {
    type Output = i32;

    #[inline]
    fn abs(self) -> Self::Output {
        self.abs()
    }
}

impl Abs for &i64 {
    type Output = i64;

    #[inline]
    fn abs(self) -> Self::Output {
        self.abs()
    }
}

impl Abs for &i128 {
    type Output = i128;

    #[inline]
    fn abs(self) -> Self::Output {
        self.abs()
    }
}

impl Abs for &isize {
    type Output = isize;

    #[inline]
    fn abs(self) -> Self::Output {
        self.abs()
    }
}
