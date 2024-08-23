//! Submodule providing the Quantize trait, which given a value V, its minimum expected value and minimum
//! expected value, and a second value Q which implements Bounded, returns the quantized value of V as Q.
//!
//! For example, a type f64 with expected values between -1.0 and 1.0, when quantized to i8, should return
//! -128 for -1.0, 0 for 0.0, and 127 for 1.0.

use super::bounded::Bounded;
use super::zero::Zero;
use core::fmt::Debug;
use core::ops::{Add, Div, Mul, Neg, Sub};

mod private {
    /// Sealed Trait for primitive types.
    pub trait ToPrimitive<T> {
        /// Convert to primitive type.
        ///
        /// # Safety
        /// This function is unsafe because it may lead to undefined behavior if the value is not
        /// representable in the primitive type, and truncation may occur.
        unsafe fn to_primitive(self) -> T;
    }

    macro_rules! impl_to_primitive {
    ($source:ty, $($target:ty),*) => {
        $(
            impl ToPrimitive<$target> for $source {
                #[inline]
                #[must_use]
                #[expect(
                    clippy::cast_possible_truncation,
                    reason = "Conversion from float to integer should only be used in the context of quantization."
                )]
                unsafe fn to_primitive(self) -> $target {
                    self as $target
                }
            }
        )*
    };
}

    impl_to_primitive!(f32, i8, i16, i32, i64);
    impl_to_primitive!(f64, i8, i16, i32, i64);
}
/// Sealed Trait for non-shifted quantization.
pub trait Quantize<T> {
    /// Quantize values without shifting.
    fn quantize(value: T, factor: T) -> Self;
}

impl<T, Q> Quantize<T> for Q
where
    Q: Bounded + Copy + Zero + Neg<Output = Q> + PartialEq + Debug,
    T: From<Q>
        + Zero
        + private::ToPrimitive<Q>
        + Copy
        + Div<T, Output = T>
        + Sub<T, Output = T>
        + Add<T, Output = T>
        + Neg<Output = T>
        + PartialOrd
        + PartialEq
        + Mul<T, Output = T>,
{
    #[inline]
    #[allow(unsafe_code)]
    #[must_use]
    fn quantize(value: T, factor: T) -> Q {
        debug_assert_eq!(
            Q::UPPER_BOUND,
            -Q::LOWER_BOUND,
            "Quantization must be symmetric."
        );
        debug_assert!(factor > T::ZERO, "Factor must be greater than zero.");
        debug_assert!(
            value >= (-factor) && value <= factor,
            "Value must be within the range of minimum and maximum values."
        );

        // If the value is zero, return the zero value of Q.
        if value == T::ZERO {
            return Q::ZERO;
        }

        let value = value / factor;
        let value = value * T::from(Q::UPPER_BOUND);

        unsafe { value.to_primitive() }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    macro_rules! test_quantize {
        ($source:ty, $($target:ty),*) => {
            $(
                paste::paste! {
                    #[test]
                    fn [<test_quantize_ $source _to_ $target>]() {
                        assert_eq!(<$target as Quantize<$source>>::quantize(-1.0, 1.0), $target::LOWER_BOUND);
                        assert_eq!(<$target as Quantize<$source>>::quantize(-1.0, 5.0), $target::LOWER_BOUND / 5);
                        assert_eq!(<$target as Quantize<$source>>::quantize(-0.5, 1.0), $target::LOWER_BOUND / 2);
                        assert_eq!(<$target as Quantize<$source>>::quantize(-0.5, 5.0), $target::LOWER_BOUND / 10);
                        assert_eq!(<$target as Quantize<$source>>::quantize(0.0, 1.0), 0);
                        assert_eq!(<$target as Quantize<$source>>::quantize(0.0, 5.0), 0);
                        assert_eq!(<$target as Quantize<$source>>::quantize(0.5, 1.0), $target::UPPER_BOUND / 2);
                        assert_eq!(<$target as Quantize<$source>>::quantize(0.5, 5.0), $target::UPPER_BOUND / 10);
                        assert_eq!(<$target as Quantize<$source>>::quantize(1.0, 1.0), $target::UPPER_BOUND);
                        assert_eq!(<$target as Quantize<$source>>::quantize(1.0, 5.0), $target::UPPER_BOUND / 5);
                    }
                }
            )*
        }
    }

    test_quantize!(f64, i8, i16, i32);
    test_quantize!(f32, i8, i16);
}
