//! Defines frequentist statistical tests.

pub use self::f::*;
pub use self::levenes::*;
pub use self::students_t::*;

mod f;
mod levenes;
mod students_t;
