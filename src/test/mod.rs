//! Defines frequentist statistical tests.

pub use self::f::*;
pub use self::levenes::*;
pub use self::students_t::*;
pub use self::shapiro_wilk::*;

mod f;
mod levenes;
mod students_t;
mod shapiro_wilk;
