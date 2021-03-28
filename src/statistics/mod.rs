//! Provides traits for statistical computation.

pub use self::iter_statistics_ext::*;
pub use self::ranks::*;
pub use self::statistics_ext::*;

mod iter_statistics_ext;
mod ranks;
mod statistics_ext;
