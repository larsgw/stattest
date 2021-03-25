//! Defines interfaces for creating and approximating statistical distributions.

pub use self::shapiro_wilk::*;
pub use self::signed_rank::*;

mod shapiro_wilk;
mod signed_rank;
