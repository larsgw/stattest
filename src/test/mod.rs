//! Defines frequentist statistical tests.

pub use self::f::*;
pub use self::levenes::*;
pub use self::mann_whitney_u::*;
pub use self::shapiro_wilk::*;
pub use self::students_t::*;
pub use self::welchs_t::*;
pub use self::wilcoxon_w::*;

mod f;
mod levenes;
mod mann_whitney_u;
mod shapiro_wilk;
mod students_t;
mod welchs_t;
mod wilcoxon_w;

/// Alternative hypothesis for comparing two means.
pub enum AlternativeHypothesis {
    Greater,
    Different,
    Less,
}
