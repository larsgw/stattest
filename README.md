# Stat-Test

This crate provides statistical tests and relevant utilities.

---

## Status

Currently a number of tests are implemented for comparing the means of two
independent samples, and for checking assumptions regarding the former tests.
However, these tests are all two-sample tests that do not yet support specifying
alternative hypotheses.

### Comparison of independent means

  - **Student's t-test**  
    `stattest::test::StudentsTTest`  
    *Assumptions:* normality, homogeneity of variances  

  - **Welch's t-test**  
    `stattest::test::WelchsTTest`  
    *Assumptions:* normality

  - **Mann-Whitney U-test/Wilcoxon rank-sum test**  
    `stattest::test::MannWhitneyUTest`  
    *Assumptions:* –

### Comparison of paired observations

  - **Student's t-test**  
    `stattest::test::StudentsTTest`  
    *Assumptions:* normality, homogeneity of variances  

  - **Wilcoxon signed rank test**  
    `stattest::test::WilcoxonWTest`  
    *Assumptions:* –

### Assumption tests

  - **Levene's test**  
    `stattest::test::LevenesTest`  
    *Tests:* homogeneity of variances  

  - **F-test of equality of variances**  
    `stattest::test::FTest`  
    *Tests:* homogeneity of variances  

  - **Shapiro-Wilk test**  
    `stattest::test::ShapiroWilkTest`  
    *Tests:* normality  
