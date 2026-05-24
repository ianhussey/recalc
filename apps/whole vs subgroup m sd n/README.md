# ANCHOR: numerical consistency between whole sample and subgroups

A Shiny front-end to the `{recalc}` package's E1 / E2 checks on the
relationship between a paper's whole-sample statistics and its subgroup
breakdown. Two tabs:

## Tab 1 — Whole vs subgroups

Compares the reported overall **N**, **mean**, and **SD** against the values
recalculated from any number of subgroups. Each row reports the rounding-aware
recalculated interval and a `consistent` verdict (`TRUE` iff the reported and
recalculated intervals overlap). Backed by `recalc::recalc_total_from_subgroups()`.

Identities used:

- **N** = Σ n_g
- **M** = Σ n_g M_g / N
- **SD²** = [ Σ (n_g − 1) s_g² + Σ n_g (M_g − M)² ] / (N − 1)

The SD identity is the **total-variance decomposition** (within + between
sums of squares). Earlier versions of this app used the textbook *pooled SD*
formula √( Σ (n_g − 1) s_g² / Σ (n_g − 1) ), which is equal to the total SD
**only when subgroup means coincide**; using it as a check would false-flag
any paper whose subgroup means differ. The current implementation uses the
correct identity via the `{recalc}` package.

## Tab 2 — Implied missing subgroup

When a paper reports overall statistics and only *k − 1* of *k* subgroups,
this tab inverts the aggregation identities to recover what the unreported
subgroup would have to be:

- *n*★ = N − Σ n_g
- *M*★ = (N M − Σ n_g M_g) / n★
- *s*★² = [(N − 1) SD² − Σ (n_g − 1) s_g² − n★ (M★ − M)² − Σ n_g (M_g − M)²] / (n★ − 1)

Optionally accepts **scale endpoints** (`scale_min`, `scale_max`) for the
measurement (e.g. 1 and 7 for a 7-point Likert). When supplied:

- the implied mean is checked against the scale,
- the implied SD is checked against the Bhatia-Davis upper bound
  √(n★ / (n★ − 1) · (M★ − a)(b − M★)) — the largest sample SD any
  n★-sample on [a, b] with mean M★ can attain (Bhatia & Davis, 2000;
  Bessel-corrected).

The SD row also fails when the implied variance is strictly negative at any
rounding corner — *no real subgroup* can reconcile the reported values.

Backed by `recalc::recalc_missing_subgroup()`.

## Cite as

Hussey, I. (2024). Assessing Numerical Consistency between wHOle sample and
subgRoups (ANCHOR): A method to check for inconsistencies in reported summary
statistics for the full sample vs. subgroups. https://github.com/ianhussey/recalc

This method of error-checking is well-established; I make no claim to having
invented it. I came across it in Wilkinson et al.'s draft of INSPECT-SR
(see [here](https://doi.org/10.1136/bmjopen-2024-084164) for the development
protocol). The contribution here is an accessible web app exposing it.

## References

Bhatia, R., & Davis, C. (2000). A Better Bound on the Variance.
*The American Mathematical Monthly*, **107**(4), 353-357.
