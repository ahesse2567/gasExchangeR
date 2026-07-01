# gasExchangeR — Claude Code Session Guide

## About this package

**gasExchangeR** is an R package developed as part of a PhD dissertation. It provides tools for:
- Processing raw gas exchange data (averaging, outlier removal, interpolation)
- Detecting ventilatory thresholds using published algorithms
- Visualizing gas exchange data and threshold results

GitHub: https://github.com/ahesse2567/gasexchanger

The author is not a formally trained software engineer. Please calibrate explanations accordingly — assume familiarity with the science but not necessarily with R packaging conventions or software engineering best practices.

---

## ⚠️ Critical constraint: do not alter algorithm logic

The threshold detection algorithms are intentionally implemented to match specific published scientific papers. The goal of this package is **research reproducibility**, not algorithmic improvement.

- **Do not suggest or make changes to the mathematical logic, formula implementations, or algorithmic steps** inside threshold detection functions, even if a "better" approach exists.
- The key source papers are in the `articles/` folder. Treat these as the ground truth for algorithm behavior.
- If you notice something that *looks* like a bug in an algorithm, flag it as a question ("this differs from what I'd expect — is this intentional?") rather than changing it.
- Data *processing* steps (averaging, interpolation, outlier removal) are fair game for code quality improvements, as long as the behavior is preserved.

---

## Goals for this review

Work through these roughly in order, but flag anything urgent (e.g., `R CMD check` errors) immediately regardless of category.

1. **CRAN readiness** — flag anything that would cause `R CMD check` warnings or errors: DESCRIPTION/NAMESPACE issues, missing documentation, non-standard file structure, etc.
2. **Code quality & readability** — consistent style, clear variable names, reduced repetition, better error messages.
3. **Documentation** — complete and accurate roxygen2 docs for all exported functions; a useful README; vignettes that reflect current behavior.
4. **Tests** — testthat coverage for key functions, especially data processing steps and threshold detection outputs.

---

## How to handle changes

**Small, mechanical changes** (linting, whitespace, naming consistency, roxygen formatting): go ahead and make them, but briefly note what you changed and why.

**Larger changes** (restructuring a function, changing behavior, adding new dependencies, modifying test logic): explain the reasoning and the trade-offs *before* making any edits. Wait for confirmation if the change touches algorithm logic or anything statistically meaningful.

**When in doubt**, ask rather than assume. The author would rather understand a change than just accept it.

---

## Useful commands

```r
# Check the package
devtools::check()

# Load and test interactively
devtools::load_all()

# Build documentation
devtools::document()

# Run tests
devtools::test()
```

---

## Style preferences

- Follow the [tidyverse style guide](https://style.tidyverse.org/) where it doesn't conflict with existing patterns
- Prefer base R or tidyverse; avoid adding heavy new dependencies without discussion
- Keep function interfaces stable — this package may already be cited or used by others
