## Test environments

* local macOS (aarch64-apple-darwin20), R 4.5.2
* win-builder (R-devel and R-release)
* GitHub Actions: ubuntu-latest (R-devel, R-release, R-oldrel-1),
  macOS-latest (R-release), windows-latest (R-release)

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission.

## Notes

* This is the first CRAN submission of this package.
* Some domain terms flagged by the spell checker are intentional and are
  listed in `inst/WORDLIST` (e.g. "recalculate", "multiverse", "Welch",
  "disattenuation").
* The package contains no compiled code and writes only to `tempdir()`.
