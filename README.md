# RDPOWER

The ``rdpower`` package provides Stata and R implementations of power and sample size calculations using robust bias-corrected local polynomial inference methods.

This work was supported by the National Science Foundation through grant [SES-1357561](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1357561).

## Website

https://rdpackages.github.io/rdpower

## Stata Implementation

To install/update in Stata type:
```
net install rdpower, from(https://raw.githubusercontent.com/rdpackages/rdpower/master/stata) replace
```

- Help: [rdpower](stata/rdpower.pdf), [rdsampsi](stata/rdsampsi.pdf).

- Replication: [do-file](stata/rdpower_illustration.do), [data-senate](stata/rdpower_senate.dta).

## R Implementation

To install/update in R type:
```
install.packages('rdpower')
```
- Help: [R Manual](https://cran.r-project.org/web/packages/rdpower/rdpower.pdf), [CRAN repository](https://cran.r-project.org/package=rdpower).

- Replication files: [R-script](R/rdpower_illustration.R), [data-senate](R/rdpower_senate.csv).

## References

For overviews and introductions, see [rdpackages website](https://rdpackages.github.io).

### Software and Implementation

- Cattaneo, Titiunik and Vazquez-Bare (2019): [Power Calculations for Regression Discontinuity Designs](https://rdpackages.github.io/references/Cattaneo-Titiunik-VazquezBare_2019_Stata.pdf), _Stata Journal_ 19(1): 210-245.

### Technical and Methodological

- Calonico, Cattaneo and Titiunik (2014): [Robust Nonparametric Confidence Intervals for Regression-Discontinuity Designs](https://rdpackages.github.io/references/Calonico-Cattaneo-Titiunik_2014_ECMA.pdf), _Econometrica_ 82(6): 2295-2326. [Supplemental Appendix](https://rdpackages.github.io/references/Calonico-Cattaneo-Titiunik_2014_ECMA--Supplemental.pdf).

- Calonico, Cattaneo, Farrell and Titiunik (2019): [Regression Discontinuity Designs Using Covariates](https://rdpackages.github.io/references/Calonico-Cattaneo-Farrell-Titiunik_2019_RESTAT.pdf), _Review of Economics and Statistics_ 101(3): 442-451. [Supplemental Appendix](https://rdpackages.github.io/references/Calonico-Cattaneo-Farrell-Titiunik_2019_RESTAT--Supplement.pdf).

- Calonico, Cattaneo and Farrell (2020): [Optimal Bandwidth Choice for Robust Bias Corrected Inference in Regression Discontinuity Designs](https://rdpackages.github.io/references/Calonico-Cattaneo-Farrell_2020_ECTJ.pdf), _Econometrics Journal_ 23(2): 192-210. [Supplemental Appendix](https://rdpackages.github.io/references/Calonico-Cattaneo-Farrell_2020_ECTJ--Supplement.pdf).

<br><br>
