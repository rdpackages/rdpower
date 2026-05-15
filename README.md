# Power and Sample Size Calculations for RD designs

The package `rdpower` implements power, sample size, and minimum detectable effect calculations for Regression Discontinuity (RD) designs using local polynomial methods.

- `rdpower`: ex-post power calculations for RD treatment effects.
- `rdsampsi`: required sample size calculations for target power.
- `rdmde`: minimum detectable effect calculations.


## Python Implementation

To install/update in Python type:
```
pip install rdpower
```

- Help: [PYPI repository](https://pypi.org/project/rdpower/).

- Replication: [rdpower illustration](Python/rdpower_illustration.py), [senate data](Python/rdpower_senate.csv).

## R Implementation

To install/update in R type:
```
install.packages('rdpower')
```

When the latest `rdrobust` is available on GitHub before CRAN, install it first:
```
install.packages('remotes')
remotes::install_github('rdpackages/rdrobust', subdir = 'R/rdrobust')
```

- Help: [R Manual](https://cran.r-project.org/web/packages/rdpower/rdpower.pdf), [CRAN repository](https://cran.r-project.org/package=rdpower).

- Examples/data: [rdpower illustration](R/rdpower_illustration.R), [senate data](R/rdpower_senate.csv).

## Stata Implementation

To install/update in Stata type:
```
net install rdrobust, from(https://raw.githubusercontent.com/rdpackages/rdrobust/main/stata) replace
net install rdpower, from(https://raw.githubusercontent.com/rdpackages/rdpower/main/stata) replace
```

- Help: [rdpower](stata/rdpower.pdf), [rdsampsi](stata/rdsampsi.pdf), [rdmde](stata/rdmde.pdf).

- Replication: [do-file](stata/rdpower_illustration.do), [senate data](stata/rdpower_senate.dta).


## References

For overviews and introductions, see the [rdpackages website](https://rdpackages.github.io).

### Software and Implementation

- Cattaneo, Titiunik and Vazquez-Bare (2019): [Power Calculations for Regression Discontinuity Designs](https://rdpackages.github.io/references/Cattaneo-Titiunik-VazquezBare_2019_Stata.pdf). _Stata Journal_ 19(1): 210-245.

### Technical and Methodological

- Calonico, Cattaneo and Titiunik (2014): [Robust Nonparametric Confidence Intervals for Regression-Discontinuity Designs](https://rdpackages.github.io/references/Calonico-Cattaneo-Titiunik_2014_ECMA.pdf). _Econometrica_ 82(6): 2295-2326. [Supplemental Appendix](https://rdpackages.github.io/references/Calonico-Cattaneo-Titiunik_2014_ECMA--Supplemental.pdf).
- Calonico, Cattaneo, Farrell and Titiunik (2019): [Regression Discontinuity Designs Using Covariates](https://rdpackages.github.io/references/Calonico-Cattaneo-Farrell-Titiunik_2019_RESTAT.pdf). _Review of Economics and Statistics_ 101(3): 442-451. [Supplemental Appendix](https://rdpackages.github.io/references/Calonico-Cattaneo-Farrell-Titiunik_2019_RESTAT--Supplement.pdf).
- Calonico, Cattaneo and Farrell (2020): [Optimal Bandwidth Choice for Robust Bias Corrected Inference in Regression Discontinuity Designs](https://rdpackages.github.io/references/Calonico-Cattaneo-Farrell_2020_ECTJ.pdf). _Econometrics Journal_ 23(2): 192-210. [Supplemental Appendix](https://rdpackages.github.io/references/Calonico-Cattaneo-Farrell_2020_ECTJ--Supplement.pdf).

## Funding

This work was supported in part by the National Science Foundation through grant [SES-1357561](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1357561).
