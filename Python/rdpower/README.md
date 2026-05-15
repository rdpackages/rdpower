# RDPOWER

The package `rdpower` implements power, sample size, and minimum detectable effect calculations for Regression Discontinuity (RD) designs using local polynomial methods.

- `rdpower`: ex-post power calculations for RD treatment effects.
- `rdsampsi`: required sample size calculations for target power.
- `rdmde`: minimum detectable effect calculations.

See references for methodological and practical details.

Website: [https://rdpackages.github.io/](https://rdpackages.github.io/).

Source code: [https://github.com/rdpackages/rdpower](https://github.com/rdpackages/rdpower).

## Authors

Matias D. Cattaneo (<matias.d.cattaneo@gmail.com>)

Ricardo Masini (<ricardo.masini@gmail.com>)

Rocio Titiunik (<rocio.titiunik@gmail.com>)

Gonzalo Vazquez-Bare (<gvazquezbare@gmail.com>)


## Installation

To install/update use pip:
```
pip install rdpower
```

## Usage
```python
import pandas as pd

from rdpower import rdmde, rdpower, rdsampsi

data = pd.read_csv("rdpower_senate.csv")
z = data[["demvoteshfor2", "demmv"]]

power = rdpower(data=z, tau=5)
sample_size = rdsampsi(data=z, tau=5)
mde = rdmde(data=z)
```

- Replication: [rdpower illustration](https://github.com/rdpackages/rdpower/blob/main/Python/rdpower_illustration.py), [senate data](https://github.com/rdpackages/rdpower/blob/main/Python/rdpower_senate.csv).


## Dependencies

- numpy
- pandas
- scipy
- rdrobust
- matplotlib

## References

For overviews and introductions, see [rdpackages website](https://rdpackages.github.io).

### Software and Implementation

- Cattaneo, Titiunik and Vazquez-Bare (2019): [Power Calculations for Regression Discontinuity Designs](https://rdpackages.github.io/references/Cattaneo-Titiunik-VazquezBare_2019_Stata.pdf).<br>
_Stata Journal_ 19(1): 210-245.

### Technical and Methodological

- Calonico, Cattaneo and Titiunik (2014): [Robust Nonparametric Confidence Intervals for Regression-Discontinuity Designs](https://rdpackages.github.io/references/Calonico-Cattaneo-Titiunik_2014_ECMA.pdf).<br>
_Econometrica_ 82(6): 2295-2326.<br>
[Supplemental Appendix](https://rdpackages.github.io/references/Calonico-Cattaneo-Titiunik_2014_ECMA--Supplemental.pdf).

- Calonico, Cattaneo, Farrell and Titiunik (2019): [Regression Discontinuity Designs Using Covariates](https://rdpackages.github.io/references/Calonico-Cattaneo-Farrell-Titiunik_2019_RESTAT.pdf).<br>
_Review of Economics and Statistics_ 101(3): 442-451.<br>
[Supplemental Appendix](https://rdpackages.github.io/references/Calonico-Cattaneo-Farrell-Titiunik_2019_RESTAT--Supplement.pdf).

- Calonico, Cattaneo and Farrell (2020): [Optimal Bandwidth Choice for Robust Bias Corrected Inference in Regression Discontinuity Designs](https://rdpackages.github.io/references/Calonico-Cattaneo-Farrell_2020_ECTJ.pdf).<br>
_Econometrics Journal_ 23(2): 192-210.<br>
[Supplemental Appendix](https://rdpackages.github.io/references/Calonico-Cattaneo-Farrell_2020_ECTJ--Supplement.pdf).


<br><br>
