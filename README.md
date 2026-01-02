# 🐟 fishFingers

fishFingers is a non-targeted screening prioritisation tool for high-resolution mass spectrometry data. It allows for the prediction of bioconcentration factors (BCF) in over 100 species of fish using structural fingerprints derived from tandem mass spectrometry data (MS2).

## Installation
Install the latest version from GitHub:
```
# install.packages("remotes")
remotes::install_github("drewszabo/fishFingers")
```

## Package overview

fishFingers provides:
* Automated generation of structural fingerprints from SMILES
* Import predicted fingerprints from SIRIUS CSI:FingerID
* A pre-trained machine-learning model for BCF prediction
* A simple interface for single or multiple chemicals
* Internal reference data stored in inst/extdata
* Development scripts and results are stored in dev/

## Quick start

```
library(fishFingers)

# Example SMILES string
smiles <- c(
  "CCOC(=O)C1=CC=CC=C1",
  "CCN(CC)CCOC(=O)C1=CC=CC=C1"
)

# Predict logBCF
pred <- predict_bcf(smiles, "smiles", "Cyprinus carpio")

pred
```

## Input requirements
* SMILES should be valid and parseable by rcdk (use `webchem::is.smiles()` to check your strings beforehand)
* Predictions are most reliable for chemicals withing the applicability domain from the training data
* Only species in the training and test sets are avialable for prediction at this time (use `fishFingers::check_species(list = TRUE)` for complete list)

## Model details
* Target variable: log10 bioconcentration factor (BCF)
* Descriptors: binary structural fingerprints & one-hot encoded species variables
* Model type: XGBoost
* Training data size: >1000 chemicals

## Contributing

Contributions are welcome, including:
* Bug reports
* Feature requests
* Model improvements
* Documentation enhancements

Please open an issue or submit a pull request on GitHub.
