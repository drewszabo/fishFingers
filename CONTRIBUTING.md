# Contributing to fishFingers

Thank you for your interest in contributing to **fishFingers**.\
We welcome community contributions that improve functionality,
usability, reproducibility, and documentation.

Because fishFingers contains a **published and peer-reviewed machine
learning model**, certain components are governed by stricter
contribution rules (see below).

------------------------------------------------------------------------

## Ways You Can Contribute

We encourage contributions in the following areas:

-   Bug reports and reproducible examples\
-   Feature requests and enhancements\
-   Documentation improvements\
-   Additional tests and validation cases\
-   Performance improvements\
-   Improvements to data preprocessing or fingerprint handling\
-   Cross-platform compatibility

If you are unsure whether your idea fits within scope, please open a
discussion or issue first.

------------------------------------------------------------------------

## Model Governance Policy (Important)

The XGBoost model file (`inst/extdata/fishFingers.json`) is a
**published model** associated with a peer-reviewed scientific
publication.

### Direct modification of the `fishFingers.json` model file is not permitted via pull request.

This includes:

-   Editing model weights\
-   Re-training and overwriting the model\
-   Changing hyperparameters inside the serialized model\
-   Modifying feature ordering inside the saved model

The integrity of this file must remain consistent with the published
version to ensure scientific reproducibility.

------------------------------------------------------------------------

## Proposing Model Updates

If you believe the model should be:

-   Retrained\
-   Extended\
-   Improved\
-   Recalibrated\
-   Updated with additional training data

Please open an issue describing:

1.  The scientific rationale\
2.  The dataset to be used\
3.  Validation approach\
4.  Expected performance impact\
5.  Reproducibility plan

Model updates require:

-   Independent validation\
-   Documentation of training workflow\
-   Version incrementing\
-   Explicit reference to the associated publication\
-   Maintainer approval

Substantive model changes will be handled under a formal versioned
release process.

------------------------------------------------------------------------

## Submitting Changes

1.  Fork the repository\
2.  Create a feature branch\
3.  Make your changes\
4.  Add or update tests where appropriate\
5.  Ensure the package builds cleanly\
6.  Submit a pull request with a clear description

Please reference any related issues in your pull request.

------------------------------------------------------------------------

## Coding Standards

-   Follow tidyverse-style R conventions\
-   Use descriptive variable names\
-   Avoid breaking public APIs without discussion\
-   Include roxygen documentation for exported functions\
-   Keep functions modular and testable

------------------------------------------------------------------------

## Reproducibility and Scientific Integrity

fishFingers is intended for scientific use. Contributions should:

-   Maintain reproducibility\
-   Avoid undocumented algorithmic changes\
-   Clearly state assumptions\
-   Preserve consistency between fingerprints and model input order

------------------------------------------------------------------------

## Questions?

If you are unsure about anything:

-   Open an issue\
-   Start a discussion\
-   Reach out to the maintainers

We appreciate contributions that strengthen the scientific and technical
robustness of fishFingers.
