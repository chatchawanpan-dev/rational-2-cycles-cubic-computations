# Rational 2-Cycles for Cubic Polynomials: Computations

This repository contains the SageMath and Magma files used for the manuscript

> Rational 2-Cycles for \(x^3+bx+a\) and the Elliptic Family
> \(Y^2=X^3+4X^2+16t^2\)

The files are included for reproducibility of the exact symbolic, algebraic,
descent, finite-field, and Chabauty computations cited in the manuscript.

## Contents

- `sage_verification/`
  - `verify_table_dictionary.sage`
  - `verify_table_dictionary_output.txt`
  - Generates the corrected sample rows in Table 1 and verifies the cycle
    equations in exact rational arithmetic.
  - `verify_sage_identities.sage`
  - `verify_sage_identities_output.txt`
  - Verifies the auxiliary elliptic curve in the 3-torsion argument, the
    duplication formulas, and the division-polynomial evaluations used for the
    infinite-order theorem.
- `magma_5torsion_verification/`
  - Magma scripts and outputs for the order-5 fiber-product, normalization,
    Prym, torsion, local-bound, and fixed-locus certificates.
- `magma_7torsion_verification/`
  - Magma scripts and outputs for the order-7 algebraic identities, rational
    point search, two-cover descent, and elliptic Chabauty certificates.
- `short_manuscript_sage/`
  - SageMath script and outputs imported from the earlier short manuscript.
    These files support the division-polynomial checks and sample tables used
    in the combined manuscript.
- `software_versions_feedback1_revision.txt`
  - Software and command records from the Feedback 1 revision computation bundle.

## Reproducing the Outputs

Run the SageMath verification:

```bash
sage sage_verification/verify_table_dictionary.sage \
  > sage_verification/verify_table_dictionary_output.txt

sage sage_verification/verify_sage_identities.sage \
  > sage_verification/verify_sage_identities_output.txt
```

Run the Magma order-5 certificates:

```bash
magma magma_5torsion_verification/scripts/c5_x0_normalization_degeneracy_certificate.m \
  > magma_5torsion_verification/results/c5_x0_normalization_degeneracy_certificate_output.txt

magma magma_5torsion_verification/prym_chabauty_route/scripts/c5_ritzenthaler_romagny_prym_certificate.m \
  > magma_5torsion_verification/prym_chabauty_route/results/c5_ritzenthaler_romagny_prym_certificate_output.txt

magma magma_5torsion_verification/prym_chabauty_route/scripts/c5_prym_torsion_fixed_locus_certificate.m \
  > magma_5torsion_verification/prym_chabauty_route/results/c5_prym_torsion_fixed_locus_certificate_output.txt
```

Run the Magma order-7 certificates:

```bash
magma magma_7torsion_verification/scripts/verify_7torsion_basic.m \
  > magma_7torsion_verification/results/verify_7torsion_basic_output.txt

magma magma_7torsion_verification/reviewer_chabauty_route/scripts/c7_reviewers_elliptic_chabauty.m \
  > magma_7torsion_verification/reviewer_chabauty_route/results/c7_reviewers_elliptic_chabauty_output.txt
```

Run the earlier short-manuscript Sage computations:

```bash
sage short_manuscript_sage/cubic_2periodic_supplementary.sage \
  > short_manuscript_sage/cubic_2periodic_supplementary_output.txt
```

## Software

The Feedback 1 table-verification output was archived using SageMath 10.9.
The Magma certificate outputs were archived with Magma V2.29-6. See
`software_versions_feedback1_revision.txt` for the recorded executables and
commands.
