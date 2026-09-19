# Native joint shift inference: adoption protocol

Status: research implementation; **not approved as a GeneGalleon default**.
Protocol recorded before collecting the new procedure's calibration outcomes.
This study does not inherit validation of the older single-trait calibrated
backend. Continuous covariance optimization, ordinary ML, missing observations,
and adaptive screening define a different procedure.

The candidate families are nested by `locations + free regime offsets`.
For each family, compare its best likelihood with the unrestricted candidate
family. Simulate from the fitted restricted model and repeat the complete
configured search, including screening and covariance estimation. Test families
in increasing order and stop at the first non-rejection. Monte Carlo p-values
use `(1 + exceedances)/(B + 1)`, including ties. This is a plug-in parametric
bootstrap, with no claim of uniform composite-null control or global continuous
optimization. Failed fits terminate the analysis; never discard failed draws.

Freeze seeds, search budgets, alpha range, optimizer settings and simulation
truths in an executable study specification before running a study. Changing any
of those after inspecting results requires a new independent confirmation set.
Known-error inputs and observation masks remain fixed in every bootstrap draw.
Report finite-alpha and Brownian/independent boundary truths separately, with
fixed and stationary roots, 1/4 traits, zero/heterogeneous known errors, and
complete/20% missing observations. Include shared and distinct shifted regimes,
nested returns, weak/strong effects, and balanced/pectinate trees. The stationary
root has no Brownian zero-alpha truth. Additional estimated error requires its
own strata. Record all failures and actual null/power denominators.

Adoption gates, all required:

1. Independent dense Gaussian likelihood/GLS checks agree within numerical
   tolerance across covariance boundaries and layouts. Small-tree exhaustive
   enumeration contains the existing candidate space and generalized layouts.
2. At nominal 0.05, a one-sided 95% binomial upper confidence bound on false
   selection is at most 0.075 in **each** predeclared null stratum, using at least
   1,000 independent datasets per stratum and 999 bootstrap draws per dataset.
   This is an operational gate, not a proof of uniform error control.
3. Under comparable models on the same simulated datasets, the one-sided 95%
   paired-bootstrap lower bound for detection-power difference against the
   existing production method exceeds -0.05 in each predeclared alternative
   stratum. Count correct branch and correct grouping recovery separately;
   require their corresponding lower differences also to exceed -0.05. Use at
   least 1,000 datasets and publish the baseline invocation and output mapping.
4. Report timings, peak memory, convergence/failure rates and heuristic versus
   exhaustive likelihood gaps. Measure balanced and pectinate 32/128/512/1,000-tip
   trees with 1/4 traits and up to 10 shifts. Performance claims must separate
   fixed-layout GLS, search, calibration and support costs. Never imply a
   whole-workflow speedup from the fixed-layout kernel alone.
5. Verify GeneGalleon replicate aggregation, missingness, original units, JSON
   provenance/restart rejection, TSV consumers and plotting in its container.
   Remove kfl1ou and switch defaults only after the statistical gates pass.

Bootstrap selection frequencies describe stability under the selected fitted
model, not posterior probabilities or hypothesis-test p-values. Alpha limits
remain first-class results; unsupported optimum parameters are NA. A failed or
incomplete gate keeps the new method research-only and the production dependency
in place. Do not narrow the study or relax thresholds to obtain acceptance.
