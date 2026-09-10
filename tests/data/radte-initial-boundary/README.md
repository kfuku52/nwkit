# Collapsed branch-only initializer

External AliSim input from default-profile validation, internal family 99:
four species, six gene tips, 300 codons, independent log-rate SD 0.3,
simulation seed 830099, AliSim seed 10830099, optimizer seed 20830099.
Generator: `GY{0.5,2}+FQ+G4{0.7}` (omega 0.5, kappa 2, gamma shape 0.7).
Species ages are fixed. The true duplication age is 7.5.

Before the repair, the branch-only initializer puts duplication D at almost
the parent age 10. The exact sequence start has an age gradient about 2.96e9;
all three SLSQP starts report incompatible inequalities despite feasible ages.
The regression reproduces that failure on the original implementation.

Resetting only this near-minimum-duration sequence initializer to the existing
chronology interior recovers a feasible exact objective minimum (age about
9.3410, objective about 2038.58376). Independent restarts of the same objective
agree. This repairs numerical initialization; it does not assert that the
estimated age equals truth or calibrate the profile confidence interval.

The source-study failure remains counted in its original 100-family denominator.
These input files are simulated data, not empirical gene sequences.
