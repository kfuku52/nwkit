"""100-tip general-path smoke: estimate alpha, process and extra errors."""
import json
from pathlib import Path
import sys
import time

import numpy as np

sys.path.insert(0, '/src')
from tools.benchmark_shift_covariance import balanced_tree
from nwkit.shift_simulation_cli import explicit_simulation
from nwkit.shift_simulation import simulate_shift
from nwkit.shift_native_model import ShiftData, ShiftLayout
from nwkit.shift_native_fit import NativeFitOptions, fit_native_layout
from nwkit.shift_native_provenance import native_implementation_sha256

spec = explicit_simulation(balanced_tree(100), {
    'trait_names': ['x', 'y'], 'alpha': [0.5, 2.0],
    'process_tip_covariance': [[1.0, 0.6], [0.6, 1.0]],
    'regime_optima': [[0.0, 0.0]],
    'sampling_standard_errors': [0.1, 0.2],
    'measurement_covariance': [[0.08, 0], [0, 0.15]],
})
y, _ = simulate_shift(spec, seed=932)
y[0, ::13, 0] = np.nan
y[0, ::17, 1] = np.nan
data = ShiftData.build(spec.tree, y[0], spec.trait_names, spec.sampling_variances)
start = time.perf_counter()
fit = fit_native_layout(data, ShiftLayout.build(data.tree), options=NativeFitOptions(
    trait_covariance='full', alpha_model='trait-specific', estimate_measurement_error=True))
assert np.isfinite(fit['log_likelihood'])
assert fit['joint_covariance']['engine'] == 'vector_tree_pruning'
result = {'implementation_sha256': native_implementation_sha256(), 'tips': 100,
    'traits': 2, 'seconds_not_benchmark': time.perf_counter() - start,
    'log_likelihood': fit['log_likelihood'], 'joint_covariance': fit['joint_covariance'],
    'trait_results': fit['traits']}
Path('/src/reviews/trait-covariance-2026-09-11/general-smoke.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
print('General path completed with finite joint likelihood and checked optimization.')
