# Reproduction

Run from the same Apple Silicon host with the frozen local GeneGalleon image available. The image ID and Python/NumPy/SciPy/ETE versions are recorded in the manifests. `source/nwkit-source.tar.gz` preserves the native Python sources; its 226 source hashes were checked against the running image. The `source/` directory also preserves the benchmark scripts and helpers.

The commands below write new reproduction directories. Published results are in `accuracy/` and `timing/`. All native inference runs use the frozen installed NWKIT in the container, not a host checkout.

```bash
docker run --rm --entrypoint python \
  -e OPENBLAS_NUM_THREADS=1 -e OMP_NUM_THREADS=1 -e MKL_NUM_THREADS=1 \
  -v /Users/kf/repos/nwkit/reviews/ten-shifts-errors-missing-2026-09-11/source:/bench/tools:ro \
  -v /Users/kf/repos/nwkit/reviews/ten-shifts-errors-missing-2026-09-11:/bench/results \
  sha256:c8c69075b5d591499bf3195ec912ef5b5e1aa1cd163f03864639ed6a669e9429 \
  /bench/tools/run_ten_shift_benchmark.py \
  --output /bench/results/reproduction-accuracy --part accuracy \
  --replicates 5 --workers 8 --timeout 1800 --traits 2
```

After all competing benchmark processes finish, run the same command with `--output /bench/results/reproduction-timing --part timing --replicates 1 --workers 1`. This replays both alpha models on replicate 0 of both generating-alpha scenarios. A failure is retained as a failure, not substituted by another seed.

For the likelihood-kernel experiment, use the same Docker flags and image, but run `/bench/tools/benchmark_ten_shift_kernel.py`. It writes `kernel.json`; copy existing kernel results aside before rerunning this experiment. Kernel timings contain three warmups followed by seven batches of ten evaluations for each tested condition.

For profiling, run `/bench/tools/benchmark_ten_shifts_missing.py` with:

```text
--job '{"truth":"shared","mode":"shared","replicate":999,"traits":2,"missing_rate":0.2,"timeout":600,"profile":true}' --output /bench/results/preflight/reproduction-shared.json
```

Use `"mode":"trait-specific"` and a different output path for the other profile. These profiles use a separate seed and are excluded from formal accuracy and isolated timing.

`inputs/*.npz` holds raw observed values (including NaN), known sampling variances, and tip/trait names. Matching JSON files contain complete generating parameters, the true branches, simulation seed and raw-data hash. The worker reconstructs the same 100-tip balanced tree. Its generative data hash is independently checked during input export.

Summaries and the figure can be rebuilt on the host with:

```bash
python source/summarize_ten_shift_benchmark.py .
python source/plot_ten_shift_benchmark.py .
python source/validate_ten_shift_benchmark.py .
```

Run the summary and plotting commands from this report directory after all 20 accuracy outcomes are present. Run validation after the four isolated timing outcomes are also present. Plotting uses Matplotlib; it does not execute native inference. Docker-backed validation does not establish SIF compatibility.
