# Optional external IQ-TREE 3 library worker

NWKIT does not bundle IQ-TREE source, libraries or executables. Its Python
package stays MIT-licensed and contains only its own adapter source and a setup
command. An optional, separately installed worker links to an **unmodified
IQ-TREE 3 library** and keeps one alignment, fixed model and topology loaded
throughout a dating fit. NWKIT sends branch lengths and receives likelihoods
and derivatives over a text protocol.

This optional interface is available in NWKIT 0.43.15 and later. Check that your
installation provides `python -m nwkit.iqtree_library --help`.

## Runtime selection

The IQ-TREE sequence engine requires the ordinary `iqtree3` executable for its
initial unclocked model fit. The library worker is an additional optional program,
not a replacement for that executable in the current implementation.

| `nwkit radte` option | Behavior |
| --- | --- |
| `--sequence-engine iqtree --iqtree-interface auto` | Default IQ-TREE interface: use a valid installed worker when available; otherwise use standard CLI evaluations |
| `--sequence-engine iqtree --iqtree-interface library` | Require the external library worker; fail with setup instructions when absent |
| `--sequence-engine iqtree --iqtree-interface cli` | Use standard IQ2MC CLI exports for every evaluation |
| `--iqtree-worker /absolute/path/nwkit-iqtree-worker` | Select a worker explicitly |

Discovery checks an explicit `--iqtree-worker`, then `NWKIT_IQTREE_WORKER`, then
`nwkit-iqtree-worker` on `PATH`. An explicitly selected missing worker or an
installed incompatible/broken worker is an error. A failed library evaluation
never switches to CLI or to the native sequence engine.

**No download, library build or adapter compilation occurs during an analysis.**
Without the optional worker, other NWKIT features remain available. `auto` can
still use the existing IQ-TREE CLI interface; `library` is unavailable.

## Build the official library once

These instructions use CMake's Makefile/Ninja compile database. Linux and macOS
are supported by the setup script; the complete build and runtime tests were
performed in a Linux ARM64 GeneGalleon Docker environment. A native Windows
worker build is not currently provided. Docker tests do not establish SIF or
native macOS compatibility.

Install a C/C++ toolchain, CMake, Boost and Eigen headers, and zlib development
files. For Debian/Ubuntu, the build dependencies are:

```sh
sudo apt-get install build-essential cmake libboost-dev libeigen3-dev zlib1g-dev
```

For other systems, follow the
[official compilation guide](https://iqtree.github.io/doc/Compilation-Guide).
Clone the official moving default branch and initialize its submodules:

```sh
git clone --recursive https://github.com/iqtree/iqtree3.git "$HOME/src/iqtree3"
cmake -S "$HOME/src/iqtree3" -B "$HOME/src/iqtree3/build-lib" \
  -DBUILD_LIB=ON \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
cmake --build "$HOME/src/iqtree3/build-lib" --parallel 2
```

The relevant output is `build-lib/libiqtree.a`. Keep the matching source tree,
`CMakeCache.txt`, `iqtree_config.h` and `compile_commands.json` in place for the
next step. An ordinary `iqtree3` executable alone is not a linkable library.
Do not mix a library from one build with headers or compiler settings from
another build.

If the ordinary executable is not already installed, it can be built from the
same source into a different build directory:

```sh
cmake -S "$HOME/src/iqtree3" -B "$HOME/src/iqtree3/build-cli" \
  -DCMAKE_BUILD_TYPE=Release
cmake --build "$HOME/src/iqtree3/build-cli" --parallel 2
export PATH="$HOME/src/iqtree3/build-cli:$PATH"
iqtree3 --version
```

No IQ-TREE source patch is needed for either build.

## Build and install the NWKIT adapter once

With NWKIT installed in your Python environment:

```sh
python -m nwkit.iqtree_library build \
  --build-dir "$HOME/src/iqtree3/build-lib" \
  --prefix "$HOME/.local"
export PATH="$HOME/.local/bin:$PATH"
python -m nwkit.iqtree_library check
```

The setup command compiles only NWKIT's adapter and links it with the existing
`libiqtree.a`. It inherits the matching library build's compiler, include paths
and compile definitions. It does not download or rebuild IQ-TREE. It installs
`$HOME/.local/bin/nwkit-iqtree-worker`, outside the NWKIT Python package, after
checking its protocol and version. A failed rebuild preserves the old worker.

For a Conda environment, use `--prefix "$CONDA_PREFIX"` to install the worker in
that environment's `bin` directory. Keep the environment active while building
and running it so its compiler/runtime libraries can be found. Alternatively:

```sh
export NWKIT_IQTREE_WORKER="$HOME/.local/bin/nwkit-iqtree-worker"
python -m nwkit.iqtree_library check
```

The check command prints the IQ-TREE library version, its SHA-256, the adapter
source SHA-256 and worker identity. `check --interface auto` prints `null` when
there is no optional worker. It never compiles anything.

Add the following to an otherwise complete sequence-dating command:

```sh
--sequence-engine iqtree --iqtree-interface library
```

The run manifest records the selected interface, CLI executable identity, library
build identity and worker identity. The worker is started once per likelihood
object; bootstrap replicates have their own refitted models and workers. Closing
or releasing the likelihood object reaps the worker. It is not a machine-wide
service that has to remain running between NWKIT invocations.

## Updating

Updating NWKIT does not download or replace IQ-TREE. To update IQ-TREE, update the
external source checkout and submodules, rebuild the library and the ordinary
executable, then rerun the adapter setup command. Rebuild the adapter after a
change to its source or an incompatible protocol update. The version/hash output
lets an analysis identify the actual installed build.

The adapter calls IQ-TREE C++ interfaces; their compatibility is checked by the
build and numerical tests, not guaranteed across all future upstream changes.
Do not patch IQ-TREE or silently change likelihood models to conceal an upstream
failure. Validate the new build before using it for production analyses.

## Conda packaging without bundling IQ-TREE in NWKIT

Conda can install distinct packages together while keeping their file ownership
and license metadata separate. A proposed layout is:

| Package | Contents and license |
| --- | --- |
| `nwkit` | Existing noarch Python package and original MIT notice; no IQ-TREE binaries/libraries |
| `iqtree` | Ordinary IQ-TREE 3 CLI, GPL-2.0-or-later |
| IQ-TREE library/worker package | Official library build and/or independently installed worker; linked worker distributed under GPL conditions |
| Optional `nwkit-iqtree` metapackage | Depends on the above, providing a one-install experience |

The last two package names/layouts are proposals, **not currently published
installation commands**. Alternatively, the Conda `nwkit` recipe can declare the
external worker as a runtime dependency once it is published for its supported
platforms. The PyPI NWKIT package can remain pure Python and dependency-optional.

The current [Bioconda IQ-TREE recipe](https://github.com/bioconda/bioconda-recipes/tree/master/recipes/iqtree)
builds the normal executable, not this library worker. Simply adding `iqtree`
as a dependency does not provide persistent likelihood evaluation. A library
build/worker recipe or an upstream packaging change is required first.

Declaring a Conda dependency does not itself relicense NWKIT's original MIT
source. It also does not exempt the worker from GPL requirements: a linked
binary must retain notices and be accompanied by the appropriate corresponding
source and build materials. Process/package separation alone is not an automatic
legal exemption; the degree of integration matters. Before publishing a combined
distribution, assess the actual artifacts and dependencies under the
[IQ-TREE license](https://github.com/iqtree/iqtree3/blob/master/LICENSE) and the
[GNU guidance on aggregation and linking](https://www.gnu.org/licenses/gpl-faq.en.html#MereAggregation).
Do not label an IQ-TREE-containing worker binary as MIT-only.
