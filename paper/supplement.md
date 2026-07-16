---
title: "Supplementary Material for NWKIT"
---

# Supplementary Methods

The files in `paper/analysis` generate the inventory, independent-output checks, predefined input corpus, worked PEPC example, benchmarks, and figures. Randomized analyses use seed 20260716. Raw TSV and JSON results are retained in `paper/results`; this document contains human-readable summaries.

## Table S1. Software versions and source provenance

| software | version_or_commit | role | provenance | accessed |
| --- | --- | --- | --- | --- |
| NWKIT | 0.27.0; f71dc345ac83 | focal command-line toolkit | https://github.com/kfuku52/nwkit | 2026-07-16 |
| ETE | 4.4.0 | NWKIT tree representation and comparison library | https://github.com/etetoolkit/ete | 2026-07-16 |
| Biopython | 1.87 | NWKIT dependency and independent tree parser | https://biopython.org/ | 2026-07-16 |
| Gotree | 0.5.2; 8ce3ce6eb5b14261b0be7ed2b1b02bcf22dfc305 | command-line comparator | https://github.com/evolbioinfo/gotree | 2026-07-16 |
| PhyKIT | 2.3.0 | command-line comparator | https://github.com/JLSteenwyk/PhyKIT | 2026-07-16 |
| Newick Utilities | 1.6; da121155a977197cab9fbb15953ca1b40b11eb87 | command-line comparator | https://github.com/tjunier/newick_utils | 2026-07-16 |
| DendroPy | 5.0.10 | programming-library comparator | https://github.com/jeetsukumaran/DendroPy | 2026-07-16 |
| TreeSwift | 1.1.45 | programming-library comparator | https://github.com/niemasd/TreeSwift | 2026-07-16 |
| ape | 5.8-1 | programming-library comparator | https://cran.r-project.org/package=ape | 2026-07-16 |

## Table S2. Task-level command-line capability classifications

Native denotes a directly documented CLI workflow; Partial denotes a narrower or distributed implementation; — means that an equivalent was not identified in the reviewed version.

| task | NWKIT_0.27.0 | Gotree_0.5.2 | PhyKIT_2.3.0 | Newick_Utilities_1.6 |
| --- | --- | --- | --- | --- |
| Standard-stream composition | Native | Native | Partial | Native |
| Explicit Newick interpretation and ambiguity rejection | Native | — | — | — |
| Consolidated preflight report | Native | Partial | Partial | Partial |
| Routine tree transformations | Native | Native | Native | Native |
| Rooting strategies | Native | Partial | Partial | Partial |
| Topology-aware annotation transfer | Native | Partial | Native | — |
| Tree-collection summaries | Native | Native | Native | Partial |
| Taxonomy-derived topology | Native | Partial | — | — |
| Trait-aware selection and monophyly | Native | Partial | Partial | — |
| Categorical ancestral-state analysis | Native | Partial | Native | — |
| MCMCtree calibration preparation | Native | — | — | — |
| Tree and trait visualization | Native | Native | Native | Native |

## Table S3. Capability evidence notes and decision rules

| task | software | classification and evidence note | decision rule |
| --- | --- | --- | --- |
| Standard-stream composition | NWKIT 0.27.0 | Native: all 30 functional commands read the primary input from stdin by default; 28 write the primary result to stdout | Native requires a documented tree-producing command to accept stdin and emit its primary result to stdout without a temporary file |
| Standard-stream composition | Gotree 0.5.2 | Native: commands accept stdin and emit stdout | Native requires a documented tree-producing command to accept stdin and emit its primary result to stdout without a temporary file |
| Standard-stream composition | PhyKIT 2.3.0 | Partial: most tree transformations require input file paths and may create derived files; selected reports use stdout | Native requires a documented tree-producing command to accept stdin and emit its primary result to stdout without a temporary file |
| Standard-stream composition | Newick Utilities 1.6 | Native: programs accept stdin and emit stdout | Native requires a documented tree-producing command to accept stdin and emit its primary result to stdout without a temporary file |
| Explicit Newick interpretation and ambiguity rejection | NWKIT 0.27.0 | Native: selectable ETE format, automatic detection, strict automatic mode, and quoted-name control | Count only controls that distinguish alternative interpretations of syntactically readable Newick, not generic parse errors |
| Explicit Newick interpretation and ambiguity rejection | Gotree 0.5.2 | Not identified: Newick is accepted, but no documented equivalent of parser dialect selection or strict ambiguous-numeric-label rejection | Count only controls that distinguish alternative interpretations of syntactically readable Newick, not generic parse errors |
| Explicit Newick interpretation and ambiguity rejection | PhyKIT 2.3.0 | Not identified in reviewed command help | Count only controls that distinguish alternative interpretations of syntactically readable Newick, not generic parse errors |
| Explicit Newick interpretation and ambiguity rejection | Newick Utilities 1.6 | Not identified in the v1.6 manual | Count only controls that distinguish alternative interpretations of syntactically readable Newick, not generic parse errors |
| Consolidated preflight report | NWKIT 0.27.0 | Native: one report covers parse status, duplicate or empty tips, branch-length and support problems, rooting, ultrametricity, binarity, and cross-tree consistency | Native requires one user-facing operation that reports representation, structural, and collection-level checks |
| Consolidated preflight report | Gotree 0.5.2 | Partial: parsing, rootedness statistics, monophyly, and tip-set comparison are separate operations | Native requires one user-facing operation that reports representation, structural, and collection-level checks |
| Consolidated preflight report | PhyKIT 2.3.0 | Partial: task-specific checks and statistics are available, but no equivalent consolidated preflight command was identified | Native requires one user-facing operation that reports representation, structural, and collection-level checks |
| Consolidated preflight report | Newick Utilities 1.6 | Partial: parsing and several consistency checks are distributed across utilities | Native requires one user-facing operation that reports representation, structural, and collection-level checks |
| Routine tree transformations | NWKIT 0.27.0 | Native: prune, rename, collapse, remove or fill annotations, rescale, resolve polytomies, and remove singletons | Grouped row records established overlap and is not used as a novelty claim |
| Routine tree transformations | Gotree 0.5.2 | Native: corresponding commands are documented | Grouped row records established overlap and is not used as a novelty claim |
| Routine tree transformations | PhyKIT 2.3.0 | Native: pruning, renaming, branch collapse, branch-length multiplication, and related operations are documented | Grouped row records established overlap and is not used as a novelty claim |
| Routine tree transformations | Newick Utilities 1.6 | Native: pruning, renaming, rerooting, topology editing, support handling, and related filters are documented | Grouped row records established overlap and is not used as a novelty claim |
| Rooting strategies | NWKIT 0.27.0 | Native: outgroup, transferred root, midpoint, minimal ancestor deviation, and taxonomy-based rooting | Modes count only when exposed in the reviewed command-line interface |
| Rooting strategies | Gotree 0.5.2 | Partial: outgroup and midpoint rooting | Modes count only when exposed in the reviewed command-line interface |
| Rooting strategies | PhyKIT 2.3.0 | Partial: user-specified outgroup rooting; no documented midpoint, MAD, or taxonomy-based mode was identified | Modes count only when exposed in the reviewed command-line interface |
| Rooting strategies | Newick Utilities 1.6 | Partial: explicit outgroup and longest-edge fallback | Modes count only when exposed in the reviewed command-line interface |
| Topology-aware annotation transfer | NWKIT 0.27.0 | Native: names, support, and lengths are mapped by matching clades; root transfer is also supported | Simple tip renaming or positional copying does not qualify |
| Topology-aware annotation transfer | Gotree 0.5.2 | Partial: annotation from a comparison tree and node-name/comment transfer are supported | Simple tip renaming or positional copying does not qualify |
| Topology-aware annotation transfer | PhyKIT 2.3.0 | Native: internal annotations are matched by descendant-taxon bipartitions | Simple tip renaming or positional copying does not qualify |
| Topology-aware annotation transfer | Newick Utilities 1.6 | Not identified in the v1.6 manual | Simple tip renaming or positional copying does not qualify |
| Tree-collection summaries | NWKIT 0.27.0 | Native: clade frequencies, weighted strict/majority/greedy consensus, optional reference tree, and branch-length summaries | Native requires at least consensus or branch-frequency/support calculation across a tree collection |
| Tree-collection summaries | Gotree 0.5.2 | Native: consensus plus classical and transfer bootstrap support | Native requires at least consensus or branch-frequency/support calculation across a tree collection |
| Tree-collection summaries | PhyKIT 2.3.0 | Native: strict or majority consensus and several gene-tree discordance analyses | Native requires at least consensus or branch-frequency/support calculation across a tree collection |
| Tree-collection summaries | Newick Utilities 1.6 | Partial: support calculations are available; no consensus-tree operation was identified | Native requires at least consensus or branch-frequency/support calculation across a tree collection |
| Taxonomy-derived topology | NWKIT 0.27.0 | Native: NCBI, NCBI+APG IV, or user backbones can produce constraints; NCBI/OpenTree/TimeTree can inform rooting | External taxonomy must directly determine a returned topology or root |
| Taxonomy-derived topology | Gotree 0.5.2 | Partial: NCBI taxonomy can be downloaded, pruned, and used to annotate trees, but no integrated species-label-aware constraint workflow was identified | External taxonomy must directly determine a returned topology or root |
| Taxonomy-derived topology | PhyKIT 2.3.0 | Not identified in reviewed command help | External taxonomy must directly determine a returned topology or root |
| Taxonomy-derived topology | Newick Utilities 1.6 | Not identified in the v1.6 manual | External taxonomy must directly determine a returned topology or root |
| Trait-aware selection and monophyly | NWKIT 0.27.0 | Native: grouped monophyly tests, maximum-PD or ranked representative sampling, and clade skimming with trait filters | Selection must choose tips using tree structure plus supplied metadata; monophyly-only implementations are partial |
| Trait-aware selection and monophyly | Gotree 0.5.2 | Partial: monophyly of a supplied tip set; tree sampling is over input trees rather than tips | Selection must choose tips using tree structure plus supplied metadata; monophyly-only implementations are partial |
| Trait-aware selection and monophyly | PhyKIT 2.3.0 | Partial: monophyly tests and phylogenetic-diversity calculations; no equivalent trait-filtered representative-tip selection was identified | Selection must choose tips using tree structure plus supplied metadata; monophyly-only implementations are partial |
| Trait-aware selection and monophyly | Newick Utilities 1.6 | Not identified in the v1.6 manual | Selection must choose tips using tree structure plus supplied metadata; monophyly-only implementations are partial |
| Categorical ancestral-state analysis | NWKIT 0.27.0 | Native: Mk likelihood under ER/SYM/ARD, marginal probabilities, maximum-a-posteriori states, and stochastic maps | Likelihood reconstruction, stochastic mapping, and parsimony are recorded separately in the detailed audit; this grouped row denotes any native categorical-ASR workflow |
| Categorical ancestral-state analysis | Gotree 0.5.2 | Partial: most-parsimonious ancestral characters and sequences | Likelihood reconstruction, stochastic mapping, and parsimony are recorded separately in the detailed audit; this grouped row denotes any native categorical-ASR workflow |
| Categorical ancestral-state analysis | PhyKIT 2.3.0 | Native: discrete model fitting, stochastic character mapping, and parsimony character mapping are documented | Likelihood reconstruction, stochastic mapping, and parsimony are recorded separately in the detailed audit; this grouped row denotes any native categorical-ASR workflow |
| Categorical ancestral-state analysis | Newick Utilities 1.6 | Not identified in the v1.6 manual | Likelihood reconstruction, stochastic mapping, and parsimony are recorded separately in the detailed audit; this grouped row denotes any native categorical-ASR workflow |
| MCMCtree calibration preparation | NWKIT 0.27.0 | Native: user calibrations or TimeTree point/interval estimates can be converted to MCMCtree-labelled Newick | Chronogram plotting or general node annotation does not qualify |
| MCMCtree calibration preparation | Gotree 0.5.2 | Not identified | Chronogram plotting or general node annotation does not qualify |
| MCMCtree calibration preparation | PhyKIT 2.3.0 | Not identified | Chronogram plotting or general node annotation does not qualify |
| MCMCtree calibration preparation | Newick Utilities 1.6 | Not identified | Chronogram plotting or general node annotation does not qualify |
| Tree and trait visualization | NWKIT 0.27.0 | Native: rectangular/radial tree drawing with support and categorical-trait annotation; species-image retrieval is a separate command | Grouped row records overlap and is not used as a novelty claim |
| Tree and trait visualization | Gotree 0.5.2 | Native: text, PNG, SVG, and Cytoscape outputs | Grouped row records overlap and is not used as a novelty claim |
| Tree and trait visualization | PhyKIT 2.3.0 | Native: tree, chronogram, trait-map, network, and other plots | Grouped row records overlap and is not used as a novelty claim |
| Tree and trait visualization | Newick Utilities 1.6 | Native: text, PostScript, and SVG tree displays | Grouped row records overlap and is not used as a novelty claim |

## Table S4. Independent-output checks

| check | n | passed | agreement | max_abs_error | notes |
| --- | --- | --- | --- | --- | --- |
| rooted_RF_vs_DendroPy | 50 | 50 | 1.0 | 0.0 | Rooted symmetric difference on random 16-tip binary trees |
| clade_frequency_vs_direct_count | 20 | 20 | 1.0 | 0.0 | Exact descendant-tip-set counts in collections of 25 random 12-tip trees |
| majority_consensus_vs_direct_count | 20 | 20 | 1.0 | 0.0 | Clades occurring in more than half of 31 rooted trees |
| Mk_marginals_vs_exhaustive_enumeration | 100 | 100 | 1.0 | 2.220446049250313e-16 | Two-state ER model on random four-tip trees; all internal-state assignments enumerated |
| shuffle_seed_reproducibility | 1 | 1 | 1.0 | 0.0 | Exact output equality for repeated seed and inequality for a second seed |

## Table S5. Standard-stream pipeline smoke test

| commands | passed | leaf_names | singleton_nodes | nonroot_branch_lengths | output_newick |
| --- | --- | --- | --- | --- | --- |
| sanitize\|rename\|rescale | True | Taxon_1,Taxon_2,Taxon_3,Taxon_4 | 0 | 2.0,2.0,2.0,2.0,2.0,4.0 | ((Taxon_2:2,Taxon_1:4):2,(Taxon_3:2,Taxon_4:2):2); |

## Table S6. Predefined input corpus

| case | expected_issue | reported_issues | outcome_matched | num_trees |
| --- | --- | --- | --- | --- |
| valid_default | none | none | True | 1 |
| duplicate_tip_label | duplicate_leaf_names | duplicate_leaf_names | True | 1 |
| empty_tip_label | empty_leaf_names | empty_leaf_names | True | 1 |
| negative_branch_length | negative_branch_length | negative_branch_length | True | 1 |
| singleton_node | not_binary | not_binary | True | 1 |
| multifurcation | not_binary | not_binary | True | 1 |
| missing_support | missing_support | missing_support | True | 1 |
| non_ultrametric | not_ultrametric | not_ultrametric | True | 1 |
| unrooted | not_rooted | not_rooted | True | 1 |
| ambiguous_numeric_internal_labels | format_ambiguous | format_ambiguous | True | 1 |
| quoted_node_name | quoted_node_names | quoted_node_names | True | 1 |
| leaf_set_mismatch | leaf_set_mismatch | leaf_set_mismatch | True | 2 |

## Table S7. Benchmark hardware

| key | value |
| --- | --- |
| platform | macOS-26.5.2-x86_64-i386-64bit |
| machine | x86_64 |
| processor | i386 |
| python | 3.12.13 |
| logical_cpus | 12 |
| physical_cpus | 12 |
| total_memory_gb | 64.0 |

## Table S8. Benchmark summaries

Values are medians, minima, and maxima across three process-level runs. Peak RSS includes descendant processes and was sampled every 5 ms.

| benchmark | software | shape | num_tips | num_trees | elapsed_s_median | elapsed_s_min | elapsed_s_max | peak_rss_mb_median | peak_rss_mb_min | peak_rss_mb_max | repetitions |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| consensus_tree_count_scaling | Gotree | mixed_binary | 256 | 10 | 0.05195 | 0.03864 | 0.07325 | 9.707 | 0.1172 | 16.51 | 3 |
| consensus_tree_count_scaling | Gotree | mixed_binary | 256 | 50 | 0.05125 | 0.05098 | 0.05132 | 15.46 | 14.12 | 17.49 | 3 |
| consensus_tree_count_scaling | Gotree | mixed_binary | 256 | 250 | 0.1299 | 0.1242 | 0.252 | 21.35 | 21.21 | 21.73 | 3 |
| consensus_tree_count_scaling | NWKIT | mixed_binary | 256 | 10 | 1.58 | 1.503 | 1.58 | 147.6 | 147.2 | 148.6 | 3 |
| consensus_tree_count_scaling | NWKIT | mixed_binary | 256 | 50 | 1.697 | 1.678 | 1.726 | 150.7 | 150.5 | 152.9 | 3 |
| consensus_tree_count_scaling | NWKIT | mixed_binary | 256 | 250 | 2.604 | 2.589 | 2.633 | 158.1 | 157.9 | 158.2 | 3 |
| consensus_tree_count_scaling | PhyKIT | mixed_binary | 256 | 10 | 0.6443 | 0.6301 | 0.6685 | 39.54 | 39.46 | 39.66 | 3 |
| consensus_tree_count_scaling | PhyKIT | mixed_binary | 256 | 50 | 1.478 | 1.46 | 1.514 | 45.82 | 45.07 | 46.1 | 3 |
| consensus_tree_count_scaling | PhyKIT | mixed_binary | 256 | 250 | 6.29 | 6.244 | 6.31 | 82.42 | 80.8 | 83.12 | 3 |
| validation_tip_scaling | NWKIT | balanced | 128 | 1 | 1.493 | 1.489 | 1.538 | 145 | 144.2 | 145.3 | 3 |
| validation_tip_scaling | NWKIT | balanced | 512 | 1 | 1.503 | 1.479 | 1.536 | 145.7 | 145.1 | 146.6 | 3 |
| validation_tip_scaling | NWKIT | balanced | 2048 | 1 | 1.569 | 1.542 | 1.67 | 150.8 | 150.4 | 151.1 | 3 |
| validation_tip_scaling | NWKIT | balanced | 8192 | 1 | 1.793 | 1.761 | 1.794 | 171.4 | 171.1 | 172.6 | 3 |
| validation_tip_scaling | NWKIT | balanced | 32768 | 1 | 2.694 | 2.58 | 2.777 | 260.8 | 257.8 | 261 | 3 |
| validation_tip_scaling | NWKIT | caterpillar | 128 | 1 | 1.492 | 1.479 | 1.513 | 145.2 | 144.1 | 145.3 | 3 |
| validation_tip_scaling | NWKIT | caterpillar | 512 | 1 | 1.562 | 1.494 | 1.67 | 146.8 | 146.2 | 147 | 3 |
| validation_tip_scaling | NWKIT | caterpillar | 2048 | 1 | 1.568 | 1.56 | 1.586 | 152.4 | 151.9 | 152.4 | 3 |
| validation_tip_scaling | NWKIT | caterpillar | 8192 | 1 | 1.7 | 1.692 | 1.728 | 174.7 | 174.1 | 174.9 | 3 |
| validation_tip_scaling | NWKIT | caterpillar | 32768 | 1 | 2.383 | 2.332 | 2.51 | 266.8 | 266.4 | 268.3 | 3 |

## Table S9. NWKIT command and stream inventory

Inventory for NWKIT 0.27.0 at commit `f71dc345ac83`; 597 tests were collected.

| command | stdin default | stdout default | seed option |
| --- | --- | --- | --- |
| asr | True | True | True |
| cladefreq | True | True | False |
| collapse | True | True | False |
| consensus | True | True | False |
| constrain | True | True | False |
| dist | True | True | False |
| draw | True | False | False |
| drop | True | True | False |
| image | True | False | False |
| info | True | True | False |
| intersection | True | True | False |
| label | True | True | False |
| mark | True | True | False |
| mcmctree | True | True | False |
| monophyly | True | True | False |
| nhx2nwk | True | True | False |
| nwk2table | True | True | False |
| printlabel | True | True | False |
| prune | True | True | False |
| rename | True | True | False |
| rescale | True | True | False |
| root | True | True | False |
| sample | True | True | False |
| sanitize | True | True | False |
| shuffle | True | True | True |
| skim | True | True | True |
| subtree | True | True | False |
| table2nwk | True | True | False |
| transfer | True | True | False |
| validate | True | True | False |
