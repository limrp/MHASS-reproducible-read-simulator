# Milestone 1 — MHASS reproducibility audit

## Scope

This audit evaluates the unmodified MHASS source at commit:

`93fbce910fef5088b923e18fc5681e277180f632`

The purpose of Milestone 1 was to determine whether repeated MHASS
executions could produce bitwise-identical outputs when the execution
environment, input files, parameters, and thread count were held
constant.

No functional modifications to the MHASS source code were introduced
during this milestone.

The target reproducibility criterion for the thesis benchmark is:

> Under the same frozen execution environment, identical inputs,
> parameters, and seeds should produce bitwise-identical benchmark
> outputs that can be verified using SHA-256 checksums.


## Audit strategy

A small test dataset was derived from the upstream ATCC 16S test data
to allow rapid and interpretable reproducibility testing.

The audit dataset contained:

- 3 ASVs
- 3 distinct genomes
- 1 sample
- 30 requested reads
- uniform genome abundance distribution
- lognormal `np` distribution
- 1 execution thread

The selected genomes were:

- `Bifidobacterium_adolescentis`
- `Enterococcus_faecalis`
- `Clostridium_beijerinckii`

Input SHA-256 fingerprints were:

| Input | SHA-256 |
|---|---|
| `amplicons_3.fa` | `3d0cf69e97513e7d6cfc36a689023e1e2855dee5932bc0ce26f0d39e06e8906a` |
| `amplicon_labels_3.tsv` | `8ff11b1427b41ae89a7cf965335d6609e399e2bbc1d3f8bc21da8864d6eb3e35` |
| upstream `barcodes.tsv` | `ab1d64d22be7fb827b68c4d6621b2c27cffd4b1d6fe7d8daccdcdbbe98c3e032` |


# Findings

## Finding 1 — The upstream dependency specification allowed an incompatible PBSIM3 version

The original installation resolved:

`pbsim3 3.0.4`

At the audited MHASS commit, MHASS expects PBSIM to generate:

`sim.bam`

However, PBSIM3 3.0.4 generated:

`sim.sam`

MHASS therefore attempted to pass a nonexistent `sim.bam` file to
PacBio CCS.

The original PBSIM executable had SHA-256:

`047c979c5b0b766cc109803430454f47003e8c69499453d2b13408c2a3f5a40c`

After upgrading the environment to:

`pbsim3 3.0.5 h9948957_2`

the Conda PBSIM executable and the copy used internally by MHASS both
had SHA-256:

`c58e20aaad3a39dc27a02f5840dc6eca16531e796ff79cb829c6499cedaecc30`

With PBSIM3 3.0.5, the unmodified MHASS smoke test completed
successfully.

### Potential solution

The reproducible installation should pin the tested PBSIM dependency,
including its build where appropriate:

`pbsim3=3.0.5=h9948957_2`

The installer should also verify that the PBSIM executable copied into
the MHASS resource directory is bitwise identical to the expected
Conda executable.


## Finding 2 — The installer does not locate the tested QSHMM model path

The tested PBSIM3 installation stored the model at:

`$CONDA_PREFIX/data/QSHMM-RSII.model`

with SHA-256:

`cff55ec94274a919529266561df48c3f72506e936a9c0ceb7840240ecdedb3b7`

The upstream installation script did not search this location.
Nevertheless, MHASS successfully found the model through its runtime
fallback logic.

### Potential solution

The reproducible installer should additionally search:

```bash
"$CONDA_ENV_PATH/data/QSHMM-RSII.model"
````

along with the existing candidate locations.

The selected QSHMM model should be recorded by checksum in the
reproducibility manifest.

## Finding 3 — Native metaSPARSim count generation is not reproducible

MHASS's `get_counts.R` contains:

```r
set.seed(42)
```

before calling metaSPARSim.

Despite this, two executions of `get_counts.R` using identical
inputs and parameters generated different count matrices.

One isolated experiment produced:

| Taxon           | Run 1 | Run 2 |
| --------------- | ----: | ----: |
| Bifidobacterium |    13 |     8 |
| Enterococcus    |     8 |     7 |
| Clostridium     |     9 |    15 |

The files also had different SHA-256 hashes.

The installed dependency was:

* metaSPARSim `1.1.2`
* GitLab commit `37b5931ab1ed6836ded1b605ab23e41824d0f91e`

Inspection of that exact source commit showed that the
`random_unif()` C++ helper initializes its random-number generator
with:

```cpp
std::default_random_engine rng(time(0));
```

Therefore this component is seeded from the current wall-clock time
rather than from R's `.Random.seed`.

This behavior was confirmed experimentally. Resetting
`set.seed(42)`, waiting two seconds, and calling `random_unif()` again
produced different values. Repeating the complete `metaSPARSim()`
simulation with the same R seed and a two-second delay also produced
different count matrices.

### Potential solutions

Possible approaches include:

1. modifying metaSPARSim's native RNG;
2. selecting another metaSPARSim version with explicit deterministic
   seed control, provided its scientific behavior remains compatible;
3. bypassing metaSPARSim for the reproducible benchmark by supplying
   an externally defined count matrix.

For this thesis, option 3 is preferred.

## Finding 4 — Lognormal `np` sampling is not seeded

MHASS samples the number of passes using:

```python
np.random.lognormal(...)
```

without initializing a NumPy seed.

The two audit runs generated different `np` assignments in
`sequence_file_mapping.tsv`.

Different `np` values modify the number of simulated passes supplied
to PBSIM and CCS. Consequently, this variation can change sequencing
error profiles, predicted CCS accuracy, and the number of simulated
molecules that pass the HiFi quality threshold.

### Potential solution

Introduce a user-facing master `--seed` and derive a dedicated,
deterministic seed for the NumPy `np` sampler.

## Finding 5 — Empirical `np` sampling is not seeded

The empirical `np` mode uses Python:

```python
random.choice(...)
```

without a controlled seed.

### Potential solution

Use an explicit Python random-number generator initialized from a
component seed derived from the user-facing master `--seed`.

## Finding 6 — PBSIM randomness is not controlled by MHASS

PBSIM3 provides a seed option, but MHASS does not pass one.

Therefore PBSIM's simulated sequencing errors are not currently tied
to any user-facing MHASS seed.

### Potential solution

Derive a deterministic PBSIM seed for every template from the master
seed and a stable template identifier.

For example, conceptually:

```text
master seed
    |
    +-- SHA256(master_seed | "pbsim" | stable_template_id)
            |
            +-- PBSIM template seed
```

The derived seed must not depend on thread scheduling or task
completion order.

## Finding 7 — Final read orientation is random and uncontrolled

During FASTQ combination, MHASS uses:

```python
random.random() < 0.5
```

to determine whether a read should be reverse-complemented.

The existing `random.seed(42)` statement in the source is commented
out.

### Potential solution

Use a dedicated deterministic orientation RNG derived from the master
seed.

## Finding 8 — Some stages were already deterministic

Across the two audit runs:

* `counts_params.tsv` was identical;
* `counts_meta.tsv` was identical;
* `sample_barcode_map.tsv` was identical.

This indicates that the parameter calculations and barcode assignment
were deterministic in the tested configuration.

The audit therefore localizes the reproducibility failures to specific
stochastic components rather than to the entire program.


# Reproducibility test result

Two complete MHASS executions with identical inputs and parameters
both finished successfully but generated:

* 24 final CCS reads in run 1;
* 23 final CCS reads in run 2.

The following retained outputs differed:

| Output                      | Result    |
| --------------------------- | --------- |
| `counts.tsv`                | DIFFERENT |
| `counts_meta.tsv`           | SAME      |
| `counts_params.tsv`         | SAME      |
| `sample_barcode_map.tsv`    | SAME      |
| `sequence_file_mapping.tsv` | DIFFERENT |
| `combined_reads.fastq`      | DIFFERENT |

Therefore the audited upstream MHASS configuration does not satisfy
the thesis requirement for bitwise reproducibility.

# Blocker summary

| Blocker                         | Cause                                       | Proposed action                                  |
| ------------------------------- | ------------------------------------------- | ------------------------------------------------ |
| PBSIM3 interface mismatch       | unpinned dependency allowed PBSIM3 3.0.4    | pin tested PBSIM3 3.0.5 build                    |
| QSHMM path mismatch             | installer does not search tested Conda path | update installer lookup and checksum model       |
| Count simulation nondeterminism | metaSPARSim C++ RNG uses `time(0)`          | bypass with fixed count matrix in benchmark mode |
| Lognormal `np` nondeterminism   | unseeded NumPy RNG                          | control through master seed                      |
| Empirical `np` nondeterminism   | unseeded Python RNG                         | control through master seed                      |
| PBSIM nondeterminism            | no PBSIM seed passed                        | deterministic per-template seeds                 |
| Orientation nondeterminism      | unseeded `random.random()`                  | deterministic orientation RNG                    |
| Environment drift               | broad dependency specifications             | pin/lock environment and later containerize      |

# Milestone 1 conclusion

Milestone 1 successfully demonstrated that the audited upstream MHASS
configuration is not bitwise reproducible.

More importantly, the audit identified the principal causes rather
than simply observing that the final FASTQ files differed.

The recommended strategy is not to redesign metaSPARSim. Instead, the
reproducible benchmark configuration should accept a fixed external
count matrix, while the remaining sequencing-related random processes
should be controlled through a master `--seed`.

This allows native MHASS behavior to remain available while providing
a separate, controlled configuration suitable for the thesis
benchmark.

# The key scientific point

chatGPT particularly recommend keeping the **“intended counts versus realized HiFi reads”** distinction. It will prevent a methodological problem later.

My ground truth will really have layers:

```text
Fixed biological/community truth
counts.tsv
        ↓
intended templates / molecules
        ↓
PBSIM sequencing simulation
        ↓
CCS quality filtering
        ↓
realized simulated HiFi reads
        ↓
DADA2 / VSEARCH / UNOISE3
        ↓
inferred ASVs / OTUs
````

That gives me two useful comparison points:

```text
counts.tsv
→ What community did I intend to simulate?

read-level manifest
→ What reads did the simulator actually deliver?

pipeline result
→ What did DADA2/VSEARCH/UNOISE3 recover?
```

For a methodological benchmarking thesis, that is stronger than simply saying *“I generated random simulated reads and used them as ground truth.”* 
I will know precisely what truth exists at each stage.
