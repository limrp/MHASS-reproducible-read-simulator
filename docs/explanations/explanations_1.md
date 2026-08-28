Absolutely. Think of MHASS as a **factory for synthetic sequencing reads**.

You give it the DNA sequences that could exist in your mock community, tell it which organism each sequence belongs to, tell it how many samples you want, and give it barcodes. MHASS then decides how many copies of each sequence should exist, simulates PacBio subreads, runs CCS, and finally gives you simulated HiFi-like FASTQ reads.

For your thesis, the key improvement is very simple:

> Instead of letting MHASS randomly decide the abundance of the ASVs every time, **you give it that abundance table yourself**.

That is what the fixed `counts.tsv` does.

---

# 1. The big picture

Imagine that we want to simulate a tiny microbial community with three ASVs:

```text
ASV_A = Bifidobacterium
ASV_B = Enterococcus
ASV_C = Clostridium
```

MHASS currently behaves roughly like this:

```text
amplicon sequences
        +
ASV → genome mapping
        +
number of samples / reads
        |
        v
┌────────────────────────────┐
│ 1. metaSPARSim             │
│                            │
│ Decide HOW MANY copies     │
│ of each ASV exist          │
└─────────────┬──────────────┘
              |
              v
         counts.tsv
              |
              v
┌────────────────────────────┐
│ 2. Template creation       │
│                            │
│ For every simulated copy:  │
│ - choose number of passes  │
│ - attach barcode           │
└─────────────┬──────────────┘
              |
              v
       template FASTAs
              |
              v
┌────────────────────────────┐
│ 3. PBSIM3                  │
│                            │
│ Simulate PacBio subreads   │
│ and sequencing errors      │
└─────────────┬──────────────┘
              |
              v
          sim.bam
              |
              v
┌────────────────────────────┐
│ 4. PacBio CCS              │
│                            │
│ Combine subreads into      │
│ HiFi-like consensus reads  │
└─────────────┬──────────────┘
              |
              v
         ccs.fastq
              |
              v
┌────────────────────────────┐
│ 5. FASTQ combiner          │
│                            │
│ Combine reads and randomly │
│ choose orientation         │
└─────────────┬──────────────┘
              |
              v
     combined_reads.fastq
```

The **fixed-counts modification changes only box 1**.

Everything after `counts.tsv` can stay essentially the same.

---

# 2. Input 1: the amplicon FASTA

This tells MHASS:

> “These are the actual DNA sequences that are allowed to exist in my simulated community.”

A simplified toy FASTA could look like:

```fasta
>ASV_A
AGAGTTTGATCCTGGCTCAG...ACGGCTACCTTGTTACGACTT

>ASV_B
AGAGTTTGATCCTGGCTCAG...ACGGCTACCTTGTTACGACTT

>ASV_C
AGAGTTTGATCCTGGCTCAG...ACGGCTACCTTGTTACGACTT
```

In your real test, the header is much longer, for example something like:

```text
>NC_008618.1 Bifidobacterium adolescentis ... coordinates=1553664-1555155...
```

The header is the **identifier** of that sequence.

You can think of the FASTA as a box of Lego pieces:

```text
FASTA

ASV_A ───── DNA sequence A
ASV_B ───── DNA sequence B
ASV_C ───── DNA sequence C
```

MHASS cannot invent a new ASV here. It works from the sequences you provide.

---

# 3. Input 2: `amplicon_labels.tsv`

This answers another question:

> “Which genome/organism does each amplicon belong to?”

Simplified:

```tsv
asvid	genomeid
ASV_A	Bifidobacterium_adolescentis
ASV_B	Enterococcus_faecalis
ASV_C	Clostridium_beijerinckii
```

So MHASS knows:

```text
ASV_A
   ↓
Bifidobacterium

ASV_B
   ↓
Enterococcus

ASV_C
   ↓
Clostridium
```

Your actual file contains the entire FASTA header in the `asvid` column.

That is important:

```text
FASTA identifier
        =
amplicon_labels.tsv ASVID
```

They must match.

---

# 4. Why does native MHASS need the genome mapping?

Because metaSPARSim does not immediately start by saying:

```text
ASV_A = 10 reads
ASV_B = 10 reads
ASV_C = 10 reads
```

First MHASS defines the **genome-level abundance distribution**.

For example:

```text
--genome-distribution uniform
```

means approximately:

```text
Bifidobacterium     1/3
Enterococcus        1/3
Clostridium         1/3
```

Then MHASS links the ASVs to those genomes.

Conceptually:

```text
genome abundance
       ↓
ASV → genome mapping
       ↓
ASV expected abundances
       ↓
metaSPARSim
       ↓
actual integer counts
```

For your tiny test, MHASS calculated:

```text
ASV_A expected intensity = 10
ASV_B expected intensity = 10
ASV_C expected intensity = 10
```

because:

```text
30 total reads / 3 ASVs
= 10
```

But **10 is not necessarily the final count**.

It is the expected abundance used by the stochastic model.

---

# 5. What metaSPARSim was doing

This is a crucial distinction.

You told MHASS:

```text
30 total reads
3 equally represented ASVs
```

You might imagine that MHASS would simply produce:

```text
ASV_A    10
ASV_B    10
ASV_C    10
```

But that is not what metaSPARSim is designed to do.

It adds biological/sequencing variability.

So one execution can produce:

```text
ASV_A    13
ASV_B     8
ASV_C     9
```

Total:

```text
13 + 8 + 9 = 30
```

Another execution can produce:

```text
ASV_A     8
ASV_B     7
ASV_C    15
```

again:

```text
8 + 7 + 15 = 30
```

The total library size is still 30, but the community composition has changed.

Think of it like this.

You say:

> “There are three kinds of candies and they should be roughly equally common.”

metaSPARSim then randomly fills a bag:

```text
Bag 1:
13 red
 8 green
 9 blue
```

Another time:

```text
Bag 2:
 8 red
 7 green
15 blue
```

Both bags contain 30 candies.

But they are **not the same experimental community**.

---

# 6. Why metaSPARSim does this

That randomness is not inherently bad.

If your scientific question were:

> “What might realistic microbial communities look like if abundance varies naturally?”

then stochastic count simulation can be useful.

metaSPARSim is trying to model variability rather than produce a perfectly fixed artificial community.

But that is **not quite your benchmarking question**.

You are asking something much closer to:

> “If I know exactly what biological sequences should be present, how well do DADA2, VSEARCH and UNOISE3 recover them?”

For that question, changing the starting community every time is unnecessary noise.

---

# 7. What is `counts.tsv`?

It is simply a table answering:

> “How many copies of every ASV should exist in every sample?”

For our three-ASV example:

```tsv
ASVID	Sample1
ASV_A	10
ASV_B	10
ASV_C	10
```

This means:

```text
Sample1
│
├── 10 copies of ASV_A
├── 10 copies of ASV_B
└── 10 copies of ASV_C

TOTAL = 30
```

Very simple.

---

# 8. A more realistic multi-sample counts file

Suppose later you want three samples:

```tsv
ASVID	Sample1	Sample2	Sample3
ASV_A	100	20	5
ASV_B	50	100	5
ASV_C	10	20	100
ASV_D	0	10	20
```

Read it by columns.

### Sample 1

```text
ASV_A = 100
ASV_B = 50
ASV_C = 10
ASV_D = 0
```

So ASV_D is absent.

### Sample 2

```text
ASV_A = 20
ASV_B = 100
ASV_C = 20
ASV_D = 10
```

### Sample 3

```text
ASV_A = 5
ASV_B = 5
ASV_C = 100
ASV_D = 20
```

Now you know **exactly** what you intended each sample to contain.

---

# 9. How MHASS uses `counts.tsv`

This is the really important programming part.

Imagine the table says:

```tsv
ASVID	Sample1
ASV_A	3
ASV_B	2
```

MHASS's template-creation step interprets that as:

```text
ASV_A
→ make 3 simulated molecules

ASV_B
→ make 2 simulated molecules
```

So logically:

```text
counts.tsv

ASV_A = 3
ASV_B = 2

        ↓

copy ASV_A #1
copy ASV_A #2
copy ASV_A #3

copy ASV_B #1
copy ASV_B #2
```

Each copy is then assigned a simulated number of PacBio passes, called `np`.

Maybe:

```text
ASV_A copy 1 → np = 12
ASV_A copy 2 → np = 30
ASV_A copy 3 → np = 12

ASV_B copy 1 → np = 45
ASV_B copy 2 → np = 20
```

---

# 10. What is `np`?

You can think of PacBio sequencing as reading the same circular molecule several times.

If a molecule is read:

```text
12 times
```

then:

```text
np = 12
```

If it is read:

```text
45 times
```

then:

```text
np = 45
```

More passes usually give CCS more evidence to build an accurate consensus.

So MHASS creates:

```text
intended molecule
       ↓
choose np
       ↓
PBSIM simulates those passes
       ↓
CCS combines them
       ↓
HiFi-like consensus read
```

---

# 11. Why did 30 counts create only 26 template FASTA files?

You already observed this, and it makes sense once we understand the grouping.

Suppose:

```text
ASV_A copy 1 → np12
ASV_A copy 2 → np30
ASV_A copy 3 → np12
```

Copies 1 and 3 have:

```text
same ASV
same sample
same barcode
same np
```

So MHASS can put them together in one template file:

```text
template1_..._np12.fasta
```

containing two records.

Then:

```text
template1_..._np30.fasta
```

contains the other one.

Therefore:

```text
30 molecules
≠ necessarily
30 template files
```

because molecules with the same properties can be grouped.

---

# 12. Where do the barcodes come in?

You also provide:

```text
barcodes.tsv
```

Conceptually it contains something like:

```tsv
BarcodeID	Forward	Reverse
BC01	ACGTACGT	TGCATGCA
BC02	GGTTAACC	CCAATTGG
```

Those sequences are only simplified examples.

For each sample, MHASS assigns a barcode.

Imagine:

```text
Sample1 → BC01
Sample2 → BC02
```

Then a simulated template is built approximately like:

```text
forward barcode
      +
amplicon
      +
reverse barcode
```

The source actually constructs something conceptually like:

```text
A + forward_barcode + amplicon + reverse_complement(reverse_barcode) + A
```

So the template sent to PBSIM resembles what the sequencing library would contain rather than containing only the biological amplicon.

---

# 13. The template stage

After reading:

```text
amplicons.fa
counts.tsv
barcodes.tsv
```

MHASS constructs the physical-like molecules to simulate.

So this:

```text
ASV_A = 3 copies
```

becomes:

```text
copy 1:
barcode + ASV_A + barcode
np = 12

copy 2:
barcode + ASV_A + barcode
np = 30

copy 3:
barcode + ASV_A + barcode
np = 12
```

These become template FASTA records.

---

# 14. PBSIM comes next

PBSIM's job is different.

PBSIM does **not decide which taxa are present**.

That decision has already been made.

PBSIM asks:

> “If PacBio sequenced this molecule, what noisy subreads might I observe?”

Conceptually:

```text
Perfect template:

ACGTACGTACGTACGT
       ↓
PBSIM
       ↓
subread 1:
ACGTACGTTCGTACGT

subread 2:
ACGTACGTACGTACGT

subread 3:
ACGTACGTACG-ACGT
...
```

Those are simplified examples, but the point is:

```text
counts.tsv
→ determines abundance

PBSIM
→ determines sequencing errors
```

These are two completely different problems.

---

# 15. Then CCS

PBSIM generates multiple noisy passes.

CCS tries to combine them:

```text
subread 1  ACGTACGTTCGT
subread 2  ACGTACGTACGT
subread 3  ACGTACGTACGT
                 ↓
                CCS
                 ↓
          ACGTACGTACGT
```

If CCS is confident enough, the molecule becomes a HiFi read.

If not, it can be rejected.

That is why:

```text
counts.tsv = 30 molecules
```

does **not guarantee**:

```text
final FASTQ = 30 reads
```

Your smoke test demonstrated this nicely:

```text
30 intended counts
      ↓
28 grouped template files
      ↓
21 final CCS reads
```

Some molecules simply did not meet the CCS quality criterion.

---

# 16. Finally MHASS combines the reads

After CCS, there may be many small:

```text
ccs.fastq
```

files.

MHASS combines them into:

```text
combined_reads.fastq
```

and currently randomly reverse-complements roughly half of them.

So the final product is:

```text
combined_reads.fastq
```

which is what you would then send into your downstream benchmarking pipelines.

---

# 17. Now let's replace metaSPARSim

Current MHASS:

```text
amplicons.fa
     +
genome labels
     +
--num-reads 10000
     +
--num-samples 3
     |
     v
metaSPARSim
     |
     | RANDOM
     v
counts.tsv
     |
     v
template creation
     |
     v
PBSIM
     |
     v
CCS
```

Our proposed benchmark mode:

```text
amplicons.fa
      +
FIXED counts.tsv
      +
barcodes.tsv
      |
      v
template creation
      |
      v
PBSIM
      |
      v
CCS
```

We simply remove:

```text
metaSPARSim
```

from the reproducible path.

---

# 18. What code behavior do we want?

Conceptually, MHASS's Step 1 would become:

```python
if counts_file_was_given:
    use_counts_file_exactly()
else:
    run_metasparsim()
```

So:

```bash
mhass \
    --amplicon-fasta amplicons.fa \
    --counts-file counts.tsv \
    ...
```

means:

> “Do not invent my community. I have already defined it.”

Whereas somebody could still run native MHASS without:

```text
--counts-file
```

and get its normal metaSPARSim behavior.

That is good software design because we **add a capability without destroying the original one**.

---

# 19. Why is the fixed count matrix especially good for your thesis?

Suppose you want to compare:

```text
DADA2
VSEARCH
UNOISE3
```

Imagine metaSPARSim produces this for experiment A:

```text
ASV_A = 80
ASV_B = 15
ASV_C = 5
```

and this for experiment B:

```text
ASV_A = 20
ASV_B = 60
ASV_C = 20
```

Now suppose DADA2 appears better in A and VSEARCH appears better in B.

You have a problem.

Was the difference because:

```text
DADA2 vs VSEARCH?
```

or because:

```text
community A vs community B?
```

You changed **two things at once**.

That weakens an experiment.

---

# 20. A good benchmark changes one important thing at a time

Instead, define:

```text
Ground truth

ASV_A = 80
ASV_B = 15
ASV_C = 5
```

once.

Then:

```text
                    SAME READ DATA
                         |
           +-------------+-------------+
           |             |             |
           v             v             v
        DADA2         VSEARCH       UNOISE3
           |             |             |
           v             v             v
        result A       result B       result C
```

Now differences in recovery can much more confidently be attributed to the pipelines.

That's the core logic of benchmarking.

---

# 21. Think of it like testing three students

Suppose you want to compare three students' math skills.

Bad experiment:

```text
Student A → easy exam
Student B → hard exam
Student C → medium exam
```

Then comparing scores isn't fair.

Good experiment:

```text
Student A ─┐
Student B ─┼→ SAME EXAM
Student C ─┘
```

Your fixed `counts.tsv` plays the role of making sure the **starting test is the same**.

---

# 22. It also gives you a cleaner ground truth

Suppose your fixed file says:

```tsv
ASVID	Sample1
ASV_A	500
ASV_B	250
ASV_C	50
ASV_D	5
```

Then you already know:

```text
Expected richness = 4 ASVs

ASV_A = abundant
ASV_B = moderately abundant
ASV_C = low abundance
ASV_D = very rare
```

Now you can ask:

```text
Did DADA2 recover ASV_D?

Did VSEARCH lose it?

Did UNOISE3 merge it with another ASV?

Did any method invent an ASV that didn't exist?
```

Those are excellent benchmarking questions because the answer key is known.

---

# 23. You could intentionally make difficult benchmark communities

This is another major advantage.

Instead of letting metaSPARSim randomly decide what happens, you can deliberately create tests.

For example:

```text
Scenario 1 — Easy

ASV_A = 1000
ASV_B = 1000
ASV_C = 1000
```

Then:

```text
Scenario 2 — Uneven

ASV_A = 2500
ASV_B = 450
ASV_C = 50
```

Then:

```text
Scenario 3 — Rare ASV challenge

ASV_A = 2990
ASV_B = 9
ASV_C = 1
```

Now you can ask:

> At what abundance does each bioinformatics method stop reliably recovering a taxon?

That is potentially much more scientifically informative than one randomly generated community.

You do not necessarily need all these scenarios for the thesis, but the architecture enables them.

---

# 24. There are actually three different kinds of “truth”

This is important for your thesis.

### Truth level 1: intended biological community

```text
counts.tsv
```

Example:

```text
ASV_A = 100
ASV_B = 50
```

This means:

> “I asked the simulator to generate 100 molecules from A and 50 from B.”

---

### Truth level 2: simulated sequencing outcome

Maybe:

```text
ASV_A intended = 100
ASV_B intended = 50
```

but after PBSIM + CCS:

```text
ASV_A emitted HiFi = 94
ASV_B emitted HiFi = 41
```

because some molecules failed CCS.

This isn't necessarily an error.

It's part of the sequencing simulation.

---

### Truth level 3: inferred pipeline result

DADA2 might finally report:

```text
ASV_A = 92
ASV_B = 39
ASV_X = 2
```

Now you can compare:

```text
intended community
       ↓
realized HiFi reads
       ↓
pipeline inference
```

That is a very powerful benchmarking framework.

---

# 25. This is why I want a read-level manifest later

You should eventually be able to trace something like:

```text
read_000001
│
├── Sample1
├── source ASV = ASV_A
├── intended copy = 37
├── np = 24
├── PBSIM seed = ...
├── CCS passed = yes
└── orientation = reverse
```

Then if a pipeline makes a surprising mistake, you can trace what happened.

That is much stronger than having only a FASTQ file whose origin you cannot fully reconstruct.

---

# 26. What inputs would the reproducible mode ultimately need?

Conceptually:

```text
amplicons.fa
      │
      │ What sequences exist?
      │
      v

counts.tsv
      │
      │ How many copies exist in each sample?
      │
      v

barcodes.tsv
      │
      │ Which barcode identifies each sample?
      │
      v

--seed 12345
      │
      │ How should all remaining randomness behave?
      │
      v

MHASS
      │
      v

simulated HiFi FASTQ
```

That's a very clean experimental definition.

---

# 27. What happens to `amplicon_labels.tsv`?

This is an interesting design question.

In **native MHASS mode**, it is definitely needed because MHASS uses:

```text
ASV
↓
genome
↓
genome abundance model
↓
metaSPARSim counts
```

But in fixed-count mode:

```text
counts.tsv
```

already says exactly how many copies every ASV gets.

So technically, the genome mapping is no longer required to calculate counts.

I would probably design it so that:

```text
native mode:
--amplicon-genome-labels required

fixed-count mode:
--amplicon-genome-labels optional
```

You might still keep the mapping for metadata and validation, but it should not determine abundance anymore.

We should decide that carefully when implementing `--counts-file`.

---

# 28. What happens to `--num-reads`?

This is another important consequence.

With metaSPARSim:

```bash
--num-reads 10000
```

means roughly:

> “Construct a count matrix whose library size is 10,000 per sample.”

But if you supply:

```text
counts.tsv
```

then the file itself contains the library size.

For example:

```tsv
ASVID	Sample1
ASV_A	6000
ASV_B	3000
ASV_C	1000
```

already tells MHASS:

```text
6000 + 3000 + 1000
= 10000
```

So in fixed-count mode, I would **not let `--num-reads` override the matrix**.

The matrix should be authoritative.

MHASS could simply validate and report:

```text
Sample1 library size = 10000
```

If both are supplied, I'd prefer either:

```text
ignore --num-reads with a clear warning
```

or, even better:

```text
reject the ambiguous combination
```

depending on how we design the interface.

---

# 29. What happens to `--num-samples`?

Same idea.

If:

```tsv
ASVID	Sample_A	Sample_B	Sample_C
...
```

then the matrix already tells us:

```text
3 samples
```

Therefore:

```text
--num-samples
```

should not independently control sample number in fixed-count mode.

The matrix becomes the source of truth.

---

# 30. So the modes become conceptually very clean

### Native MHASS

```text
amplicons.fa
amplicon_labels.tsv
--num-samples
--num-reads
--genome-distribution
        |
        v
metaSPARSim
        |
        v
counts.tsv
```

### Thesis benchmark mode

```text
amplicons.fa
counts.tsv
        |
        v
validate exact ASV/sample counts
        |
        v
template creation
```

That is simpler.

And simpler experimental design is usually easier to defend.

---

# What I want you to remember most

If you remember only four ideas, remember these:

1. **`amplicons.fa` says WHAT sequences exist.**
2. **`counts.tsv` says HOW MANY copies of each sequence exist in each sample.**
3. **PBSIM + CCS say WHAT sequencing those molecules looks like.**
4. **DADA2/VSEARCH/UNOISE3 are then tested on those simulated sequencing reads.**

So your thesis benchmark becomes:

```text
KNOWN SEQUENCES
      +
KNOWN ABUNDANCES
      ↓
SIMULATED SEQUENCING
      ↓
UNKNOWN TO PIPELINE
      ↓
PIPELINE TRIES TO RECONSTRUCT TRUTH
      ↓
COMPARE INFERENCE WITH KNOWN TRUTH
```

That is the essence of a **ground-truth benchmarking experiment**.

And this is why I think the fixed counts matrix is not merely an escape from the metaSPARSim reproducibility problem. For your specific thesis question, it gives you **better experimental control and a clearer answer key**.

### English Corrections

* **Mistake 1**

  * Original: “explain to me the data flow in MHASS and how are the inputs used.”
  * Corrected: “explain to me the data flow in MHASS and how the inputs are used.”
  * Mistake: In an embedded question, English uses statement word order: `how the inputs are used`.
  * Type: Word order.

* **Mistake 2**

  * Original: “Explain how they a fixed counts file would look”
  * Corrected: “Explain what a fixed counts file would look like”
  * Mistake: `how they` appears to be a typo, and the natural construction is **what X would look like**.
  * Type: Grammar / typo.

* **Mistake 3**

  * Original: “why this is good for my thesis vs what metaSPARSim used to do.”
  * Corrected: “why this is better for my thesis than what metaSPARSim used to do.”
  * Mistake: `vs` is understandable, but **better ... than ...** is more natural in formal English.
  * Type: Style / comparison.

