#!/usr/bin/env python3

import argparse
import csv
import os
import sys
import subprocess
import shutil
import tempfile
from pathlib import Path
import multiprocessing

# Import local modules (when running as a package)
try:
    from mhass.fastq_combiner import combine_fastqs
    from mhass.pbsim_runner import run_pbsim
    from mhass.template_creator import create_per_sequence_templates
except ImportError:
    # For development/testing
    from fastq_combiner import combine_fastqs
    from pbsim_runner import run_pbsim
    from template_creator import create_per_sequence_templates


def get_package_path():
    """Get the path to the package resources."""
    return Path(__file__).resolve().parent


def get_resource_path(resource_name):
    """Get the path to a resource file bundled with the package."""
    return get_package_path() / "resources" / resource_name


def prepare_external_counts(counts_input, amplicon_fasta, output_counts, output_meta):
    """
    Validate an external ASV-by-sample count matrix and copy it
    unchanged into the MHASS output directory.
    """
    counts_input = Path(counts_input)
    amplicon_fasta = Path(amplicon_fasta)
    output_counts = Path(output_counts)
    output_meta = Path(output_meta)

    if not counts_input.is_file():
        raise ValueError(f"Counts file not found: {counts_input}")

    # ------------------------------------------------------------
    # Read all FASTA identifiers.
    # MHASS uses the complete FASTA header (without '>') as ASVID.
    # ------------------------------------------------------------
    fasta_ids = set()

    with amplicon_fasta.open() as fasta:
        for line in fasta:
            if line.startswith(">"):
                fasta_ids.add(line[1:].strip())

    if not fasta_ids:
        raise ValueError(f"No FASTA records found in: {amplicon_fasta}")

    # ------------------------------------------------------------
    # Validate the count matrix.
    # ------------------------------------------------------------
    seen_asvs = set()
    sample_totals = None
    sample_names = None

    with counts_input.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")

        try:
            header = next(reader)
        except StopIteration:
            raise ValueError("Counts file is empty")

        if len(header) < 2:
            raise ValueError(
                "Counts file must contain an ASVID column "
                "and at least one sample column"
            )

        if header[0] != "ASVID":
            raise ValueError("First column of counts file must be named 'ASVID'")

        sample_names = header[1:]

        if any(not name.strip() for name in sample_names):
            raise ValueError("Sample names cannot be empty")

        if len(sample_names) != len(set(sample_names)):
            raise ValueError("Sample names must be unique")

        sample_totals = [0] * len(sample_names)

        for line_number, row in enumerate(reader, start=2):

            if len(row) != len(header):
                raise ValueError(
                    f"Line {line_number}: expected {len(header)} columns, "
                    f"found {len(row)}"
                )

            asv_id = row[0]

            if not asv_id:
                raise ValueError(f"Line {line_number}: empty ASVID")

            if asv_id in seen_asvs:
                raise ValueError(f"Line {line_number}: duplicate ASVID: {asv_id}")

            seen_asvs.add(asv_id)

            if asv_id not in fasta_ids:
                raise ValueError(
                    f"Line {line_number}: ASVID not found in amplicon FASTA: "
                    f"{asv_id}"
                )

            for i, value in enumerate(row[1:]):

                try:
                    count = int(value)
                except ValueError:
                    raise ValueError(
                        f"Line {line_number}, sample {sample_names[i]}: "
                        f"count must be an integer, found '{value}'"
                    )

                if count < 0:
                    raise ValueError(
                        f"Line {line_number}, sample {sample_names[i]}: "
                        "count cannot be negative"
                    )

                sample_totals[i] += count

    if not seen_asvs:
        raise ValueError("Counts file contains no ASV rows")

    # Check for ASVs that exist in the FASTA but were not found in the count matrix.
    missing_asvs = fasta_ids - seen_asvs

    if missing_asvs:
        missing_sorted = sorted(missing_asvs)
        preview = ", ".join(missing_sorted[:3])

        if len(missing_sorted) > 3:
            preview += f", ... and {len(missing_sorted) - 3} more"

        raise ValueError(
            f"Counts file is missing {len(missing_asvs)} FASTA ASV(s): "
            f"{preview}. Every FASTA ASV must be represented explicitly "
            "in the count matrix; use a count of 0 when an ASV is "
            "intentionally absent from a sample."
        )
    # ------------------------------------------------------------
    # Preserve the user's count matrix exactly.
    # ------------------------------------------------------------
    if counts_input.resolve() != output_counts.resolve():
        shutil.copyfile(counts_input, output_counts)

    # ------------------------------------------------------------
    # Create metadata equivalent to the useful part of
    # metaSPARSim's counts_meta.tsv.
    # ------------------------------------------------------------
    with output_meta.open("w") as handle:
        handle.write("Sample\tLibrarySize\n")

        for sample, total in zip(sample_names, sample_totals):
            handle.write(f"{sample}\t{total}\n")

    print(f"Using external count matrix: {counts_input}")
    print(f"Validated {len(seen_asvs)} ASVs across {len(sample_names)} samples")

    for sample, total in zip(sample_names, sample_totals):
        print(f"  {sample}: {total} counts")

    print(f"Copied count matrix to: {output_counts}")
    print(f"Wrote metadata to: {output_meta}")


def run_r_count_script(
    fasta_path,
    genome_map,
    num_samples,
    num_reads,
    var_intercept,
    var_slope,
    genome_distribution,
    output_counts,
    output_meta,
):
    """Run the R script to generate count matrix."""
    r_script_path = get_package_path() / "get_counts.R"

    cmd = [
        "Rscript",
        str(r_script_path),
        "-f",
        str(fasta_path),
        "-n",
        str(num_samples),
        "-r",
        str(num_reads),
        "--var-intercept",
        str(var_intercept),
        "--var-slope",
        str(var_slope),
        "-G",
        str(genome_map),
        "--genome-distribution",
        genome_distribution,
        "-o",
        str(output_counts),
        "-m",
        str(output_meta),
    ]

    print(f"Running command: {' '.join(cmd)}")
    result = subprocess.run(cmd, check=True)
    return result.returncode == 0


def run_simulation(args):
    """Run the entire simulation pipeline."""
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Set up file paths
    counts_file = output_dir / "counts.tsv"
    meta_file = output_dir / "counts_meta.tsv"
    templates_dir = output_dir / "templates"
    barcode_mapping = output_dir / "sample_barcode_map.tsv"
    sim_reads_dir = output_dir / "sim_reads"
    combined_reads = output_dir / "combined_reads.fastq"

    # Get default resources
    default_barcode_file = get_resource_path("barcodes.txt")
    default_np_dist = get_resource_path("np_distribution.tsv")

    # Use defaults or provided values
    barcode_file = args.barcode_file or default_barcode_file
    np_distribution = args.np_distribution or default_np_dist

    # Step 1: Obtain count matrix
    if args.counts_file:

        print("\n==> STEP 1: Using external count matrix <==\n")

        try:
            prepare_external_counts(
                args.counts_file, args.amplicon_fasta, counts_file, meta_file
            )
        except (ValueError, OSError) as exc:
            print(f"Error: {exc}", file=sys.stderr)
            return False

    else:

        print("\n==> STEP 1: Simulating count matrix <==\n")

        success = run_r_count_script(
            args.amplicon_fasta,
            args.amplicon_genome_labels,
            args.num_samples,
            args.num_reads,
            args.var_intercept,
            args.var_slope,
            args.genome_distribution,
            counts_file,
            meta_file,
        )

        if not success:
            print("Error: Failed to generate count matrix")
            return False

    # Step 2: Create templates
    print("\n==> STEP 2: Creating templates <==\n")
    templates_dir.mkdir(parents=True, exist_ok=True)

    # Prepare np distribution parameters
    np_params = {
        "distribution_type": args.np_distribution_type,
        "empirical_file": (
            np_distribution if args.np_distribution_type == "empirical" else None
        ),
        "lognormal_mu": args.lognormal_mu,
        "lognormal_sigma": args.lognormal_sigma,
        "np_min": args.np_min,
        "np_max": args.np_max,
    }

    create_per_sequence_templates(
        args.amplicon_fasta,
        counts_file,
        templates_dir,
        barcode_file,
        barcode_mapping,
        np_params,  # Pass the np_params dict instead of just the file
    )

    # Step 3: Simulate and process reads
    print("\n==> STEP 3: Simulating and processing reads in parallel <==\n")
    # Get the path to PBSIM and QSHMM from package resources
    pbsim_path = get_resource_path("pbsim3/src/pbsim")
    qshmm_path = get_resource_path("pbsim3/data/QSHMM-RSII.model")

    # If running in development mode, try to find PBSIM in system PATH
    if not pbsim_path.exists():
        pbsim_path = shutil.which("pbsim") or "pbsim3/src/pbsim"
        qshmm_path = Path(os.environ.get("PBSIM_QSHMM", "pbsim3/data/QSHMM-RSII.model"))

    sim_reads_dir.mkdir(parents=True, exist_ok=True)
    run_pbsim(
        templates_dir,
        sim_reads_dir,
        pbsim_path,
        qshmm_path,
        args.threads,
        args.subread_accuracy,
        args.difference_ratio,
    )

    # Step 4: Combine and relabel CCS reads
    print("\n==> STEP 4: Combining and relabeling CCS reads <==\n")
    combine_fastqs(sim_reads_dir, combined_reads)

    print("\n==> STEP 5: Cleaning up intermediate files <==\n")
    # Move sequence mapping file up one directory
    sequence_mapping_src = templates_dir / "sequence_file_mapping.tsv"
    sequence_mapping_dest = output_dir / "sequence_file_mapping.tsv"

    if sequence_mapping_src.exists():
        shutil.move(str(sequence_mapping_src), str(sequence_mapping_dest))
        print(f"Moved sequence mapping file to: {sequence_mapping_dest}")
    else:
        print("Warning: sequence_file_mapping.tsv not found in templates directory")

    # Remove templates directory
    if templates_dir.exists():
        shutil.rmtree(templates_dir)
        print(f"Removed templates directory: {templates_dir}")

    # Remove sim_reads directory
    if sim_reads_dir.exists():
        shutil.rmtree(sim_reads_dir)
        print(f"Removed sim_reads directory: {sim_reads_dir}")

    print("\n==> FINISHED <==\n")
    print(f"Final output: {combined_reads}")
    return True


def main():
    parser = argparse.ArgumentParser(
        description="Microbiome Sequencing Simulator - A complete pipeline for simulating PacBio amplicon data",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Required arguments
    parser.add_argument(
        "--amplicon-fasta", required=True, help="Input FASTA file of amplicons"
    )

    parser.add_argument(
        "--amplicon-genome-labels",
        default=None,  # this is now OPTIONAL
        help=(
            "TSV file mapping amplicons to genomes (asvid<TAB>genomeid). "
            "Required when --counts-file is not supplied."
        ),
    )
    
    parser.add_argument(
        "--output-dir", required=True, help="Output directory for all files"
    )

    # Optional arguments
    parser.add_argument(
        "--counts-file",
        default=None,
        help=(
            "Optional precomputed ASV-by-sample count matrix (TSV). "
            "When supplied, skip metaSPARSim count generation and use "
            "these counts exactly."
        ),
    )

    parser.add_argument(
        "--num-samples", type=int, default=10, help="Number of samples to simulate"
    )
    parser.add_argument(
        "--num-reads", type=int, default=10000, help="Number of reads per sample"
    )
    parser.add_argument(
        "--var-intercept",
        type=float,
        default=1.47565981333483,
        help="Intercept for ASV variability model (controls baseline variation between samples)",
    )
    parser.add_argument(
        "--var-slope",
        type=float,
        default=-0.909890963463704,
        help="Slope for ASV variability model (how variation changes with abundance)",
    )
    parser.add_argument(
        "--genome-distribution",
        default="uniform",
        help="Distribution for genome abundances (uniform, lognormal, powerlaw, or empirical:<file>)",
    )
    parser.add_argument(
        "--barcode-file",
        default=None,
        help="TSV file with barcodes (id, forward, reverse)",
    )
    parser.add_argument(
        "--subread-accuracy",
        type=float,
        default=0.65,
        help="Mean subread accuracy used in PBSIM (default: 0.65)",
    )
    parser.add_argument(
        "--difference-ratio",
        type=str,
        default="6:55:39",
        help="Difference (error) ratio for PBSIM (substitution:insertion:deletion). "
        "Default 6:55:39 is for PacBio RS II. Use 22:45:33 for PacBio Sequel, "
        "39:24:36 for ONT",
    )
    parser.add_argument(
        "--np-distribution-type",
        default="empirical",
        choices=["empirical", "lognormal"],
        help="Type of np distribution: empirical (from file) or lognormal",
    )
    parser.add_argument(
        "--lognormal-mu",
        type=float,
        default=3.88,
        help="μ parameter for lognormal distribution of Number of Passes (default from empirical fit)",
    )
    parser.add_argument(
        "--lognormal-sigma",
        type=float,
        default=1.22,
        help="σ parameter for lognormal distribution of Num Passes (default from empirical fit)",
    )
    parser.add_argument(
        "--np-min",
        type=int,
        default=2,
        help="Minimum np value when using lognormal distribution",
    )
    parser.add_argument(
        "--np-max",
        type=int,
        default=59,
        help="Maximum np value when using lognormal distribution",
    )
    parser.add_argument(
        "--np-distribution",
        default=None,
        help="TSV file with empirical num-passes distribution",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=multiprocessing.cpu_count(),
        help="Number of threads to use for parallel processing",
    )
    args = parser.parse_args()

    if args.counts_file is None and args.amplicon_genome_labels is None:
        parser.error(
            "--amplicon-genome-labels is required unless --counts-file is supplied"
        )

    success = run_simulation(args)
    return 0 if success else 1


if __name__ == "__main__":
    sys.exit(main())
