#!/usr/bin/env python3
import os
import random
from pathlib import Path
from tqdm import tqdm

from mhass.seed_utils import derive_seed

def reverse_complement(seq):
    """Generate reverse complement of a DNA sequence."""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 
                  'a': 't', 't': 'a', 'c': 'g', 'g': 'c',
                  'N': 'N', 'n': 'n'}
    
    # Handle any non-standard characters by keeping them as-is
    return ''.join(complement.get(base, base) for base in reversed(seq))

def combine_fastqs(
    input_dir,
    output_fastq,
    master_seed=None,
):
    """Combine FASTQ files and add template prefix to headers."""
    input_dir = Path(input_dir)
    total_records = 0
    total_files = 0
    total_revcomp = 0

    orientation_rng = None

    if master_seed is not None:
        orientation_seed = derive_seed(
            master_seed,
            "orientation",
        )
        orientation_rng = random.Random(
            orientation_seed
        )

    # First, count total files to process for progress display
    subdirs = [d for d in input_dir.iterdir() if d.is_dir()]
    total_subdirs = len(subdirs)
    
    print(f"Found {total_subdirs} template directories to process")
    
    with open(output_fastq, 'w') as outfile:
        for subdir in tqdm(sorted(subdirs), desc="Combining FASTQ files"):
            template_name = subdir.name
            fastq_file = subdir / "ccs.fastq"
            
            if not fastq_file.exists():
                print(f"[Warning] Missing FASTQ: {fastq_file}")
                continue
            
            total_files += 1
            records_in_file = 0
            
            with open(fastq_file, 'r') as f:
                while True:
                    header = f.readline().strip()
                    if not header:
                        break  # end of file
                    
                    seq = f.readline().strip()
                    plus = f.readline().strip()
                    qual = f.readline().strip()

                    # Making the choice of reversing the read explicit
                    if orientation_rng is None:
                        # original MHASS behavior
                        reverse_this_read = random.random() < 0.5
                    else:
                        # deterministic local RNG
                        reverse_this_read = orientation_rng.random() < 0.5

                    if reverse_this_read:
                        seq = reverse_complement(seq)
                        qual = qual[::-1]  # Reverse quality string
                        total_revcomp += 1

                    # Rewrite header to prefix with template name
                    if header.startswith('@'):
                        new_header = f"@{template_name}/{header[1:]}"
                    else:
                        new_header = f"{template_name}/{header}"
                    
                    # Write to output with proper FASTQ format
                    outfile.write(f"{new_header}\n{seq}\n{plus}\n{qual}\n")
                    
                    records_in_file += 1
                    total_records += 1
    
    print(f"[Done] Combined {total_records} records from {total_files} FASTQ files into: {output_fastq}")
    print(f"       Reverse complemented {total_revcomp} sequences ({total_revcomp/total_records*100:.1f}%)")
    return output_fastq
