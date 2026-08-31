#!/usr/bin/env python3

import os
import random
import numpy as np
from collections import defaultdict
from pathlib import Path
from mhass.seed_utils import derive_seed

def reverse_complement(seq):
    """Return the reverse complement of a DNA sequence."""
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 
                  'a': 't', 'c': 'g', 'g': 'c', 't': 'a',
                  'N': 'N', 'n': 'n', 'Y': 'R', 'R': 'Y',
                  'S': 'S', 'W': 'W', 'K': 'M', 'M': 'K',
                  'B': 'V', 'V': 'B', 'D': 'H', 'H': 'D'}
    return ''.join(complement.get(base, base) for base in reversed(seq))

def load_np_distribution(np_file):
    """Load np distribution as a weighted list of values for random sampling."""
    weighted = []
    with open(np_file) as f:
        for line in f:
            if line.strip() and not line.startswith("#"):
                np_val, count = map(int, line.strip().split('\t'))
                weighted.extend([np_val] * count)
    return weighted

def sample_np_lognormal(mu, sigma, np_min, np_max, rng=None):
    """Sample np value from lognormal distribution with bounds."""
    while True:
        if rng is None:
            value = np.random.lognormal(mu, sigma)
        else:
            value = rng.lognormal(mu, sigma)
        # Round to integer and apply bounds
        np_val = int(round(value))
        if np_min <= np_val <= np_max:
            return np_val

def create_np_sampler(np_params, master_seed=None):
    """Create appropriate np sampler based on distribution type."""
    dist_type = np_params['distribution_type']
    
    if dist_type == 'empirical':
        # Load empirical distribution
        np_values = load_np_distribution(np_params['empirical_file'])
        if master_seed is None:
            return lambda: random.choice(np_values)

        empirical_seed = derive_seed(master_seed, "np-empirical")
        empirical_rng = random.Random(empirical_seed)

        return lambda: empirical_rng.choice(np_values)
    
    elif dist_type == 'lognormal':
        # Create lognormal sampler
        mean_np = np_params['lognormal_mu']
        sd_np = np_params['lognormal_sigma']
        np_min = np_params['np_min']
        np_max = np_params['np_max']
        
        print(f"Using lognormal distribution: mu={mean_np}, sigma={sd_np}, range=[{np_min}, {np_max}]")

        if master_seed is None:
            return lambda: sample_np_lognormal(
                mean_np,
                sd_np,
                np_min,
                np_max,
            )

        lognormal_seed = derive_seed(master_seed, "np-lognormal")
        lognormal_rng = np.random.default_rng(lognormal_seed)

        return lambda: sample_np_lognormal(
            mean_np,
            sd_np,
            np_min,
            np_max,
            rng=lognormal_rng,
        )
    
    else:
        raise ValueError(f"Unknown distribution type: {dist_type}")

def create_per_sequence_templates(
    fasta_path,
    counts_path,
    output_dir,
    barcode_file,
    barcode_mapping_file,
    np_params,
    master_seed=None,
):
    """Create per-sequence template files with barcodes and sampled np values for PacBio simulation."""
    # Load FASTA
    seqs = {}
    with open(fasta_path) as f:
        current_id = None
        seq_lines = []
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id:
                    seqs[current_id] = ''.join(seq_lines)
                current_id = line[1:]
                seq_lines = []
            else:
                seq_lines.append(line)
        if current_id:
            seqs[current_id] = ''.join(seq_lines)
    print(f"Loaded {len(seqs)} sequences from {fasta_path}")

    # Load counts
    with open(counts_path) as f:
        lines = [line.strip().split('\t') for line in f if line.strip()]
    sample_names = lines[0][1:]
    count_data = {row[0]: list(map(int, row[1:])) for row in lines[1:]}

    # Create np value sampler
    sample_np = create_np_sampler(
        np_params,
        master_seed=master_seed,
    )

    # Load barcodes
    barcodes = []
    with open(barcode_file) as f:
        next(f)  # skip header
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                barcodes.append({'id': parts[0], 'forward': parts[1], 'reverse': parts[2]})
    if len(barcodes) < len(sample_names):
        raise ValueError(f"Only {len(barcodes)} barcodes for {len(sample_names)} samples")

    os.makedirs(output_dir, exist_ok=True)

    # Write barcode mapping file
    with open(barcode_mapping_file, 'w') as mapf:
        mapf.write("SampleID\tBarcodeID\tForwardBarcode\tReverseBarcode\tRevCompReverse\n")
        for i, sample in enumerate(sample_names):
            bc = barcodes[i]
            revcomp = reverse_complement(bc['reverse'])
            mapf.write(f"{sample}\t{bc['id']}\t{bc['forward']}\t{bc['reverse']}\t{revcomp}\n")

    # Buffers for grouped output by (template, barcode, np)
    output_buffers = defaultdict(list)
    sequence_mapping = []
    template_counter = 1

    for asv_id, seq in seqs.items():
        if asv_id not in count_data:
            print(f"Skipping {asv_id}: not found in count matrix")
            continue

        template_prefix = f"template{template_counter}"
        template_counter += 1
        counts = count_data[asv_id]

        for i, count in enumerate(counts):
            if count == 0:
                continue
            sample = sample_names[i]
            bc = barcodes[i]
            revcomp = reverse_complement(bc['reverse'])
            for j in range(count):
                np_val = sample_np()  # Use the appropriate sampler
                barcoded_seq = f"A{bc['forward']}{seq}{revcomp}A"
                filename = f"{template_prefix}_{bc['id']}_np{np_val}.fasta"
                header = f">{asv_id}::{sample}::Copy{j+1}::Barcode{bc['id']}"
                output_buffers[filename].append((header, barcoded_seq))
                sequence_mapping.append((asv_id, sample, bc['id'], f"{filename}"))

    # Write all grouped FASTA files
    for filename, records in output_buffers.items():
        full_path = os.path.join(output_dir, filename)
        with open(full_path, 'w') as f:
            for header, seq in records:
                f.write(f"{header}\n{seq}\n")

    # Write sequence mapping file
    mapping_file = os.path.join(output_dir, "sequence_file_mapping.tsv")
    with open(mapping_file, 'w') as mf:
        mf.write("ASVID\tSampleID\tBarcodeID\tTemplateFile\n")
        for row in sequence_mapping:
            mf.write('\t'.join(row) + '\n')

    print(f"Wrote ASV×Sample×np-to-template mapping to {mapping_file}")
    print(f"Created {len(output_buffers)} template files in {output_dir}")
    
    return mapping_file
