from mhass.template_creator import create_per_sequence_templates


def make_template_inputs(base_dir):
    base_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    fasta = base_dir / "amplicons.fa"
    counts = base_dir / "counts.tsv"
    barcodes = base_dir / "barcodes.tsv"
    output_dir = base_dir / "templates"
    barcode_mapping = base_dir / "sample_barcode_map.tsv"

    fasta.write_text(
        ">ASV1\n"
        "ACGTACGTACGT\n"
        ">ASV2\n"
        "TTTTCCCCAAAA\n"
    )

    counts.write_text(
        "ASVID\tSample1\tSample2\n"
        "ASV1\t5\t3\n"
        "ASV2\t4\t2\n"
    )

    barcodes.write_text(
        "id\tforward\treverse\n"
        "BC01\tAAAA\tCCCC\n"
        "BC02\tGGGG\tTTTT\n"
    )

    np_params = {
        "distribution_type": "lognormal",
        "empirical_file": None,
        "lognormal_mu": 3.88,
        "lognormal_sigma": 1.22,
        "np_min": 2,
        "np_max": 59,
    }

    return {
        "fasta": fasta,
        "counts": counts,
        "barcodes": barcodes,
        "output_dir": output_dir,
        "barcode_mapping": barcode_mapping,
        "np_params": np_params,
    }


def run_template_creation(base_dir, seed):
    inputs = make_template_inputs(base_dir)

    create_per_sequence_templates(
        inputs["fasta"],
        inputs["counts"],
        inputs["output_dir"],
        inputs["barcodes"],
        inputs["barcode_mapping"],
        inputs["np_params"],
        master_seed=seed,
    )

    template_contents = {
        path.name: path.read_text()
        for path in sorted(
            inputs["output_dir"].glob("*.fasta")
        )
    }

    sequence_mapping = (
        inputs["output_dir"]
        / "sequence_file_mapping.tsv"
    ).read_text()

    return {
        "templates": template_contents,
        "sequence_mapping": sequence_mapping,
    }


def test_template_creation_same_seed_is_identical(tmp_path):
    first = run_template_creation(
        tmp_path / "run_a",
        12345,
    )

    second = run_template_creation(
        tmp_path / "run_b",
        12345,
    )

    assert first == second


def test_template_creation_different_seed_changes_result(tmp_path):
    seed_12345 = run_template_creation(
        tmp_path / "run_12345",
        12345,
    )

    seed_67890 = run_template_creation(
        tmp_path / "run_67890",
        67890,
    )

    assert seed_12345 != seed_67890


def test_template_creation_maps_every_intended_copy(tmp_path):
    result = run_template_creation(
        tmp_path / "run",
        12345,
    )

    mapping_lines = (
        result["sequence_mapping"]
        .strip()
        .splitlines()
    )

    assert len(mapping_lines) == 15


def test_template_creation_unseeded_call_still_works(tmp_path):
    result = run_template_creation(
        tmp_path / "unseeded",
        None,
    )

    assert result["templates"]
    assert result["sequence_mapping"]
