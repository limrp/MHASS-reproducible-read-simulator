from mhass.template_creator import create_per_sequence_templates


def test_headerless_barcode_file_preserves_first_barcode(tmp_path):
    fasta = tmp_path / "amplicons.fa"
    counts = tmp_path / "counts.tsv"
    barcodes = tmp_path / "barcodes.tsv"
    templates_dir = tmp_path / "templates"
    barcode_mapping = tmp_path / "sample_barcode_map.tsv"

    fasta.write_text(
        ">ASV1\n"
        "ACGTACGTACGT\n"
    )

    counts.write_text(
        "ASVID\tSample1\tSample2\n"
        "ASV1\t1\t1\n"
    )

    # Deliberately headerless.
    #
    # A third barcode is included so that the current bug does not merely
    # trigger a "not enough barcodes" error. Instead, it silently shifts
    # the assignments:
    #
    # expected: Sample1 -> A1, Sample2 -> A2
    # buggy:    Sample1 -> A2, Sample2 -> A3
    barcodes.write_text(
        "A1\tAAAA\tCCCC\n"
        "A2\tGGGG\tTTTT\n"
        "A3\tACAC\tGTGT\n"
    )

    np_params = {
        "distribution_type": "lognormal",
        "empirical_file": None,
        "lognormal_mu": 3.88,
        "lognormal_sigma": 1.22,
        "np_min": 2,
        "np_max": 59,
    }

    create_per_sequence_templates(
        fasta,
        counts,
        templates_dir,
        barcodes,
        barcode_mapping,
        np_params,
        master_seed=12345,
    )

    mapping_lines = barcode_mapping.read_text().strip().splitlines()

    sample1 = mapping_lines[1].split("\t")
    sample2 = mapping_lines[2].split("\t")

    assert sample1[0] == "Sample1"
    assert sample1[1] == "A1"

    assert sample2[0] == "Sample2"
    assert sample2[1] == "A2"


def test_headered_barcode_file_skips_header_and_preserves_first_barcode(tmp_path):
    fasta = tmp_path / "amplicons.fa"
    counts = tmp_path / "counts.tsv"
    barcodes = tmp_path / "barcodes.tsv"
    templates_dir = tmp_path / "templates"
    barcode_mapping = tmp_path / "sample_barcode_map.tsv"

    fasta.write_text(
        ">ASV1\n"
        "ACGTACGTACGT\n"
    )

    counts.write_text(
        "ASVID\tSample1\tSample2\n"
        "ASV1\t1\t1\n"
    )

    barcodes.write_text(
        "BarcodeID\tForwardBarcode\tReverseBarcode\n"
        "A1\tAAAA\tCCCC\n"
        "A2\tGGGG\tTTTT\n"
    )

    np_params = {
        "distribution_type": "lognormal",
        "empirical_file": None,
        "lognormal_mu": 3.88,
        "lognormal_sigma": 1.22,
        "np_min": 2,
        "np_max": 59,
    }

    create_per_sequence_templates(
        fasta,
        counts,
        templates_dir,
        barcodes,
        barcode_mapping,
        np_params,
        master_seed=12345,
    )

    mapping_lines = barcode_mapping.read_text().strip().splitlines()

    sample1 = mapping_lines[1].split("\t")
    sample2 = mapping_lines[2].split("\t")

    assert sample1[0] == "Sample1"
    assert sample1[1] == "A1"

    assert sample2[0] == "Sample2"
    assert sample2[1] == "A2"
