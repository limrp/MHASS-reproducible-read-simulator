from mhass.fastq_combiner import combine_fastqs


MASTER_SEED_A = 12345
MASTER_SEED_B = 67890


def write_fastq(
    path,
    template_label,
    n_reads,
):
    with open(path, "w") as handle:
        for i in range(1, n_reads + 1):
            seq = (
                f"ACGTTGCAACGT{i % 10}A"
                .replace("0", "A")
                .replace("1", "C")
                .replace("2", "G")
                .replace("3", "T")
                .replace("4", "A")
                .replace("5", "C")
                .replace("6", "G")
                .replace("7", "T")
                .replace("8", "A")
                .replace("9", "C")
            )

            qual = "".join(
                chr(33 + ((i + j) % 40))
                for j in range(len(seq))
            )

            handle.write(
                f"@{template_label}_read{i}\n"
                f"{seq}\n"
                "+\n"
                f"{qual}\n"
            )


def create_input(input_dir):
    template_b = input_dir / "templateB"
    template_a = input_dir / "templateA"

    # Deliberately create B before A.
    template_b.mkdir(parents=True)
    template_a.mkdir(parents=True)

    write_fastq(
        template_a / "ccs.fastq",
        "A",
        20,
    )

    write_fastq(
        template_b / "ccs.fastq",
        "B",
        20,
    )


def test_fastq_combiner_same_seed_is_bitwise_identical(
    tmp_path,
):
    input_dir = tmp_path / "input"
    create_input(input_dir)

    first_output = tmp_path / "first.fastq"
    second_output = tmp_path / "second.fastq"

    combine_fastqs(
        input_dir,
        first_output,
        master_seed=MASTER_SEED_A,
    )

    combine_fastqs(
        input_dir,
        second_output,
        master_seed=MASTER_SEED_A,
    )

    assert (
        first_output.read_bytes()
        == second_output.read_bytes()
    )


def test_fastq_combiner_different_seed_changes_output(
    tmp_path,
):
    input_dir = tmp_path / "input"
    create_input(input_dir)

    seed_a_output = tmp_path / "seed_a.fastq"
    seed_b_output = tmp_path / "seed_b.fastq"

    combine_fastqs(
        input_dir,
        seed_a_output,
        master_seed=MASTER_SEED_A,
    )

    combine_fastqs(
        input_dir,
        seed_b_output,
        master_seed=MASTER_SEED_B,
    )

    assert (
        seed_a_output.read_bytes()
        != seed_b_output.read_bytes()
    )


def test_fastq_combiner_unseeded_call_still_works(
    tmp_path,
):
    input_dir = tmp_path / "input"
    create_input(input_dir)

    output = tmp_path / "unseeded.fastq"

    combine_fastqs(
        input_dir,
        output,
        master_seed=None,
    )

    assert output.exists()
    assert output.stat().st_size > 0

    lines = output.read_text().splitlines()

    assert len(lines) == 40 * 4
