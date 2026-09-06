from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from mhass.pbsim_runner import run_pbsim
from mhass.seed_utils import derive_pbsim_seed


MASTER_SEED = 12345


def create_templates(template_dir, n=12):
    template_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    filenames = []

    for i in range(1, n + 1):
        np_value = 5 + i
        filename = (
            f"template{i}_BC01_np{np_value}.fasta"
        )

        path = template_dir / filename

        path.write_text(
            f">template{i}\n"
            "ACGTACGTACGT\n"
        )

        filenames.append(filename)

    return filenames


def capture_seed_mapping(
    template_dir,
    output_dir,
    threads,
    master_seed,
):
    pbsim_commands = []

    def fake_subprocess_run(
        cmd,
        *args,
        **kwargs,
    ):
        if "--strategy" in cmd:
            pbsim_commands.append(
                list(cmd)
            )

        return SimpleNamespace(
            returncode=0,
            stderr="",
        )

    with patch(
        "mhass.pbsim_runner.subprocess.run",
        side_effect=fake_subprocess_run,
    ):
        run_pbsim(
            template_dir=template_dir,
            output_base=output_dir,
            pbsim_path="/bin/true",
            qshmm_path="/dev/null",
            threads=threads,
            subread_accuracy=0.65,
            difference_ratio="6:55:39",
            master_seed=master_seed,
        )

    mapping = {}

    for cmd in pbsim_commands:
        template_path = cmd[
            cmd.index("--template") + 1
        ]

        template_name = Path(
            template_path
        ).name

        if "--seed" in cmd:
            seed = int(
                cmd[
                    cmd.index("--seed") + 1
                ]
            )
        else:
            seed = None

        mapping[template_name] = seed

    return mapping


def test_pbsim_seed_mapping_is_thread_invariant(
    tmp_path,
):
    template_dir = tmp_path / "templates"

    filenames = create_templates(
        template_dir
    )

    threads_1 = capture_seed_mapping(
        template_dir,
        tmp_path / "threads_1",
        threads=1,
        master_seed=MASTER_SEED,
    )

    threads_4 = capture_seed_mapping(
        template_dir,
        tmp_path / "threads_4",
        threads=4,
        master_seed=MASTER_SEED,
    )

    assert threads_1 == threads_4

    assert set(threads_1) == set(
        filenames
    )


def test_pbsim_seed_mapping_matches_direct_derivation(
    tmp_path,
):
    template_dir = tmp_path / "templates"

    filenames = create_templates(
        template_dir
    )

    actual = capture_seed_mapping(
        template_dir,
        tmp_path / "output",
        threads=4,
        master_seed=MASTER_SEED,
    )

    expected = {
        filename: derive_pbsim_seed(
            MASTER_SEED,
            filename,
        )
        for filename in filenames
    }

    assert actual == expected


def test_pbsim_different_master_seed_changes_mapping(
    tmp_path,
):
    template_dir = tmp_path / "templates"

    create_templates(
        template_dir
    )

    seed_12345 = capture_seed_mapping(
        template_dir,
        tmp_path / "seed_12345",
        threads=4,
        master_seed=12345,
    )

    seed_67890 = capture_seed_mapping(
        template_dir,
        tmp_path / "seed_67890",
        threads=4,
        master_seed=67890,
    )

    assert seed_12345 != seed_67890


def test_pbsim_unseeded_mode_omits_seed(
    tmp_path,
):
    template_dir = tmp_path / "templates"

    filenames = create_templates(
        template_dir
    )

    mapping = capture_seed_mapping(
        template_dir,
        tmp_path / "unseeded",
        threads=4,
        master_seed=None,
    )

    assert set(mapping) == set(
        filenames
    )

    assert all(
        seed is None
        for seed in mapping.values()
    )
