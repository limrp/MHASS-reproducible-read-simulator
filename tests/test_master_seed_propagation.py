from contextlib import ExitStack
from types import SimpleNamespace
from unittest.mock import patch

from mhass.main import run_simulation


def make_args(output_dir, seed):
    return SimpleNamespace(
        output_dir=str(output_dir),
        barcode_file=None,
        np_distribution=None,
        counts_file="dummy_counts.tsv",
        amplicon_fasta="dummy_amplicons.fa",
        np_distribution_type="lognormal",
        lognormal_mu=3.88,
        lognormal_sigma=1.22,
        np_min=2,
        np_max=59,
        threads=4,
        subread_accuracy=0.65,
        difference_ratio="6:55:39",
        seed=seed,
    )


def capture_propagated_seeds(
    output_dir,
    seed,
):
    args = make_args(
        output_dir,
        seed,
    )

    with ExitStack() as stack:
        stack.enter_context(
            patch(
                "mhass.main.prepare_external_counts"
            )
        )

        mock_templates = stack.enter_context(
            patch(
                "mhass.main.create_per_sequence_templates"
            )
        )

        mock_pbsim = stack.enter_context(
            patch(
                "mhass.main.run_pbsim"
            )
        )

        mock_combiner = stack.enter_context(
            patch(
                "mhass.main.combine_fastqs"
            )
        )

        success = run_simulation(args)

    mock_templates.assert_called_once()
    mock_pbsim.assert_called_once()
    mock_combiner.assert_called_once()

    return {
        "success": success,
        "template_seed": (
            mock_templates.call_args[1]["master_seed"]
        ),
        "pbsim_seed": (
            mock_pbsim.call_args[1]["master_seed"]
        ),
        "combiner_seed": (
            mock_combiner.call_args[1]["master_seed"]
        ),
    }


def test_master_seed_propagates_to_all_components(
    tmp_path,
):
    result = capture_propagated_seeds(
        tmp_path / "seeded",
        12345,
    )

    assert result["success"] is True

    assert result["template_seed"] == 12345
    assert result["pbsim_seed"] == 12345
    assert result["combiner_seed"] == 12345


def test_unseeded_mode_propagates_none(
    tmp_path,
):
    result = capture_propagated_seeds(
        tmp_path / "unseeded",
        None,
    )

    assert result["success"] is True

    assert result["template_seed"] is None
    assert result["pbsim_seed"] is None
    assert result["combiner_seed"] is None
