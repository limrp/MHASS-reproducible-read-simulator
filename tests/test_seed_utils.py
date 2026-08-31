from mhass.seed_utils import (
    PBSIM_SEED_MAX,
    derive_pbsim_seed,
    derive_seed,
)


def test_derive_seed_matches_known_value():
    seed = derive_seed(
        12345,
        "orientation",
    )

    assert seed == 17686612631035152410


def test_derive_seed_is_deterministic():
    first = derive_seed(
        12345,
        "orientation",
    )

    second = derive_seed(
        12345,
        "orientation",
    )

    assert first == second


def test_derive_seed_separates_components():
    orientation_seed = derive_seed(
        12345,
        "orientation",
    )

    empirical_np_seed = derive_seed(
        12345,
        "np-empirical",
    )

    assert orientation_seed != empirical_np_seed


def test_derive_seed_separates_identifiers():
    template1_seed = derive_seed(
        12345,
        "pbsim",
        "template1.fasta",
    )

    template2_seed = derive_seed(
        12345,
        "pbsim",
        "template2.fasta",
    )

    assert template1_seed != template2_seed


def test_derive_seed_changes_with_master_seed():
    seed_12345 = derive_seed(
        12345,
        "orientation",
    )

    seed_67890 = derive_seed(
        67890,
        "orientation",
    )

    assert seed_12345 != seed_67890


def test_derive_pbsim_seed_is_deterministic():
    first = derive_pbsim_seed(
        12345,
        "template1.fasta",
    )

    second = derive_pbsim_seed(
        12345,
        "template1.fasta",
    )

    assert first == second


def test_derive_pbsim_seed_stays_in_valid_range():
    master_seeds = (
        0,
        1,
        12345,
        67890,
        2**63,
    )

    template_identifiers = (
        "template1.fasta",
        "template2.fasta",
    )

    for master_seed in master_seeds:
        for template_identifier in template_identifiers:
            pbsim_seed = derive_pbsim_seed(
                master_seed,
                template_identifier,
            )

            assert 1 <= pbsim_seed <= PBSIM_SEED_MAX
