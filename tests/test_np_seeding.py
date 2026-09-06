from mhass.template_creator import create_np_sampler


def make_lognormal_params():
    return {
        "distribution_type": "lognormal",
        "lognormal_mu": 3.88,
        "lognormal_sigma": 1.22,
        "np_min": 2,
        "np_max": 59,
        "empirical_file": None,
    }


def sample_lognormal(seed, n=50):
    sampler = create_np_sampler(
        make_lognormal_params(),
        master_seed=seed,
    )

    return [sampler() for _ in range(n)]


def make_empirical_params(np_file):
    return {
        "distribution_type": "empirical",
        "empirical_file": np_file,
        "lognormal_mu": 3.88,
        "lognormal_sigma": 1.22,
        "np_min": 2,
        "np_max": 59,
    }


def sample_empirical(np_file, seed, n=50):
    sampler = create_np_sampler(
        make_empirical_params(np_file),
        master_seed=seed,
    )

    return [sampler() for _ in range(n)]


def write_empirical_distribution(tmp_path):
    np_file = tmp_path / "np_distribution.tsv"

    np_file.write_text(
        "2\t1\n"
        "5\t2\n"
        "9\t1\n"
    )

    return np_file


def test_lognormal_same_seed_reproduces_sequence():
    first = sample_lognormal(12345)
    second = sample_lognormal(12345)

    assert first == second


def test_lognormal_different_seed_changes_sequence():
    seed_12345 = sample_lognormal(12345)
    seed_67890 = sample_lognormal(67890)

    assert seed_12345 != seed_67890


def test_lognormal_samples_respect_bounds():
    values = sample_lognormal(12345)

    assert all(
        2 <= value <= 59
        for value in values
    )


def test_lognormal_unseeded_sampler_still_works():
    sampler = create_np_sampler(
        make_lognormal_params(),
        master_seed=None,
    )

    value = sampler()

    assert 2 <= value <= 59


def test_empirical_same_seed_reproduces_sequence(tmp_path):
    np_file = write_empirical_distribution(tmp_path)

    first = sample_empirical(
        np_file,
        12345,
    )

    second = sample_empirical(
        np_file,
        12345,
    )

    assert first == second


def test_empirical_different_seed_changes_sequence(tmp_path):
    np_file = write_empirical_distribution(tmp_path)

    seed_12345 = sample_empirical(
        np_file,
        12345,
    )

    seed_67890 = sample_empirical(
        np_file,
        67890,
    )

    assert seed_12345 != seed_67890


def test_empirical_samples_come_from_distribution(tmp_path):
    np_file = write_empirical_distribution(tmp_path)

    values = sample_empirical(
        np_file,
        12345,
    )

    allowed_values = {
        2,
        5,
        9,
    }

    assert all(
        value in allowed_values
        for value in values
    )


def test_empirical_unseeded_sampler_still_works(tmp_path):
    np_file = write_empirical_distribution(tmp_path)

    sampler = create_np_sampler(
        make_empirical_params(np_file),
        master_seed=None,
    )

    value = sampler()

    assert value in {
        2,
        5,
        9,
    }
