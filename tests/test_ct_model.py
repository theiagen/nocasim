from nocasim.ct_model import ct_to_viral_fraction, predict_completeness


def test_ct25_gives_base_fraction():
    vf = ct_to_viral_fraction(25.0)
    assert abs(vf - 0.03) < 1e-6


def test_lower_ct_higher_fraction():
    vf_low = ct_to_viral_fraction(20.0)
    vf_high = ct_to_viral_fraction(30.0)
    assert vf_low > vf_high


def test_ct28_approx_1_2_percent():
    vf = ct_to_viral_fraction(28.0)
    assert 0.005 < vf < 0.02


def test_ct20_high_titer():
    vf = ct_to_viral_fraction(20.0)
    assert 0.05 < vf < 0.20


def test_ct35_near_detection_limit():
    vf = ct_to_viral_fraction(35.0)
    assert vf < 0.005


def test_predict_completeness_low_ct():
    call = predict_completeness(20.0, total_reads=1_000_000)
    assert call == "complete"


def test_predict_completeness_high_ct():
    call = predict_completeness(40.0, total_reads=10_000)
    assert call in ("incomplete", "low_coverage")
