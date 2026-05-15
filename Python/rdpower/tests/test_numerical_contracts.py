import contextlib
import io
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from rdpower import rdmde, rdpower, rdsampsi


def senate_data():
    data = pd.read_csv(Path(__file__).resolve().parents[2] / "rdpower_senate.csv")
    return data[["demvoteshfor2", "demmv"]]


def full_senate_data():
    return pd.read_csv(Path(__file__).resolve().parents[2] / "rdpower_senate.csv")


def quiet_call(fn, *args, **kwargs):
    stream = io.StringIO()
    with contextlib.redirect_stdout(stream):
        return fn(*args, **kwargs)


def test_rdpower_default_matches_python_baseline():
    result = quiet_call(rdpower, data=senate_data(), tau=5)

    assert result["N_l"] == 595
    assert result["N_r"] == 702
    assert result["Nh_l"] == 360
    assert result["Nh_r"] == 323
    np.testing.assert_allclose(result["power_rbc"], 0.8189905549293113, rtol=1e-10, atol=1e-10)
    np.testing.assert_allclose(result["se_rbc"], 1.7412585395308662, rtol=1e-10, atol=1e-10)
    np.testing.assert_allclose(result["samph_l"], 17.754395590143446, rtol=1e-10, atol=1e-10)


def test_rdpower_accepts_numpy_arrays():
    result = quiet_call(rdpower, data=senate_data().to_numpy(), tau=5)

    assert result["N_l"] == 595
    assert result["N_r"] == 702
    np.testing.assert_allclose(result["power_rbc"], 0.8189905549293113, rtol=1e-10, atol=1e-10)


def test_rdpower_forwards_current_rdrobust_options():
    data = full_senate_data()
    yx = data[["demvoteshfor2", "demmv"]]
    subset = data["demmv"].abs().to_numpy() <= 40
    subset_mask = pd.Series(subset, index=data.index)

    result = quiet_call(
        rdpower,
        data=yx,
        tau=5,
        subset=subset,
        cluster=data["state"].to_numpy(),
        vce="cr2",
        nnmatch=4,
    )

    assert result["N_l"] == data.loc[subset_mask & (data["demmv"] < 0), "state"].nunique()
    assert result["N_r"] == data.loc[subset_mask & (data["demmv"] >= 0), "state"].nunique()


def test_removed_python_r_vce_aliases_are_rejected():
    with pytest.raises(ValueError, match="no longer RDROBUST"):
        quiet_call(rdpower, data=senate_data(), tau=5, vce="cluster")


def test_rdsampsi_default_matches_python_baseline():
    result = quiet_call(rdsampsi, data=senate_data(), tau=5)

    assert result["N.l"] == 595
    assert result["N.r"] == 702
    assert result["Nh.l"] == 360
    assert result["Nh.r"] == 323
    np.testing.assert_allclose(result["sampsi.h.tot"], 657.0, rtol=1e-10, atol=1e-10)
    np.testing.assert_allclose(result["sampsi.h.l"], 366.0, rtol=1e-10, atol=1e-10)
    np.testing.assert_allclose(result["sampsi.h.r"], 291.0, rtol=1e-10, atol=1e-10)


def test_rdsampsi_returns_left_bias():
    result = quiet_call(
        rdsampsi,
        data=None,
        tau=1,
        samph=(1, 1),
        nsamples=(100, 50, 100, 50),
        bias=(2, 3),
        variance=(1, 1),
        init_cond=100,
    )

    assert result["bias.l"] == 2
    assert result["bias.r"] == 3


def test_rdmde_default_matches_python_baseline():
    result = quiet_call(rdmde, data=senate_data())

    assert result["N_l"] == 595
    assert result["N_r"] == 702
    assert result["Nh_l"] == 360
    assert result["Nh_r"] == 323
    np.testing.assert_allclose(result["mde"], 4.878278210831539, rtol=1e-10, atol=1e-10)
    np.testing.assert_allclose(result["se_rbc"], 62.709486906140825, rtol=1e-10, atol=1e-10)
