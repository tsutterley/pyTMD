"""Synthetic extrapolation tests without external model downloads."""

import numpy as np
import pytest
import xarray as xr

import pyTMD  # noqa: F401 - register xarray accessors

pytest.importorskip("dask.array")
from dask.callbacks import Callback


def model():
    values = np.array([[1.0, 2.0, 3.0], [4.0, np.nan, 6.0], [7.0, 8.0, 9.0]])
    return xr.Dataset(
        {
            "m2": (("y", "x"), values.astype(complex) * (1 + 2j)),
            "s2": (("y", "x"), values * 2),
        },
        coords={"x": [0.0, 1.0, 2.0], "y": [0.0, 1.0, 2.0]},
    )


@pytest.mark.parametrize("k", [1, 3])
@pytest.mark.parametrize("cutoff", [np.inf, 150.0])
@pytest.mark.parametrize("scalar", [False, True])
def test_extrapolation_is_lazy(k, cutoff, scalar):
    ds = model()
    x = xr.DataArray([0.5, 1.0, 3.0, 10.0], dims="point")
    y = xr.DataArray([0.5, 1.0, 1.0, 10.0], dims="point")
    if scalar:
        x, y = x[2], y[2]
    target = ds.interp(x=x, y=y)
    expected = ds.tmd.extrap_like(target.copy(deep=True), k=k, cutoff=cutoff)
    lazy_target = target.copy(deep=True)
    for v in lazy_target:
        lazy_target[v].data = lazy_target[v].chunk().data
    tasks = []
    with Callback(pretask=lambda key, *args: tasks.append(key)):
        actual = ds.chunk({"x": 2, "y": 2}).tmd.extrap_like(
            lazy_target, k=k, cutoff=cutoff
        )
    assert not tasks
    assert all(hasattr(actual[v].data, "__dask_graph__") for v in actual)
    xr.testing.assert_identical(actual.compute(), expected)


def test_interp_extrapolate_is_lazy():
    ds = model()
    points = xr.DataArray([0.5, 3.0], dims="point")
    tasks = []
    with Callback(pretask=lambda key, *args: tasks.append(key)):
        actual = ds.chunk({"x": 2, "y": 2}).tmd.interp(
            points, points, extrapolate=True
        )
    assert not tasks
    expected = ds.tmd.interp(points, points, extrapolate=True)
    xr.testing.assert_identical(actual.compute(), expected)


def test_extrapolation_propagates_read_errors():
    import dask.array as da
    from dask import delayed

    @delayed
    def broken():
        raise RuntimeError("model read failed")

    ds = model()
    ds["m2"].data = da.from_delayed(broken(), shape=(3, 3), dtype=complex)
    target = model().interp(
        x=xr.DataArray([3.0], dims="point"), y=xr.DataArray([1.0], dims="point")
    )
    result = ds.tmd.extrap_like(target)
    with pytest.raises(RuntimeError, match="model read failed"):
        result.compute()


@pytest.mark.parametrize("missing", [False, True])
@pytest.mark.parametrize("crs", [4326, 3857])
def test_lazy_extrapolation_masks_and_metadata(missing, crs):
    ds = model()
    ds.attrs["crs"] = crs
    ds["m2"].attrs["units"] = "m"
    ds["s2"].data[0, 0] = np.nan
    x = xr.DataArray([0.5, 1.5], dims="point")
    target = ds.interp(x=x, y=x)
    if not missing:
        target = target.fillna(42.0)
    expected = ds.tmd.extrap_like(target.copy(deep=True), cutoff=1.0)
    actual = ds.chunk({"x": 1, "y": 1}).tmd.extrap_like(
        target.copy(deep=True), cutoff=1.0
    )
    xr.testing.assert_identical(actual.compute(), expected)


def test_interp_extrapolate_grid_is_lazy():
    ds = model()
    x, y = np.meshgrid([0.5, 1.0, 3.0], [0.5, 1.0, 3.0])
    x, y = ds.tmd.coords_as(x, y, type="grid")
    tasks = []
    with Callback(pretask=lambda key, *args: tasks.append(key)):
        actual = ds.chunk({"x": 2, "y": 2}).tmd.interp(
            x, y, extrapolate=True, cutoff=150.0, k=3
        )
    assert not tasks
    expected = ds.tmd.interp(x, y, extrapolate=True, cutoff=150.0, k=3)
    xr.testing.assert_identical(actual.compute(), expected)


def test_model_reads_are_shared_and_inputs_are_preserved():
    import dask.array as da
    from dask import delayed

    reads = []
    ds = model()
    values = ds["m2"].values.copy()

    @delayed
    def read():
        reads.append(True)
        return values

    ds["m2"].data = da.from_delayed(read(), shape=values.shape, dtype=complex)
    target = model().interp(
        x=xr.DataArray([3.0], dims="point"), y=xr.DataArray([1.0], dims="point")
    )
    original = target.copy(deep=True)
    actual = ds.tmd.extrap_like(target)
    assert not reads
    actual.compute(scheduler="synchronous")
    assert reads == [True]
    xr.testing.assert_identical(target, original)


def test_nearest_neighbor_values_and_cutoff():
    ds = model().chunk({"x": 2, "y": 2})
    x = xr.DataArray([3.0, 10.0], dims="point")
    y = xr.DataArray([1.0, 10.0], dims="point")
    result = ds.tmd.interp(x, y, extrapolate=True, cutoff=150.0).compute()
    np.testing.assert_allclose(result["m2"], [6.0 + 12.0j, np.nan])
    np.testing.assert_allclose(result["s2"], [12.0, np.nan])


def test_constituents_with_different_chunk_sizes():
    ds = model()
    ds["m2"] = ds["m2"].chunk({"x": 1, "y": 2})
    ds["s2"] = ds["s2"].chunk({"x": 2, "y": 1})
    points = xr.DataArray([0.5, 3.0], dims="point")
    actual = ds.tmd.interp(points, points, extrapolate=True)
    expected = model().tmd.interp(points, points, extrapolate=True)
    xr.testing.assert_identical(actual.compute(), expected)


def test_interp_extrapolation_reads_model_once():
    """Interpolation and extrapolation must share the model read graph."""
    import dask.array as da
    from dask import delayed

    reads = []
    ds = model()
    values = ds["m2"].values.copy()

    @delayed
    def read():
        reads.append(True)
        return values

    ds["m2"].data = da.from_delayed(read(), shape=values.shape, dtype=complex)
    points = xr.DataArray([0.5, 3.0], dims="point")
    actual = ds.tmd.interp(points, points, extrapolate=True)
    assert not reads
    actual.compute(scheduler="synchronous")
    assert reads == [True]


def test_shared_graph_restores_names_and_scalar_coordinates():
    ds = model().rename({"m2": "__extrap_0_0", "s2": "__extrap_1_0"})
    ds = ds.assign_coords(source_marker="source")
    ds.attrs["description"] = "source metadata"
    points = xr.DataArray([0.5, 3.0], dims="point")
    target = ds.interp(x=points, y=points).drop_vars("source_marker")
    target = target.assign_coords(target_marker="target")
    target.attrs["description"] = "target metadata"
    for v in target:
        target[v].attrs["units"] = "m"
        target[v].encoding["custom_metadata"] = "kept"
    expected = ds.tmd.extrap_like(target.copy(deep=True))
    actual = ds.chunk({"x": 2, "y": 2}).tmd.extrap_like(target)
    xr.testing.assert_identical(actual.compute(), expected)
    for v in target:
        assert actual[v].encoding == target[v].encoding
