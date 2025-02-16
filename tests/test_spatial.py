import numpy as np
import pytest

from crunchflow.output.spatial import SpatialProfile


def test_cationexchange():
    totcon = SpatialProfile("totcon", folder="tests/data/cation_exchange")
    exchange = SpatialProfile("exchange", folder="tests/data/cation_exchange")
    correct_coords = np.load("tests/data/cation_exchange/correct_coords.npy")

    # Extract all exchange data into single array
    ex_arr = np.empty((len(totcon.times), 3, len(totcon.coords)))
    for time in exchange.times:
        for i, species in enumerate(["NaX", "CaX2", "CsX"]):
            ex_arr[time - 1, i, :] = np.log10(exchange.extract(species, time=time))
    correct_ex_arr = np.load("tests/data/cation_exchange/correct_exchange.npy")

    # Same with totcon
    conc_arr = np.empty((len(totcon.times), 3, len(totcon.coords)))
    for time in totcon.times:
        for i, species in enumerate(["Cs+", "Na+", "Tracer"]):
            conc_arr[time - 1, i, :] = np.log10(totcon.extract(species, time=time))
    correct_conc_arr = np.load("tests/data/cation_exchange/correct_totcon.npy")

    assert totcon.fmt == ".out", "Incorrectly read totcon file format"
    assert len(totcon.times) == 4, "Did not read all totcon files correctly"
    assert len(exchange.times) == 4, "Did not read all exchange files correctly"
    assert totcon.output_times == [150.0, 300.0, 450.0, 600.0], "Incorrectly read totcon output times"
    assert exchange.output_times == [150.0, 300.0, 450.0, 600.0], "Incorrectly read exchange output times"
    assert np.allclose(totcon.coords, correct_coords, atol=1e-3), "Incorrectly read totcon coordinates"
    assert np.allclose(ex_arr, correct_ex_arr, atol=1e-3), "Incorrectly read exchange data"
    assert np.allclose(conc_arr, correct_conc_arr, atol=1e-3), "Incorrectly read totcon data"


def test_wrr_floodplain():
    volume = SpatialProfile("volume", folder="tests/data/wrr_floodplain_redox")
    perm = SpatialProfile("permeability", folder="tests/data/wrr_floodplain_redox")
    correct_coords = np.load("tests/data/wrr_floodplain_redox/correct_coords.npy")
    correct_gridded_x = np.load("tests/data/wrr_floodplain_redox/correct_griddedX.npy")
    domain_shape = (volume.nx, volume.ny, volume.nz)
    correct_columns = [
        "Gypsum",
        "Calcite",
        "Mackinawite",
        "Ferrihydrite",
        "Fe(OH)3_HS",
        "Siderite",
        "S(s)",
        "C6H10O5(s)",
    ]
    calcite = volume.extract("Calcite")
    correct_calcite = np.load("tests/data/wrr_floodplain_redox/correct_calcite.npy")
    # compare permeability data
    correct_perm = np.load("tests/data/wrr_floodplain_redox/correct_perm.npy")
    perm_arr = np.empty((len(perm.columns), perm.ny, perm.nx))
    for i, col in enumerate(perm.columns):
        perm_arr[i, :, :] = perm.extract(col)

    assert domain_shape == (100, 20, 1), "Incorrectly read nx, ny, nz from volume file"
    assert (perm.nx, perm.ny, perm.nz) == (100, 20, 1), "Incorrectly read nx, ny, nz from permeability file"
    assert volume.griddedX.T.shape == tuple(n for n in domain_shape if n != 1), "Incorrectly gridded coordinates"
    assert sorted(volume.columns) == sorted(correct_columns), "Incorrectly read volume columns"
    assert volume.title == "Mineral Volumes (m^3 mineral/m^3 porous medium)", "Incorrectly read volume title"
    assert volume.times == [740], "Incorrectly read volume times"
    assert np.allclose(volume.griddedX, correct_gridded_x, atol=1e-3), "Incorrectly gridded coordinates"
    assert np.allclose(volume.coords, correct_coords, atol=1e-3), "Incorrectly read volume coordinates"
    assert np.allclose(calcite, correct_calcite, atol=1e-3), "Incorrectly read calcite data"
    assert np.allclose(perm_arr, correct_perm, atol=1e-3), "Incorrectly read permeability data"


def test_flow2d():
    rate = SpatialProfile("rate", folder="tests/data/flow2d")
    totmineral = SpatialProfile("TotMineral", folder="tests/data/flow2d")
    correct_columns = ["Calcite", "Gypsum", "Ferrihydrite", "Jarosite", "Gibbsite", "Siderite", "TracerMineral"]
    correct_oxygen = np.load("tests/data/flow2d/correct_oxygen.npy")
    oxygen = np.empty((len(totmineral.times), totmineral.ny, totmineral.nx), dtype=float)
    for i, time in enumerate(totmineral.times):
        oxygen[i] = totmineral.extract("O2(aq)", time)

    assert rate.fmt == ".tec", "Could not determine rate file format"
    assert rate.times == [1, 2, 3, 4, 5], "Incorrectly read number of rate files"
    assert totmineral.times == [1, 2, 3, 4, 5], "Incorrectly read number of TotMineral files"
    assert (rate.nx, rate.ny, rate.nz) == (31, 41, 1), "Incorrectly read nx, ny, nz from rate file"
    assert rate.griddedX.shape == (41, 31), "Incorrectly gridded coordinates"
    assert sorted(rate.columns) == sorted(correct_columns), "Incorrectly read rate columns"
    assert totmineral.title == "Total Component Concentration in Minerals (mol/m^3 PM)", (
        "Incorrectly read TotMineral title"
    )
    assert np.allclose(np.log10(oxygen), np.log10(correct_oxygen), atol=1e-3), "Incorrectly read TotMineral data"
    for col in rate.columns:
        data = rate.extract(col)
        assert data.dtype == np.float64, "Incorrect data type in rate file"


def test_phread():
    ph = SpatialProfile("pH", folder="tests/data/output_files")

    assert ph.fmt == ".out", "Could not determine pH file format"
    assert ph.times == [1, 2, 3, 4, 5], "Incorrectly read number of pH files"
    assert ph.columns == ["Distance", "pH"], "Incorrectly read columns from pH files"


def test_readdat():
    velocity = SpatialProfile("velocity", folder="tests/data/output_files")
    correct_coords = np.load("tests/data/output_files/correct_coords.npy")
    correct_velocity = np.load("tests/data/output_files/correct_velocity.npy")
    correct_columns = ["X Velocity", "Y Velocity"]

    assert velocity.fmt == ".dat", "Could not determine velocity1.dat file format"
    assert velocity.nx == 100, "Incorrectly read nx from velocity1.dat file"
    assert velocity.ny == 60, "Incorrectly read ny from velocity1.dat file"
    assert velocity.times == [1], "Incorrectly read number of velocity1.dat files"
    assert velocity.columns == correct_columns, "Incorrectly read columns from velocity1.dat file"
    assert np.allclose(velocity.coords, correct_coords, atol=1e-3), "Incorrectly read velocity1.dat coordinates"
    assert np.allclose(velocity.extract("X Velocity", 1), correct_velocity, atol=1e-3), (
        "Incorrectly read velocity1.dat data"
    )


if __name__ == "__main__":
    pytest.main()
