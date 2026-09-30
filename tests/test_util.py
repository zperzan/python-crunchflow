import os
import shutil

import numpy as np
import pytest

from crunchflow.output import SpatialProfile
from crunchflow.util import clear_output, correct_exponent


def test_correctexponent():
    shutil.copy("tests/data/wrr_floodplain_redox/volume740.tec", "tmp1.tec")
    correct_exponent("tmp1.tec", verbose="low")

    volume = SpatialProfile("tmp")
    data = volume.extract("Gypsum")
    os.remove("tmp1.tec")

    assert data.dtype == np.float64, "Incorrectly read data type"


# Files that clear_output should delete, including a three-digit time step
REMOVABLE = [
    "conc1.out",
    "pH12.out",
    "volume740.tec",
    "velocityx3.tec",
    "permeability1.dat",
    "CrunchJunk2.out",
    "fort.123",
    "initial_condition_Richards.tec",
]

# Files that clear_output should leave alone: input files, time_series output,
# and names that look like output but carry no time-step index
PRESERVED = [
    "cation_exchange.in",
    "ObservationWell01.txt",
    "notes.txt",
    "conc.out",
    "velocityx.out",
    "myconc1.tec",
]


def _populate(folder, names):
    for name in names:
        (folder / name).write_text("")


def test_clearoutput(tmp_path):
    _populate(tmp_path, REMOVABLE + PRESERVED)

    deleted = clear_output(folder=str(tmp_path), verbose=False)

    assert sorted(os.path.basename(p) for p in deleted) == sorted(REMOVABLE)
    assert sorted(p.name for p in tmp_path.iterdir()) == sorted(PRESERVED)


def test_clearoutput_dryrun(tmp_path):
    _populate(tmp_path, REMOVABLE)

    deleted = clear_output(folder=str(tmp_path), dry_run=True, verbose=False)

    assert sorted(os.path.basename(p) for p in deleted) == sorted(REMOVABLE)
    assert sorted(p.name for p in tmp_path.iterdir()) == sorted(REMOVABLE)


def test_clearoutput_suffixes(tmp_path):
    _populate(tmp_path, ["conc1.out", "conc1.tec", "conc1.dat"])

    clear_output(folder=str(tmp_path), suffixes=[".tec"], verbose=False)

    assert sorted(p.name for p in tmp_path.iterdir()) == ["conc1.dat", "conc1.out"]


def test_clearoutput_pestcontrol(tmp_path):
    _populate(tmp_path, ["myrun.out", "other.out"])
    (tmp_path / "PestControl.ant").write_text("myrun.in\n")

    clear_output(folder=str(tmp_path), verbose=False)

    assert sorted(p.name for p in tmp_path.iterdir()) == ["PestControl.ant", "other.out"]


def test_clearoutput_skips_subfolders(tmp_path):
    subfolder = tmp_path / "run2"
    subfolder.mkdir()
    _populate(subfolder, ["conc1.tec"])

    assert clear_output(folder=str(tmp_path), verbose=False) == []
    assert (subfolder / "conc1.tec").exists()


def test_clearoutput_missing_folder(tmp_path):
    with pytest.raises(NotADirectoryError):
        clear_output(folder=str(tmp_path / "nope"))


if __name__ == "__main__":
    pytest.main()
