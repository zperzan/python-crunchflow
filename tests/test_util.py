import os
import shutil

import numpy as np
import pytest

from crunchflow.output import SpatialProfile
from crunchflow.util import correct_exponent


def test_correctexponent():
    shutil.copy("tests/data/wrr_floodplain_redox/volume740.tec", "tmp1.tec")
    correct_exponent("tmp1.tec", verbose="low")

    volume = SpatialProfile("tmp")
    data = volume.extract("Gypsum")
    os.remove("tmp1.tec")

    assert data.dtype == np.float64, "Incorrectly read data type"


if __name__ == "__main__":
    pytest.main()
