import pytest

from crunchflow.output.main_output import MainOutputFile


def test_short_course_ten():
    ex = MainOutputFile("tests/data/output_files/ShortCourse10.out")
    first_condition_name = list(ex.speciation_blocks.keys())[0]
    second_cond = ex.speciation_blocks["initial"]

    assert first_condition_name == "mingliang_boundary"

    assert second_cond.temperature == 25.0
    assert second_cond.porosity == 0.35
    assert second_cond.liquid_saturation == 1.0
    assert second_cond.liquid_density == 997.075
    assert second_cond.solid_density is None
    assert second_cond.solid_solution_ratio is None
    assert second_cond.ionic_strength == 0.003
    assert second_cond.pH == 7.0
    assert second_cond.pe == 10.751
    assert second_cond.eh == 0.636
    assert second_cond.total_charge == -7.007e-18


if __name__ == "__main__":
    pytest.main()
