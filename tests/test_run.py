import os

import pytest

from crunchflow.input import Condition, InputFile


def test_add_condition():
    crunchrun = InputFile()
    crunchrun.set_block('condition', {'name': 'groundwater', 'units': 'mol/l', 'temperature': '25'})

    assert 'groundwater' in crunchrun.conditions
    assert crunchrun.conditions['groundwater'].units == 'mol/l'
    assert crunchrun.conditions['groundwater'].temperature == '25'

def test_get_condition_by_name():
    crunchrun = InputFile()
    crunchrun.set_block('condition', {'name': 'groundwater', 'units': 'mol/l', 'temperature': '25'})
    condition = crunchrun.conditions['groundwater']
    assert condition.units == 'mol/l'
    assert condition.temperature == '25'

def test_condition_str():
    condition = Condition('groundwater')
    condition.set_parameters({'units': 'mol/l', 'temperature': '25'})
    condition_str = str(condition)
    assert 'units                  mol/l' in condition_str
    assert 'temperature            25' in condition_str

def test_run_str():
    crunchrun = InputFile()
    crunchrun.set_block('condition', {'name': 'groundwater', 'units': 'mol/l', 'temperature': '25'})
    run_str = str(crunchrun)
    assert 'CONDITION              groundwater' in run_str
    assert 'units                  mol/l' in run_str
    assert 'temperature            25' in run_str

def test_save_and_load():
    crunchrun = InputFile()
    crunchrun.set_block('condition', {'name': 'groundwater', 'units': 'mol/l', 'temperature': '25'})
    crunchrun.save('test_output.in', path='.')

    loaded_run = InputFile.load('test_output.in', path='.')
    if os.path.exists('test_output.in'):
        os.remove('test_output.in')

    assert 'groundwater' in loaded_run.conditions
    assert loaded_run.conditions['groundwater'].units == 'mol/l'
    assert loaded_run.conditions['groundwater'].temperature == '25'

def test_load_cationexchange():
    crunchrun = InputFile.load('tests/data/cation_exchange.in')
    primary_species = crunchrun.primary_species.species
    correct_species = ['Cs+', 'Na+', 'Ca++', 'Tracer', 'H+', 'NO3-']
    conditions = list(crunchrun.conditions.keys())

    assert len(crunchrun.__dict__) == 22, 'Did not read all blocks in cation_exchange.in'
    assert sorted(primary_species) == sorted(correct_species), 'Did not read primary species correctly'
    assert conditions == ['flush', 'initial_condition']

def test_load_flow2d():
    crunchrun = InputFile.load('tests/data/flow2d.in')
    crunchrun.runtime.Benchmark = True
    run_str = str(crunchrun.runtime)

    assert len(crunchrun.secondary_species.species) == 41, 'Incorrect number of secondary species'
    assert crunchrun.flow.distance_units == 'meters', 'Incorrectly read distance_units in FLOW'
    assert 'Benchmark              true' in run_str, 'Unable to properly set Benchmark to true'

def test_load_speciation():
    crunchrun = InputFile.load('tests/data/speciation_exercise.in')

    conditions = list(crunchrun.conditions.keys())
    primary_species = crunchrun.primary_species.species
    cond_species = crunchrun.conditions['speciate1'].species
    correct_species = ['pH', 'CO2(aq)', '13CO2(aq)', 'Na+', 'Ca++', 'Cl-']
    na_conc = crunchrun.conditions['speciate1'].concentrations['Na+']

    assert conditions == ['speciate1', 'speciate2', 'speciate3'], 'Did not read conditions correctly'
    assert len(primary_species) == 6, 'Did not read primary species correctly'
    assert sorted(cond_species) == sorted(correct_species), 'Did not read species in speciate1 correctly'
    assert float(na_conc) == 0.01, 'Could not read Na+ concentration in speciate1'

if __name__ == "__main__":
    pytest.main()
