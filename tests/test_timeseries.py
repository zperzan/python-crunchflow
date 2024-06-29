import os
import numpy as np
import pytest
from crunchflow.output import TimeSeries

def test_timeseries_load_well():
    ts = TimeSeries('Well01.txt', folder='tests/data/wrr_floodplain_redox')
    logged_conc = np.log10(ts.data[:, 1:])
    saved_conc = np.load('tests/data/wrr_floodplain_redox/correct_conc.npy')

    assert np.allclose(logged_conc, saved_conc, atol=1e-3), 'Logged concentrations do not match saved concentrations'
    assert len(ts.columns) == 48, 'Did not read in all columns'
    assert ts.coords == (100, 1, 1), 'Did not read in correct coordinates'
    assert ts.unit == 'mol/L', 'Did not read in unit correctly'
    assert ts.timeunit == 'day', 'Did not read in time unit correctly'

def test_timeseries_convert_mgL():
    ts = TimeSeries('Well01.txt', folder='tests/data/wrr_floodplain_redox')
    ts.convert_mgL(database='datacom.dbs', folder='tests/data/wrr_floodplain_redox',
                   warnings=False)
    logged_conc_mgl = np.log10(ts.data[:, 1:])
    saved_conc_mgl = np.load('tests/data/wrr_floodplain_redox/correct_conc_mgl.npy')

    assert ts.unit == 'mg/L', "Unit should be converted to mg/L"
    assert np.allclose(logged_conc_mgl, saved_conc_mgl, atol=1e-3), 'Logged concentrations do not match saved concentrations'

def test_timeseries_load_breakthrough():
    ts = TimeSeries('breakthrough.out', folder='tests/data/cation_exchange')
    logged_conc = np.log10(ts.data[:, 1:])
    saved_conc = np.load('tests/data/cation_exchange/correct_breakthrough.npy')
    correct_species = ['Tracer', 'Cs+', 'Ca++', 'Na+']
    read_species = ts.species

    assert sorted(read_species) == sorted(correct_species), 'Did not read in species correctly'
    assert np.allclose(logged_conc, saved_conc, atol=1e-3), 'Logged concentrations do not match saved concentrations'
    assert ts.coords == (30, 1, 1), 'Did not read in correct coordinates'
    assert ts.unit == 'mol/L', 'Did not read in unit correctly'
    assert ts.timeunit == 'yrs', 'Did not read in time unit correctly'

if __name__ == "__main__":
    pytest.main()
