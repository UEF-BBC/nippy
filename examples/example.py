# A simple example of semi-automic preprocessing of NIR data using NIRPREP.
#
# jtorniainen // UEF
# 2018, MIT License

import nippy
import numpy as np


if __name__ == '__main__':
    # 1. Load configuration
    pipelines = nippy.read_configuration('example.ini')

    # 2. Load data
    data = np.genfromtxt('nir_data.csv', delimiter=',')
    wavelength = data[0, :]
    spectra = data[1:, :].T  # Rows = wavelength, Columns = samples

    # 3. Dataset through all pipelines
    datasets = nippy.nippy(wavelength, spectra, pipelines)

    # 4. Export the preprocessed data (showcasing three variants)
    nippy.export_pipelines_to_csv('/tmp/export-csv', datasets, pipelines, mkdir=True)
    nippy.export_pipelines_to_mat('/tmp/export-mat', datasets, pipelines, mkdir=True)
    nippy.export_pipelines_to_pickle('/tmp/export-pickle.p', datasets, pipelines)
