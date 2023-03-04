# A simple example of how to use the iterator functionality of nippy.
#
# jtorniainen // UEF
# 2019, MIT License

import nippy
import numpy as np


if __name__ == '__main__':

    # 1. Load data
    data = np.genfromtxt('intern.csv', delimiter=',')
    wavelength = data[0, :]
    spectra = data[1:, :].T  # Rows = wavelength, Columns = samples

    # 2. Initialize iterator
    iterator = nippy.Preprocessor(wavelength, spectra, 'example.ini')

    for wavelength_, spectra_, pipe_ in iterator:
        print(pipe_)


