# Semi-automatic preprocessing script for NIR data. This script contains the preprocessing functions and some utility
# functions (like data export).
#
# jtorniainen, ioafara // Department of Applied Physics, University of Eastern Finland
# 2018, MIT License


import scipy.signal
import scipy.io as io
import scipy.ndimage as nd
import numpy as np
from sklearn.preprocessing import normalize, scale
from . import handler
import pickle
import os


class Preprocessor(object):
    """ Preprocessor object can be used to run nippy as an iterator (see documentation for examples). """

    def __init__(self, wavelength, spectra, configuration_file):
        """
        Args:
            wavelength <numpy.ndarray>: Vector of wavelengths.
            spectra <numpy.ndarray>: NIRS data matrix.
            configuration_file <str>: A path to the configuration file.
        """
        self.wavelength = wavelength
        self.spectra = spectra
        self.configuration = handler.read_configuration(configuration_file)
        self.current_pipe_idx = 0

    def __iter__(self):
        return self

    def __next__(self):
        """ Returns the next preprocessed dataset and a summary of preprocessing operations. """
        if self.current_pipe_idx >= len(self.configuration):
            raise StopIteration
        else:
            this_idx = self.current_pipe_idx
            wavelength_, spectra_ = run_pipeline(self.wavelength.copy(),
                                                 self.spectra.copy(),
                                                 self.configuration[this_idx])
            self.current_pipe_idx += 1
            return wavelength_, spectra_, self.configuration[this_idx]


# PREPROCESSING FUNCTIONS
def baseline(spectra):
    """ Removes baseline (mean) from each spectrum.

    Args:
        spectra <numpy.ndarray>: NIRS data matrix.

    Returns:
        spectra <numpy.ndarray>: Mean-centered NIRS data matrix
    """

    return spectra - np.mean(spectra, axis=0)

def snv(spectra):
    """ Perform scatter correction using the standard normal variate.

    Args:
        spectra <numpy.ndarray>: NIRS data matrix.

    Returns:
        spectra <numpy.ndarray>: NIRS data with (S/R)NV applied.
    """

    return (spectra - np.mean(spectra, axis=0)) / np.std(spectra, axis=0)


def rnv(spectra, iqr=[75, 25]):
    """ Perform scatter correction using robust normal variate.

    Args:
        spectra <numpy.ndarray>: NIRS data matrix.
        iqr <list>: IQR ranges [lower, upper] for robust normal variate.

    Returns:
        spectra <numpy.ndarray>: NIRS data with (S/R)NV applied.
    """

    return (spectra - np.median(spectra, axis=0)) / np.subtract(*np.percentile(spectra, iqr, axis=0))


def lsnv(spectra, num_windows=10):
    """ Perform local scatter correction using the standard normal variate.

    Args:
        spectra <numpy.ndarray>: NIRS data matrix.
        num_windows <int>: number of equispaced windows to use (window size (in points) is length / num_windows)

    Returns:
        spectra <numpy.ndarray>: NIRS data with local SNV applied.
    """

    parts = np.array_split(spectra, num_windows, axis=0)
    for idx, part in enumerate(parts):
        parts[idx] = snv(part)

    return np.concatenate(parts, axis=0)


def savgol(spectra, filter_win=11, poly_order=3, deriv_order=0, delta=1.0):
    """ Perform Savitzky–Golay filtering on the data (also calculates derivatives). This function is a wrapper for
    scipy.signal.savgol_filter.

    Args:
        spectra <numpy.ndarray>: NIRS data matrix.
        filter_win <int>: Size of the filter window in samples (default 11).
        poly_order <int>: Order of the polynomial estimation (default 3).
        deriv_order <int>: Order of the derivation (default 0).

    Returns:
        spectra <numpy.ndarray>: NIRS data smoothed with Savitzky-Golay filtering
    """
    return scipy.signal.savgol_filter(spectra, filter_win, poly_order, deriv_order, delta=delta, axis=0)


def trim(wavelength, spectra, bins):
    """ Trim spectra to a specified wavelength bin (or bins).

    Args:
        wavelength <numpy.ndarray>: Vector of wavelengths.
        spectra <numpy.ndarray>: NIRS data matrix.
        bins <list>: A bin or a list of bins defining the trim operation.

    Returns:
        spectra <numpy.ndarray>: NIRS data smoothed with Savitzky-Golay filtering
    """
    if type(bins[0]) != list:
        bins = [bins]

    spectra_trim = np.array([]).reshape(0, spectra.shape[1])
    wavelength_trim = np.array([])
    for wave_range in bins:
        mask = np.bitwise_and(wavelength >= wave_range[0], wavelength <= wave_range[1])
        spectra_trim = np.vstack((spectra_trim, spectra[mask, :]))
        wavelength_trim = np.hstack((wavelength_trim, wavelength[mask]))
    return wavelength_trim, spectra_trim


def resample(wavelength, spectra, resampling_ratio):
    """ Resample spectra according to the resampling ratio.

    Args:
        wavelength <numpy.ndarray>: Vector of wavelengths.
        spectra <numpy.ndarray>: NIRS data matrix.
        resampling_ratio <float>: new length with respect to original length

    Returns:
        wavelength_ <numpy.ndarray>: Resampled wavelengths.
        spectra_ <numpy.ndarray>: Resampled NIR spectra
    """

    new_length = int(np.round(wavelength.size * resampling_ratio))
    spectra_, wavelength_ = scipy.signal.resample(spectra, new_length, wavelength)
    return wavelength_, spectra_


def norml(spectra, udefined=True, imin=0, imax=1):
    """ Perform spectral normalisation with user-defined limits.

    Args:
        spectra <numpy.ndarray>: NIRS data matrix.
        udefined <bool>: use user defined limits
        imin <float>: user defined minimum
        imax <float>: user defined maximum

    Returns:
        spectra <numpy.ndarray>: Normalized NIR spectra
    """
    if udefined:
        f = (imax - imin)/(np.max(spectra) - np.min(spectra))
        n = spectra.shape
        arr = np.empty((0, n[0]), dtype=float) #create empty array for spectra
        for i in range(0, n[1]):
            d = spectra[:,i]
            dnorm = imin + f*d
            arr = np.append(arr, [dnorm], axis=0)
        return np.transpose(arr)
    else:
        return spectra / np.linalg.norm(spectra, axis=0)


def detrend(spectra, bp=0):
    """ Perform spectral detrending to remove linear trend from data.

    Args:
        spectra <numpy.ndarray>: NIRS data matrix.
        bp <list>: A sequence of break points. If given, an individual linear fit is performed for each part of data
        between two break points. Break points are specified as indices into data.

    Returns:
        spectra <numpy.ndarray>: Detrended NIR spectra
    """
    return scipy.signal.detrend(spectra, bp=bp)


def msc(spectra):
    """ Performs multiplicative scatter correction to the mean.

    Args:
        spectra <numpy.ndarray>: NIRS data matrix.

    Returns:
        spectra <numpy.ndarray>: Scatter corrected NIR spectra.
    """

    spectra = scale(spectra, with_std=False, axis=0) # Demean
    reference = np.mean(spectra, axis=1)

    for col in range(spectra.shape[1]):
        a, b = np.polyfit(reference, spectra[:, col], deg=1)
        spectra[:, col] = (spectra[:, col] - b) / a

    return spectra


def clip(wavelength, spectra, threshold, substitute=None):
    """ Removes or substitutes values above the given threshold.

    Args:
        wavelength <numpy.ndarray>: Vector of wavelengths.
        spectra <numpy.ndarray>: NIRS data matrix.
        threshold <float>: threshold value for rejection
        substitute <float>: substitute value for rejected values (None removes values from the spectra)

    Returns:
        wavelength <numpy.ndarray>: Vector of wavelengths.
        spectra <numpy.ndarray>: NIR spectra with threshold exceeding values removed.
    """

    if substitute == None:  # remove threshold violations
        mask = np.any(spectra > threshold, axis=1)
        spectra = spectra[~mask, :]
        wavelength = wavelength[~mask]
    else:  # substitute threshold violations with a value
        spectra[spectra > threshold] = substitute
    return wavelength, spectra

    return wavelength, spectra


def smooth(spectra, filter_win, window_type='flat', mode='reflect'):
    """ Smooths the spectra using convolution.

    Args:
        spectra <numpy.ndarray>: NIRS data matrix.
        filter_win <float>: length of the filter window in samples.
        window_type <str>: filtering window to use for convolution (see scipy.signal.windows)
        mode <str>: convolution mode

    Returns:
        spectra <numpy.ndarray>: Smoothed NIR spectra.
    """

    if window_type == 'flat':
        window = np.ones(filter_win)
    else:
        window = scipy.signal.windows.get_window(window_type, filter_win)
    window = window / np.sum(window)

    for column in range(spectra.shape[1]):
        spectra[:, column] = nd.convolve(spectra[:, column], window, mode=mode)

    return spectra


def derivate(spectra, order=1, delta=1):
    """ Computes Nth order derivates with the desired spacing using numpy.gradient.
    Args:
        spectra <numpy.ndarray>: NIRS data matrix.
        order <float>: Order of the derivation.
        delta <int>: Delta of the derivate (in samples).

    Returns:
        spectra <numpy.ndarray>: Derivated NIR spectra.
    """
    for n in range(order):
        spectra = np.gradient(spectra, delta, axis=0)
    return spectra


# UTILITY FUNCTIONS
def export_pipelines_to_csv(output_path, datasets, pipelines, mkdir=False):
    """ Exports all datasets and the related pipelines to csv files.

    Args:
        filename <str> output directory.
        datasets <list> list of datasets processed by nippy.
        pipelines <list> list of nippy pipelines.
        mkdir <bool> create output directory if it does not exist.
    """

    if mkdir and not os.path.isdir(output_path):
        os.mkdir(output_path)

    for idx, dataset in enumerate(datasets):
        filename = os.path.join(output_path, '{}.csv'.format(idx + 1))
        np.savetxt(filename, np.hstack((dataset[0].reshape(-1, 1), dataset[1])), delimiter=',')

    with open(os.path.join(output_path, 'pipelines.log'), 'w') as f:
        for idx, pipe in enumerate(pipelines):
            f.write('{};{}\n'.format(idx + 1, str(pipe)))


def export_pipelines_to_mat(output_path, datasets, pipelines, mkdir=False):
    """ Exports all datasets and the related pipelines to csv files.

    Args:
        filename <str> output directory.
        datasets <list> list of datasets processed by nippy.
        pipelines <list> list of nippy pipelines.
        mkdir <bool> create output directory if it does not exist.
    """

    if mkdir and not os.path.isdir(output_path):
        os.mkdir(output_path)

    new_datasets = []
    for idx, data, pipe in zip(range(len(datasets)), datasets, pipelines):
        dataset = {'data': data[1], 'wave': data[0], 'params': str(pipe)}
        io.savemat(os.path.join(output_path, '{}.mat'.format(idx + 1)), dataset)

    with open(os.path.join(output_path, 'pipelines.log'), 'w') as f:
        for idx, pipe in enumerate(pipelines):
            f.write('{};{}\n'.format(idx + 1, str(pipe)))


def export_pipelines_to_pickle(filename, datasets, pipelines):
    """ Exports all datasets and the related pipelines to a pickle file.

    Args:
        filename <str> output filepath.
        datasets <list> list of datasets processed by nippy.
        pipelines <list> list of nippy pipelines.
    """
    data  = {'datasets': datasets, 'pipelines': pipelines}
    pickle.dump(data, open(filename, 'wb'))


def run_pipeline(wavelength_, spectra_, pipeline):

        if 'CLIP' in pipeline.keys() and pipeline['CLIP'] != None:
            wavelength_, spectra_ = clip(wavelength_, spectra_, **pipeline['CLIP'])

        if 'BASELINE' in pipeline.keys() and pipeline['BASELINE'] != None:
            spectra_ = baseline(spectra_, **pipeline['BASELINE'])

        if 'SNV' in pipeline.keys() and pipeline['SNV'] != None:
            spectra_ = snv(spectra_, **pipeline['SNV'])

        if 'RNV' in pipeline.keys() and pipeline['RNV'] != None:
            spectra_ = rnv(spectra_, **pipeline['RNV'])

        if 'LSNV' in pipeline.keys() and pipeline['LSNV'] != None:
            spectra_ = lsnv(spectra_, **pipeline['LSNV'])

        if 'MSC' in pipeline.keys() and pipeline['MSC'] != None:
            spectra_ = msc(spectra_)

        if 'NORML' in pipeline.keys() and pipeline['NORML'] != None:
            spectra_ = norml(spectra_, **pipeline['NORML'])

        if 'SAVGOL' in pipeline.keys() and pipeline['SAVGOL'] != None:
            spectra_ = savgol(spectra_, **pipeline['SAVGOL'])

        if 'SMOOTH' in pipeline.keys() and pipeline['SMOOTH'] != None:
            spectra_ = smooth(spectra_, **pipeline['SMOOTH'])

        if 'DERIVATE' in pipeline.keys() and pipeline['DERIVATE'] != None:
            wavelength_, spectra_ = derivate(spectra_, **pipeline['DERIVATE'])

        if 'DETREND' in pipeline.keys() and pipeline['DETREND'] != None:
            spectra_ = detrend(spectra_, **pipeline['DETREND'])

        if 'RESAMPLE' in pipeline.keys() and pipeline['RESAMPLE'] != None:
            wavelength_, spectra_ = resample(wavelength_, spectra_, **pipeline['RESAMPLE'])

        if 'TRIM' in pipeline.keys() and pipeline['TRIM'] != None:
            wavelength_, spectra_ = trim(wavelength_, spectra_, **pipeline['TRIM'])


        return wavelength_, spectra_


def nippy(wavelength, spectra, pipelines):
    """ Main processing script of nippy. Applies operations specified in the 'pipelines' parameter to the given spectra.

    Args:
        wavelength <numpy.ndarray>: Vector of wavelengths.
        spectra <numpy.ndarray>: NIRS data matrix.
        pipelines <list>: list of nippy pipelines.

    Returns:
        datasets <list>: a list containing different preprocessed versions of the original spectra and wavelength.
    """

    datasets = []
    for idx, pipeline in enumerate(pipelines):
        wavelength_, spectra_ = run_pipeline(wavelength.copy(), spectra.copy(), pipeline)
        print('Running pipe {}:\n{}\n'.format(idx + 1, pipeline))
        datasets.append((wavelength_, spectra_))
    return datasets
