![alt text](https://raw.githubusercontent.com/UEF-BBC/nippy/master/nippy.png?token=AIFYREKKYl0silMboodhUkS4orBUeJJLks5b5AcHwA%3D%3D "Semi-automic NIRS preprocessor")
# nippy
Semi-automated preprocessing Python module for near infrared spectroscopic (NIRS) data.

## Introduction
`nippy` is a Python (3.6+) module for rapid exploration of different NIRS preprocessing methods. `nippy` collects and wraps the most common preprocessing methods and provides tools for quickly constructing preprocessing pipes with alternative preprocessing combinations. Aim of this module is to enable the user to quickly test multiple alternativ preprocessing techniques and test how that affects the performance of the NIRS model.

## Usage
Comprehensive manual is still being worked on. For a simplified example of how `nippy` works you can look into the _examples_ directory. We provide here a crash-course into how `nippy` can be used.

The typical structure of the `nippy` analysis is as follows:

1. Specify the methods you wish to try and the associated parameters by generating an INI-formatted configuration file.
(for more detailed documentation about writing configuration files please check out the [CONFIGURATION.md](CONFIGURATION.md)). For example, configuring `nippy` to test 2nd derivative Savitzky-Golay filtering (with 3rd order polynomial fit) at three different filter-lengths (7, 11 and 31 samples) can be accomplished by adding the following section to the configuration file.

```ini
[SAVGOL]
    filter_win = 7, 11, 31
    poly_order = 3
    deriv_order = 2
    also_skip = True
```

2. Load your NIR data into a `numpy` matrix (rows wavelengths, columns samples). Load your wavelengths into a numpy vector.

```python
    data = np.genfromtxt('nir_data.csv', delimiter=',')
    wavelength = data[0, :]
    spectra = data[1:, :].T  # Rows = wavelength, Columns = samples
```

3. Import `nippy` and read your protocol file using `nippy.read_configuration`.

```python
    import nippy
    pipelines = nippy.read_pipeline('example_protocol.ini')
```

4. `nippy` generates a list of all possible preprocessing permutation. Pass your data and the list of pipelines to the `nippy`-function.

```python
    datasets = nippy.nippy(wavelength, spectra, pipelines)
```

The variable `datasets` now contains a list of datasets that have been preprocessed according to the methods listed in the `pipelines` variable. Preprocessed data can be used in Python or exported for use in other applications.


## Requirements

```
numpy (1.13.1+)
scipy (0.19.1+)
sklearn (0.19.2+)
```

## Installation
```
pip install git+https://github.com/uef-bbc/nippy
```

## Repository content
### nippy
- `nippy.py`: contains all of the preprocessing operations
- `handler.py`: top-level script for generating and running multiple preprocessing pipelines
### examples
- `example.py`: example script for performing multiple preprocessing pipelines
- `example.ini`: example `nippy` protocol
- `nir_data.csv`: small NIR dataset for demonstration purposes
