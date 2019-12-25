# CONFIGURATION: How to write the configuration file?

_(Last updated: 20.11.2019)_

Here you can find the syntax and a list of parameters for the different preprocessing operations currently included in `nippy`. Additional documentation can be found in the source code comments and docstrings.

Configuration files in `nippy` use the venerable [INI-format](https://en.wikipedia.org/wiki/INI_file). Each preprocessing option is entered into the configuration file as a section and the iterable parameters for that method are given as key-value pairs. Description of each preprocessing method and related parameters are given below.

**Note:** In addition to method specific parameters, each preprocessing steps accepts `also_skip` parameter. `also_skip` is a boolean which generates an additional pipeline where the method in question is left out.

**Note:** Currently, the order in which the operations are carried out is static. For instance, if clipping is part of the pipeline being tested it will always performed as the first operation, followed by scatter correction and so on. We might make the order of operations a configurable parameter in the future. For now, however, if you want to change the order, you can do so by modifying the `run_pipeline` function of `nippy.py`:

```python
def run_pipeline(wavelength_, spectra_, pipeline):  # CLIP performed always before SNV/RNV

        if 'CLIP' in pipeline.keys() and pipeline['CLIP'] != None:
            wavelength_, spectra_ = clip(wavelength_, spectra_, **pipeline['CLIP'])

        if 'SNV' in pipeline.keys() and pipeline['SNV'] != None:
            spectra_ = snv(spectra_, **pipeline['SNV'])

        if 'RNV' in pipeline.keys() and pipeline['RNV'] != None:
            spectra_ = rnv(spectra_, **pipeline['RNV'])
 
 
 ....
 
 
def run_pipeline(wavelength_, spectra_, pipeline):  # SNV/RNV performed always before CLIP
        
        if 'SNV' in pipeline.keys() and pipeline['SNV'] != None:
            spectra_ = snv(spectra_, **pipeline['SNV'])

        if 'RNV' in pipeline.keys() and pipeline['RNV'] != None:
            spectra_ = rnv(spectra_, **pipeline['RNV'])
        
        if 'CLIP' in pipeline.keys() and pipeline['CLIP'] != None:
            wavelength_, spectra_ = clip(wavelength_, spectra_, **pipeline['CLIP'])

        
```

## CLIP
Clipping removes values from the spectra that exceed the user given threshold. 

Available parameters are:

- `threshold`: Values in the spectra that exceed this parameter value are either rejected or subsituted depending on the value of the `substitute` parameter.
- `substitute`: This parameter determines what is done to the spectral values that exceed the `threshold` paramater. If substitute is `None`, the values are simply removed from the spectrum. If the parameter is a float, the exceeding values are subsituted with this value.

Example: Single pipeline where spurious spectral values exceeding 1e10 are removed and subsituted with zeroes.
```ini
[CLIP]
    threshold: 1e10
    substitute: 0
```

## TRIM
Removes unwanted wavelength ranges form the spectrum. The kept wavelength range can consist of one or several continuous ranges.

Available parameters are:

- `bins`: Determines the wavelength region (or regions) to use in the analysis. Single continuous wavelength region is specified by giving the starting and ending wavelength separated by a hyphen (`-`). The total wavelength region can consist of multiple continuous regions, different regions are separted with a comma (`,`). Finally, the delimiter for different pipelines is semicolon (`;`).

Example: single pipeline which removes everything outside of a specified wavelength range (1000-1900nm)
```ini
[TRIM]
    bins = 1000-1900
```
Example: single pipeline with a two wavelength ranges (1000-1500nm and 1900-2200nm)
```ini
[TRIM]
    bins = 1000-1500, 1900-2200
```
Example: two pipelines. Pipeline 1 has a single wavelength range (700-1800nm) and pipeline 2  has two wavelength ranges (1000-1500nm and 1900-2200nm)
```ini
[TRIM]
    bins = 700-1800; 1000-1500, 1900-2200
```

## SAVGOL
Performs Savitzky-Golay [Savitzky-Golay](https://en.wikipedia.org/wiki/Savitzky%E2%80%93Golay_filter) filtering on the spectrum.

Available parameters are:

- `filter_win`: Length of the filtering window in samples.
- `poly_order`: Order of the polynomial approximation in the filtering process. Defaults to 3.
- `deriv_order`: Savitzky-Golay filtering can also return the smooth derivate of the spectrum. This parameter determines the order of the derivate. Defaults to 0 (i.e., no derivate).
- `delta`: Spacing used when calculating derivates. Defaults to 1.

Example: Produces six pipelines with two filter window sizes (31 and 61) performed for the original, 1st derivative and 2nd derivative of the spectrum.

```ini
[SAVGOL]
    filter_win = 31, 61
    poly_order = 3
    deriv_order = 0, 1, 2
```

## SMOOTH
Smoothing spectrum via convolution according to a specific filter window length and windowing function.

Available parameters are:

- `filter_win`: Length of the filtering window in samples.
- `window_type`: Windowing function to use in the filtering (see for options see `scipy.signal.windows`). Defaults to `flat`.
- `mode`: How to treat the edges of the spectrum during filtering. Defaults to `'reflect'`

Example: Three pipelines, each with a different filtering window length (31, 61, and 91).

```ini
[SMOOTH]
    filter_win = 31, 61, 91
    window_type = hamming
```

## SNV
Performs standard (or robust) normal variate on the spectrum.

Available parameters are:

- `snv_type`: Type of normal variate to use. Options are `snv` for _standard normal variate_ and `rnv` for _robust normal variate_.
- `iqr`: In the case of `rnv`, this parameter can be used to define the inter-quartile range used for normalization. Defaults to `[75, 25]`.

Example: Three pipelines, two with scatter correction (standard and robust) and one without (by utilizing the `also_skip` option).

```ini
[SNV]
    snv_type = snv, rnv
    also_skip = True
```

## LSNV
Performs local standard (or robust) normal variate on the spectrum. Instead of performing the opration over the entire spectrum, the local (S/R)NV variant performs scatter correction in equispaced non-overlapping windows (as specified by the `num_windows` parameter).

Available parameters are:

- `snv_type`: Type of normal variate to use. Options are `snv` for _standard normal variate_ and `rnv` for _robust normal variate_.
- `num_windows`: Number of windows to split the spectrum to (defaults to `10`)
- `iqr`: In the case of `rnv`, this parameter can be used to define the inter-quartile range used for normalization. Defaults to `[75, 25]`.

Example: Three pipelines, two with scatter correction (standard and robust) and one without (by utilizing the `also_skip` option).

```ini
[SNV]
    snv_type = snv, rnv
    num_windows = 3, 6
    also_skip = True

## MSC
Performs multiplicative scatter correction to the mean of spectrum.

Available parameters are:

- No configurable parameters

Example: Two pipelines, with and without MSC.

```ini
[MSC]
    also_skip = True
```

## DETREND
Performs linear detrending on the spectrum.

Available parameters are:

- `bp`: A sequence of break points. If given, an individual linear fit is performed for each part of data between two break points. Break points are specified as indices into data.

Example: Two pipelines, with and without detrending.

```ini
[DETREND]
    also_skip = True
```

## NORML
Performs normalization on the spectrum.

Available parameters are:
- `udefined`: use user defined limits
- `imin`: user defined minimum
- `imax`: user defined maximum

Example: Two pipelines, with and without normalization.

```ini
[NORML]
    also_skip = True
```
