Here you can find the syntax and a list of parameters for the different preprocessing operations currently included in `nippy`. Additional documentation can be found in the source code docstrings.

Configuration files follow the INI-format. Each preprocessing option is entered into the configuration file as a section and the iterable parameters for that method are given as key-value pairs. Description of each preprocessing method and related parameters are given below.

**Note:** In addition to method specific parameters, each preprocessing steps accepts `also_skip` parameter. `also_skip` is a boolean which generates an additional pipeline where the method in question is left out.

## CLIP
Clipping removes values from the spectra that exceed the user given threshold. 

Available parameters are:

- `threshold`: Values in the spectra that exceed this parameter value are either rejected or subsituted depending on the value of the `substitute` parameter.
- `substitute`: This parameter determines what is done to the spectral values that exceed the `threshold` paramater. If substitute is `None`, the values are simply removed from the spectrum. If the parameter is a float, the exceeding values are subsituted with this value.

Example
```ini
[CLIP]
    threshold: 1e10
    substitute: 0
```

## TRIM
Removes unwanted wavelength ranges form the spectrum. 

Available parameters are:

- `bins`: Determines the wavelength region (or regions) to use in the analysis. Single continuous wavelength region is specified by giving the starting and ending wavelength separated by a hyphen (`-`). The total wavelength region can consist of multiple continuous regions, different regions are separted with a comma (`,`). Finally, the delimiter for different pipelines is semicolon (`;`).

Example with a single pipeline with a single wavelength range (1000-1900nm)
```ini
[TRIM]
    bins = 1000-1900
```
Example with a single pipeline with a two wavelength ranges (1000-1500nm and 1900-2200nm)
```ini
[TRIM]
    bins = 1000-1500, 1900-2200
```
Example with two pipelines. One with a single wavelength range (700-1800nm) and one with a two wavelength ranges (1000-1500nm and 1900-2200nm)
```ini
[TRIM]
    bins = 700-1800; 1000-1500, 1900-2200
```

## SAVGOL
Performs Savitzky-Golay filtering on the spectrum.

Available parameters are:

- `filter_win`: Length of the filtering window in samples.
- `poly_order`: Order of the polynomial approximation in the filtering process. Defaults to 3.
- `deriv_order`: Savitzky-Golay filtering can also return the smooth derivate of the spectrum. This parameter determines the order of the derivate. Defaults to 0 (i.e., no derivate).


## SMOOTH
Smooth filtering of the spectrum.

Available parameters are:

- `filter_win`: Length of the filtering window in samples.
- `window_type`: Windowing function to use in the filtering (see for options see `scipy.signal.windows`). Defaults to `flat`.
- `mode`: How to treat the edges of the spectrum during filtering. Defaults to `'reflect'`

## SNV
Performs standard (or robust) normal variate on the spectrum.

Available parameters are:

- `snv_type`: Type of normal variate to use. Options are `snv` for _standard normal variate_ and `rnv` for _robust normal variate_.
- `iqr`: In the case of `rnv`, this parameter can be used to define the inter-quartile range used for normalization. Defaults to `[75, 25]`.

## MSC
Performs multiplicative scatter correction to the mean of spectrum.

Available parameters are:

- No configurable parameters

## DETREND
Performs detrending on the spectrum.

Available parameters are:

## NORML
Performs normalization on the spectrum.

Available parameters are:
