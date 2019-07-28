# Handler function for generating different preprocessing combinations from configuration files.
#
# jtorniainen // Department of Applied Physics, University of Eastern Finland
# 2018, MIT License

import configparser
import itertools


# ------- SECTION-SPECIFIC PARSING FUNCTIONS --------
def parse_savgol(config):
    """ Parse arguments for Savitzky-Golay filtering.

    Args:
        config <dict>: dictionary of configuration options for savgol.
    Returns:
        config <dict>: dictionary of configuration options with parsed values.

    """
    for key in config:
        config[key] = _parse_list(config[key])

    return config


def parse_resample(config):
    """ Parse arguments for Savitzky-Golay filtering.

    Args:
        config <dict>: dictionary of configuration options for resample.
    Returns:
        config <dict>: dictionary of configuration options with parsed values.

    """
    for key in config:
        config[key] = _parse_list(config[key], dtype=float)

    return config


def parse_norml(config):
    """ Parse arguments for Normalization.

    Args:
        config <dict>: dictionary of configuration options for norml.
    Returns:
        config <dict>: dictionary of configuration options with parsed values.

    """

    for key in config:
        if key == 'udefined':
            config[key] = _parse_list(config[key], dtype=bool)
        else:
            config[key] = _parse_list(config[key])

    return config


def parse_clip(config):
    """ Parse arguments for Clipping.

    Args:
        config <dict>: dictionary of configuration options for clip.
    Returns:
        config <dict>: dictionary of configuration options with parsed values.

    """

    # Note: parameter 'substitute' can contain a None and float values.
    for key in config:
        config[key] = _parse_list(config[key], dtype=float)

    return config


def parse_trim(config):
    """ Parse arguments for trim operation.

    Args:
        config <dict>: dictionary of configuration options for trim.
    Returns:
        config <dict>: dictionary of configuration options with parsed values.

    """
    config['bins'] = _parse_list_of_lists(config['bins'], delimiter_elements='-', delimiter_lists=',')
    return config



def parse_rnv(config):
    """ Parse arguments for standard normal variate scatter correction.

    Args:
        config <dict>: dictionary of configuration options for snv.
    Returns:
        config <dict>: dictionary of configuration options with parsed values.

    """

    if 'iqr' in config.keys():
        config['iqr'] = _parse_list_of_lists(config['iqr'])
    return config


def parse_lsnv(config):
    """ Parse arguments for standard normal variate scatter correction.

    Args:
        config <dict>: dictionary of configuration options for snv.
    Returns:
        config <dict>: dictionary of configuration options with parsed values.

    """
    if 'num_windows' in config.keys():
        config['num_windows'] = _parse_list(config['num_windows'], dtype=int)
    return config


def parse_detrend(config):
    config['bp'] = _parse_list_of_lists(config['bp'], dtype=int)

    return config


def parse_derivate(config):
    """ Parse arguments for N-th order derivation.

    Args:
        config <dict>: dictionary of configuration options for derivate.
    Returns:
        config <dict>: dictionary of configuration options with parsed values.

    """

    if 'order' in config.keys():
        config['order'] = _parse_list(config['order'], dtype=int)

    if 'delta' in config.keys():
        config['delta'] = _parse_list(config['delta'], dtype=int)

    return config


def parse_smooth(config):
    """ Parse arguments for moving average filtering.

    Args:
        config <dict>: dictionary of configuration options for smooth.
    Returns:
        config <dict>: dictionary of configuration options with parsed values.

    """
    config['filter_win'] = _parse_list(config['filter_win'], dtype=int)

    if 'mode' in config.keys():
        config['mode'] = _parse_list(config['mode'], dtype=str)

    if 'window_type' in config.keys():
        config['window_type'] = _parse_list(config['window_type'], dtype=str)

    return config


# ------- PARSING UTILITY FUNCTIONS --------
def _parse_list_of_lists(string, delimiter_elements=',', delimiter_lists=':', delimiter_pipelines=';', dtype=float):
    """ Parses a string that contains single or multiple lists.

    Args:
        delimiter_elements <str>: delimiter between inner elements of a list.
        delimiter_lists <str>: delimiter between lists.
        delimiter_pipelines <str>: delimiter between different pipelines.

    Returns:
        new_list <list> parsed list of configuration parameters.
    """
    new_list = []
    for sub_list in string.strip().replace(' ', '').split(delimiter_pipelines):
        if delimiter_lists in sub_list:
            new_list.append([_parse_list(item, dtype=dtype, delimiter=delimiter_elements) for item in sub_list.split(delimiter_lists)])
        else:
            new_list.append(_parse_list(sub_list, dtype=dtype, delimiter=delimiter_elements))
    return new_list


def _parse_list(string, dtype=int, delimiter=','):
    """ Converts a string to a list (of specified data type).

    Args:
        string <str>: string containing arguments.
        dtype <type>: data type of the arguments.
        delimiter <str>: delimiter between arguments.

    Returns:
        items <list>: list of arguments.

    """

    items = string.lower().strip().replace(' ', '').split(delimiter)

    if 'none' in items:
        items.pop(items.index('none'))
        contains_none = True
    else:
        contains_none = False


    if dtype == bool:
        items = [item == 'true' for item in items]
    else:
        items = [dtype(item) for item in items]

    if contains_none:
        items.append(None)

    return items
# ------------------------------------------


def parse_section(config, config_type):
    """ Parse different sections of the configuration file.

    Args:
        config <dict>: dictionary containing un-parsed configuration information for all sections.
        config_type <str>: name of the section (e.g., SAVGOL, SNV, etc) to parse.
    Returns:
        config <dict>: dictionary with parsed section.
    """


    if 'also_skip' in config:
        also_skip =  config['also_skip'].lower() == 'true' or config['also_skip'].lower() == '1'
        config.pop('also_skip')
    else:
        also_skip = False

    if config_type == 'SAVGOL':
        config = parse_savgol(config)
    elif config_type == 'SNV':
        config = {}
    elif config_type == 'RNV':
        config = parse_rnv(config)
    elif config_type == 'LSNV':
        config = parse_lsnv(config)
    elif config_type == 'TRIM':
        config = parse_trim(config)
    elif config_type == 'DETREND':
        config = parse_detrend(config)
    elif config_type == 'MSC':
        config = {}
    elif config_type == 'NORML':
        config = parse_norml(config)
    elif config_type == 'CLIP':
        config = parse_clip(config)
    elif config_type == 'SMOOTH':
        config = parse_smooth(config)
    elif config_type == 'RESAMPLE':
        config = parse_resample(config)
    else:
        raise TypeError('Preprocessing option "{}" not recognized!'.format(config_type))

    if also_skip:
        config['also_skip'] = also_skip

    return config


def construct_pipelines(config):
    """ Uses the (parsed) configuration dict to generate a preprocessing pipeline.

    Args:
        config <dict>: parsed contents of a configuration file.

    Returns:
        pipelines <list>: a list of nippy preprocessing pipelines.
    """


    def _get_argument_combinations(arguments):
        """ Utility to function to obtain all permutations of preprocessing arguments. """
        arg_names = sorted(arguments)
        combinations = itertools.product(*(arguments[arg] for arg in arg_names))
        combinations = [dict(zip(arg_names, arg_values)) for arg_values in combinations]
        return combinations

    options = {}
    for key in config.keys():
        # 1. Check if we got also_skip
        if 'also_skip' in config[key] and config[key]['also_skip']:
            config[key].pop('also_skip')
            options[key] = _get_argument_combinations(config[key])
            options[key].append(None)
        else:
            options[key] = _get_argument_combinations(config[key])

    return _get_argument_combinations(options)


def remove_incompatible_operations(pipelines):
    """ Removes preprocessing combinations that should not occur in a single pipeline.

    Args:
        pipelines <list>: list of nippy pipelines.

    Returns:
        pipelines <list>: list of nippy pipelines without pipelines containing incompatible operations.
    """

    def find_duplicates(pipelines):
        for idx in range(len(pipelines)):
            for idx_ in range(idx + 1, len(pipelines)):
                if pipelines[idx] == pipelines[idx_]:
                    return idx
        return -1

    def check_pair(arg1, arg2, pipe):
        if arg1 in pipe.keys() and arg2 in pipe.keys(): # Might be a problem
            if pipe[arg1] == None or pipe[arg2] == None: # Not a problem
                is_problem = False
            else:  # Is a problem
                is_problem = True
        else:  # Not a problem
            is_problem = False
        return is_problem

    # Remove illegal combinations
    # FIXME: Come up with a smarter system for checking illegal pairs
    bad_pairs = [('MSC', 'SNV'), ('MSC', 'RNV'), ('SNV', 'RNV'), ('SMOOTH', 'SAVGOL'), ('LSNV', 'SNV'), ('MSC', 'LSNV'), ('RNV', 'LSNV')]
    bad_idx = []
    new_pipes = []
    for bad_pair in bad_pairs:
        a, b = bad_pair
        for idx, pipeline in enumerate(pipelines):
            if check_pair(a, b, pipeline):
                pipeline_b = pipeline.copy()
                pipeline[a] = None
                pipeline_b[b] = None
                new_pipes.append(pipeline_b)

    pipelines.extend(new_pipes)
    # remove duplicates
    while find_duplicates(pipelines) != -1:
        pipelines.pop(find_duplicates(pipelines))

    return pipelines


def read_configuration(file_path):
    """ Read and parse configuration for a preprocessing pipeline from an INI-file.

    Args:
        file_path <str>: file path to the configuration file.
    Returns:
        config <dict>: dictionary of all preprocessing sections and arguments.
    """
    parser = configparser.ConfigParser()
    parser.read(file_path)

    # Parse predefined configuration sections
    config = {}
    for part in ['SAVGOL', 'TRIM', 'SNV', 'RNV', 'LSNV', 'DETREND', 'MSC', 'NORML', 'CLIP', 'SMOOTH', 'RESAMPLE']:
        if part in parser:
            config[part] = parse_section(dict(parser[part]), part)

    pipelines = construct_pipelines(config)
    pipelines = remove_incompatible_operations(pipelines)

    return pipelines
