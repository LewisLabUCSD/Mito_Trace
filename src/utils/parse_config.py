"""
Configuration handling
.. todo::
    switch ruamel.yaml to round trip loading
    to preserver order and comments?
Authors:
  Thomas A. Hopf
"""

import ruamel.yaml as yaml


class MissingParameterError(Exception):
    """Exception for missing parameters"""


class InvalidParameterError(Exception):
    """Exception for invalid parameter settings"""


def parse_config(config_str, preserve_order=False):
    """Parse a configuration string :param config_str: Configuration to be
    parsed :type config_str: str :param preserve_order: Preserve formatting of
    input configuration

        string

    Args:
        config_str:
        preserve_order:

    Returns:
        dict: Configuration dictionary
    """
    try:
        if preserve_order:
            return yaml.load(config_str, Loader=yaml.RoundTripLoader)
        else:
            return yaml.safe_load(config_str)
    except yaml.parser.ParserError as e:
        raise InvalidParameterError(
            "Could not parse input configuration. "
            "Formatting mistake in config file? "
            "See ParserError above for details."
        ) from e


def read_config_file(filename, preserve_order=False):
    """Read and parse a configuration file. :param filename: Path of
    configuration file :type filename: str

    Args:
        filename:
        preserve_order:

    Returns:
        dict: Configuration dictionary
    """
    with open(filename) as f:
        return parse_config(f, preserve_order)


def write_config_file(out_filename, config):
    """Save configuration data structure in YAML file. :param out_filename:
    Filename of output file :type out_filename: str :param config: Config data
    that will be written to file :type config: dict

    Args:
        out_filename:
        config:
    """
    if isinstance(config, yaml.comments.CommentedBase):
        dumper = yaml.RoundTripDumper
    else:
        dumper = yaml.Dumper

    with open(out_filename, "w") as f:
        f.write(
            yaml.dump(config, Dumper=dumper, default_flow_style=False)
        )


def check_required(params, keys):
    """Verify if required set of parameters is present in configuration :param
    params: Dictionary with parameters :type params: dict :param keys: Set of
    parameters that has to be present in params :type keys: list-like

    Args:
        params:
        keys:

    Raises:
        MissingParameterError
    """
    missing = [k for k in keys if k not in params]

    if len(missing) > 0:
        raise MissingParameterError(
            "Missing required parameters: {} \nGiven: {}".format(
                ", ".join(missing), params
            )
        )


def iterate_files(outcfg, subset=None):
    """Generator function to iterate a list of file items in an outconfig :param
    outcfg: Configuration to extract file items for iteration from :type outcfg:
    dict(str) :param subset: List of keys in outcfg to restrict iteration to
    :type subset: list(str)

    Args:
        outcfg:
        subset:

    Returns:
        tuple(str, str, int): Generator over tuples (file path, entry key,
        index). index will be None if this is a single file entry (i.e. ending
        with _file rather than _files).
    """
    for k, v in outcfg.items():
        # skip items if there is a subset filter and it matches
        if subset is not None and k not in subset:
            continue

        # also skip in case file has a null value
        if v is None:
            continue

        # only look at file entries, so skip everything else
        # if not (k.endswith("_file") or k.endswith("_files")):
        #    continue
        if k.endswith("_file"):
            yield (v, k, None)
        elif k.endswith("_files"):
            for i, f in enumerate(v):
                yield (f, k, i)
        else:
            # skip any other entries
            pass


# def validate_schema_files(params_f, schema_f, to_normalize=True):
#     """ Validates a parameter file and forces proper types.
#
#     This will compare params file to a schema in json format.
#     Raises an exception if the schema is not valid.
#     See cerberus for additional information on schema structure.
#     :param params_f: Parameter file in yaml/json format
#     :param schema_f: Schema file in json formate based on cerberus schema.
#     :param to_normalize: If True, will coerve to proper types when
#     applicable.
#     :raise ValidationError: If documaent is not valid.
#     :return: params
#     """
#
#     params = read_config_file(params_f)
#     schema = read_config_file(schema_f)
#     print(schema["PROJECT"]["rename_handler"])
#
#     # Execute the rename_handler, which loads a rename_handler function,
#     # and sets it.
#     exec(schema["PROJECT"]["rename_handler"])
#     exec(f'schema["PROJECT"]["rename_handler"] = {schema["PROJECT"]["rename_handler"]}')
#
#     worked, v = validate_schema(params, schema)
#     print('v',v)
#     if to_normalize:
#         params = v.document
#     if not worked:
#         print("Errors", v.errors)
#         raise(ValidationError)
#     else:
#         return params