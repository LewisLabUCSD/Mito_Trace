class Paramspace:
    """A wrapper for pandas dataframes that provides helpers for using them as a parameter
    space in Snakemake.

    This is heavily inspired by @soumitrakp work on JUDI (https://github.com/ncbi/JUDI).

    By default, a directory structure with on folder level per parameter is created
    (e.g. column1~{column1}/column2~{column2}/***).

    The exact behavior can be tweeked with two parameters:

      - ``filename_params`` takes a list of column names of the passed dataframe.
        These names are used to build the filename (separated by '_') in the order
        in which they are passed.
        All remaining parameters will be used to generate a directoty structure.
        Example for a data frame with four columns named column1 to column4:

        | ``Paramspace(df, filename_params=["column3", "column2"])`` ->
        | column1~{value1}/column4~{value4}/column3~{value3}_column2~{value2}

        If ``filename_params="*"``, all columns of the dataframe are encoded into
        the filename instead of parent directories.

      - ``param_sep`` takes a string which is used to join the column name and
        column value in the genrated paths (Default: '~'). Example:

        | ``Paramspace(df, param_sep=":")`` ->
        | column1:{value1}/column2:{value2}/column3:{value3}/column4:{value4}
    """

    def __init__(self, dataframe, filename_params=None, param_sep="~"):
        self.dataframe = dataframe
        self.param_sep = param_sep
        if filename_params is None or not filename_params:
            # create a pattern of the form {}/{}/{} with one entry for each
            # column in the dataframe
            self.pattern = "/".join([r"{}"] * len(self.dataframe.columns))
            self.ordered_columns = self.dataframe.columns
        else:
            if isinstance(filename_params, str) and filename_params == "*":
                filename_params = dataframe.columns

            if any((param not in dataframe.columns for param in filename_params)):
                raise KeyError(
                    "One or more entries of filename_params are not valid coulumn names for the param file."
                )
            elif len(set(filename_params)) != len(filename_params):
                raise ValueError("filename_params must be unique")
            # create a pattern of the form {}/{}_{} with one entry for each
            # column in the dataframe. The number of underscore-separated
            # fields is equal to the number filename_params
            self.pattern = "/".join(
                [r"{}"] * (len(self.dataframe.columns) - len(filename_params) + 1)
            )
            self.pattern = "_".join(
                [self.pattern] + [r"{}"] * (len(filename_params) - 1)
            )
            self.ordered_columns = [
                param
                for param in self.dataframe.columns
                if param not in filename_params
            ]
            self.ordered_columns.extend(list(filename_params))
        self.dataframe = self.dataframe[self.ordered_columns]

    @property
    def wildcard_pattern(self):
        """Wildcard pattern over all columns of the underlying dataframe of the form
        column1~{column1}/column2~{column2}/*** or of the provided custom pattern.
        """
        return self.pattern.format(
            *map(self.param_sep.join(("{0}", "{{{0}}}")).format, self.ordered_columns)
        )

    @property
    def instance_patterns(self):
        """Iterator over all instances of the parameter space (dataframe rows),
        formatted as file patterns of the form column1~{value1}/column2~{value2}/...
        or of the provided custom pattern.
        """
        return (
            self.pattern.format(
                *(
                    self.param_sep.join(("{}", "{}")).format(name, value)
                    for name, value in row.items()
                )
            )
            for index, row in self.dataframe.iterrows()
        )

    def instance(self, wildcards):
        """Obtain instance (dataframe row) with the given wildcard values."""
        import pandas as pd

        def convert_value_dtype(name, value):
            if self.dataframe.dtypes[name] == bool and value == "False":
                # handle problematic case when boolean False is returned as
                # boolean True because the string "False" is misinterpreted
                return False
            else:
                return pd.Series([value]).astype(self.dataframe.dtypes[name])[0]

        return {
            name: convert_value_dtype(name, value)
            for name, value in wildcards.items()
            if name in self.ordered_columns
        }


    def __getattr__(self, name):
        import pandas as pd

        ret = getattr(self.dataframe, name)
        if isinstance(ret, pd.DataFrame):
            return Paramspace(ret)
        return ret

    def __getitem__(self, key):
        import pandas as pd

        ret = self.dataframe[key]
        if isinstance(ret, pd.DataFrame):
            return Paramspace(ret)
        return ret
