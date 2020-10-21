from .utils.config import check_required


def jackhmmer_search(**kwargs):
    """
    Protocol:
    Iterative jackhmmer search against a sequence database.
    Parameters
    ----------
    Mandatory kwargs arguments:
        See list below in code where calling check_required
    .. todo::
        explain meaning of parameters in detail.
    Returns
    -------
    outcfg : dict
        Output configuration of the protocol, including
        the following fields:
        * sequence_id (passed through from input)
        * first_index (passed through from input)
        * target_sequence_file
        * sequence_file
        * raw_alignment_file
        * hittable_file
        * focus_mode
        * focus_sequence
        * segments
    """
    check_required(
        kwargs,
        [
            "prefix", "sequence_id", "sequence_file",
            "sequence_download_url", "region", "first_index",
            "use_bitscores", "domain_threshold", "sequence_threshold",
            "database", "iterations", "cpu", "nobias", "reuse_alignment",
            "checkpoints_hmm", "checkpoints_ali", "jackhmmer",
            "extract_annotation"
        ]
    )
    prefix = kwargs["prefix"]


    # make sure output directory exists
    create_prefix_folders(prefix)

    # store search sequence file here
    target_sequence_file = prefix + ".fa"
    full_sequence_file = prefix + "_full.fa"
