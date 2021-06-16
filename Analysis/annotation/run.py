from src.mlflow.clones import MTExp
import src.mlflow.mlflow_utils as mu
from mlflow.utils import mlflow_tags
import os
import mlflow
import click

@click.command()
@click.argument("entry", type=click.Choice(["main", "preprocess", "results"]))
@click.argument("configfile", type=click.Path(exists=True))
@click.option("--preprocmode", default="filters", type=click.Choice(["filters", "vireo"]))
@click.option("--datauri", default="", type=click.Path())
def main(entry,configfile, preprocmode, datauri):
    protocol = MTExp(config_f=configfile)
    protocol.initialize()
    if entry == "main":
        with mlflow.start_run() as active_run:
            git_commit = active_run.data.tags.get(
                mlflow_tags.MLFLOW_GIT_COMMIT)

            load_raw_data_run = mu._get_or_run("preprocess", {'parameters_f': configfile, 'preprocmode':preprocmode},
                                            git_commit)
            print('load_raw_data_run', load_raw_data_run)
            mtexp_uri = os.path.join(
                load_raw_data_run.info.artifact_uri, "data/intermediate")
            print('mtexp_uri', mtexp_uri)
            _ = mu._get_or_run("results",
                {"data_uri": mtexp_uri,
                 'parameters_f': configfile}, git_commit)

        protocol.is_completed=True

    elif entry == "preprocess":
        print('preprocessing')
        protocol.preprocess(preprocmode)
    elif entry == 'results':
        protocol.results(datauri)
    else:
        raise("incorrect entry point")
    return


if __name__ == "__main__":
    main()
    # main(['mtexp_a.yaml',  'main',
    #       '--preproc_mode', 'filters',
    #       '--datauri', '.'])
    # main(['mtexp_a.yaml',  'main',
    #       '--preproc_mode', 'filters',
    #       '--datauri', '.'])
