import click

@click.command()
@click.argument('params', type=click.STRING, nargs=-1)
#@click.argument('topn')
# @click.argument('minreads')
# @click.argument('hetthresh')
# @click.argument('minhetcells')
# @click.argument('hetcountthresh')
# @click.argument('bqthresh')
# @click.argument('lowcovthresh')
# @click.argument('ncellsthresh')
def main(params):
    print(params)
    #print(exec(params))
    return

if __name__ == "__main__":
    main()