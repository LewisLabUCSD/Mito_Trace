import click

# def next(hey):
#     exec("sample=" + hey)
#     print(sample)
#     return
# @click.command()
# @click.option('--hey', default="", type=click.STRING)
# def main(hey):
#     print('hey arg')
#     print(hey)
#     print(type(hey))
#     print("sample="+hey)
#     exec("sample="+hey)
#     ldic = locals()
#     exec("sample="+hey, globals(), ldic)
#     sample = ldic["sample"]
#     print(sample)  # it works! returns 2
#
#     next(hey)
#     #print(sample)
#     return
#
#

def main():
    from src.utils.data_io import load_mtx_df
    df = load_mtx_df("Analysis/multiplex/data/CHIP_april08_2021/MTblacklist/chrM/pseudo/minC200_minAF0.01/numC100000_ispropFalse/cellSNP.tag.AD.mtx")
    print('df')
    print(df)
    return
if __name__ == '__main__':
    main()
