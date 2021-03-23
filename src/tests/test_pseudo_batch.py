import click

def next(hey):
    exec("sample=" + hey)
    print(sample)
    return
@click.command()
@click.option('--hey', default="", type=click.STRING)
def main(hey):
    print('hey arg')
    print(hey)
    print(type(hey))
    print("sample="+hey)
    exec("sample="+hey)
    ldic = locals()
    exec("sample="+hey, globals(), ldic)
    sample = ldic["sample"]
    print(sample)  # it works! returns 2

    next(hey)
    #print(sample)
    return


if __name__ == '__main__':
    main()
