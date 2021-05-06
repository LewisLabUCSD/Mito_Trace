import pandas as pd


def preprocess_variants(variants, style=">"):
    print('variants')
    print(variants)
    if style == ">":
        def split(x):
            s = x.split(">")
            return [s[0][:-1], s[0][-1], x[-1]]

        curr = pd.DataFrame(list(map(split, variants)),
                            columns=["position", "ref", "alt"],
                            index=variants)
    else:
        print(f"style {style} not implemented")
        return
    return curr


def type_of_variants(variants, style=">"):
    variants = preprocess_variants(variants, style=style)
    # Get types of mutations
    def var_type(x):
        nts = set(x[["ref", "alt"]])
        if nts == {"A","G"} or nts == {"T", "C"}:
            return "Transition"
        return "Transversion"
    variants["variant type"] = variants.apply(var_type, axis=1)

    variants["variant change"] = variants["ref"]+">"+variants["alt"]
    return variants
