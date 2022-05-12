files_map = {}


class File(object):
    def __init__(self, analysis, name, f_in, f_out, params=None):
        self.analysis = analysis
        self.name = name
        self.f_in = f_in
        self.f_out = f_out
        self.params = params

        self.load()
        return

    def load(self, *args, **kwargs):
        return

class BtwnClust(File):
    def __init__(self, analysis, name, f_in, f_out, params=None, **kwargs):
        #super.__init__(kwargs)
        super().__init__(analysis, name, f_in, f_out, params)
        return

    def load(self, **kwargs):
        self.f_in = ""
        return


# Analysis for single runs:
def an_dendro_barcodes_clones():
    out = {}
    out["barcode_dendro"] = "barode.png"
    return


def an_umaps():
    return


def an_btwncluster_gsea():
    return


def an_mt_as_clones():
    return