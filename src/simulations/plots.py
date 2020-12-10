class plots:
    def __init__(self):
        return

class plot_rocs(rocs):
    def __init__(self):
        return

    def plot_multiple_rocs(self):
        ## Needed for putting all curves in same space
        base_fpr = np.linspace(0, 1, 101)
        tprs = []
        for roc in rocs:
            fpr, tpr = roc[0], roc[1]
            tpr = interp(base_fpr, fpr, tpr)
            tpr[0] = 0.0
            tprs.append(tpr)

        tprs = np.array(tprs)
        print(tprs)
        print(len(tprs))
        mean_tprs = tprs.mean(axis=0)
        std = tprs.std(axis=0)