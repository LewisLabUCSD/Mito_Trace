"""
Specifically, the Kullback–Leibler divergence from Q to P, denoted DKL(P‖Q), is
a measure of the information gained when one revises one's beliefs from the
prior probability distribution Q to the posterior probability distribution P. In
other words, it is the amount of information lost when Q is used to approximate
P.
"""
import numpy as np
from scipy.stats import entropy

def kl(p, q):
    """Kullback-Leibler divergence D(P || Q) for discrete distributions
    Parameters
    ----------
    p, q : array-like, dtype=float, shape=n
    Discrete probability distributions.
    """
    p = np.asarray(p, dtype=np.float)
    q = np.asarray(q, dtype=np.float)

    return np.sum(np.where(p != 0, p * np.log(p / q), 0))


p = [0.1, 0.9]
q = [0.1, 0.9]
assert entropy(p, q) == kl(p, q)