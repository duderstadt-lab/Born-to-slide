import numpy as np
from sklearn.utils import resample


def bootstrap(data, n_boot=10000, sample_size=1, estimator=np.mean):
    """
    :param data: array with data
    :param n_boot: number for bootstrapping iterations (default 10000)
    :param sample_size: sample coverage ]0;1] (default 1)
    :param estimator: default np.mean
    :return: list of bootstrap samples
    """
    return estimator(
        [resample(data, replace=True, n_samples=int(sample_size * len(data))) for _ in range(n_boot)],
        axis=1)


def calc_ci(data, ci=95):
    """
    Calculates values for confidence interval
    :param data: arrayed data
    :param ci: confidence interval (default 95)
    :return: lower_bound, upper_bound
    """
    return np.percentile(data, 50 - ci / 2), np.percentile(data, 50 + ci / 2)


def significance(p):
    """
    Returns significance symbol based on set alpha values
    :param p: probability of statistical test
    :return: string expression for significance
    """
    if p < 0.001:
        expression = "***"
    elif p < 0.01:
        expression = "**"
    elif p < 0.05:
        expression = "*"
    else:
        expression = "ns"
    return expression
