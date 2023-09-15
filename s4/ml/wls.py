"""Simple linear models (OLS) and weighted least squares (WLS)"""
import logging
from multiprocessing import cpu_count
from multiprocessing.pool import Pool
from typing import Optional, List

import matplotlib.pyplot as plt
import numpy
import pandas
import statsmodels.formula.api as sm
from IPython.display import display
from pandas import DataFrame
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
from tqdm import tqdm

from s4.ml.utils import weights2colors

__all__ = [
    'linear_model_analysis',
    'forward_selection_by_bic',
    'feature_importance',
]


def do_linear_model(data_frame, features, y_name, weights, do_wls):
    """Construct a linear model using statsmodels"""
    formula = f'{y_name} ~ {"+".join(features)}'
    if do_wls:
        return sm.wls(formula=formula, data=data_frame, weights=weights).fit()

    return sm.ols(formula=formula, data=data_frame).fit()


def forward_selection_by_bic(
        data_frame: DataFrame, features_sorted: List[str],
        y_name: str, weights: Optional[numpy.ndarray] = None, do_wls: bool = False):
    """
    Perform feature forward selection by optimizing BIC (Bayesian Information Criteria).
    The features are taken according to the order defined by `features_sorted`.

    :param data_frame: The training data frame.
    :param features_sorted: Features to use. Features are taken in order.
    :param y_name: Prediction target name.
    :param weights: If `do_wls=True`, use this as sample weights.
    :param do_wls: Whether to use Weighted Least Squares.
    """
    current_features = []
    current_bic = 9e99

    assert len(features_sorted) != 0, "Input feature list is empty!"

    for feature in features_sorted:
        current_features.append(feature)
        result = do_linear_model(
            data_frame=data_frame,
            features=current_features,
            y_name=y_name,
            do_wls=do_wls, weights=weights
        )
        if result.bic < current_bic:
            current_bic = result.bic
        else:
            current_features.pop(-1)

    return current_features, result


def _setup_inc_r2(data, weights):
    _setup_inc_r2.data = data
    _setup_inc_r2.weights = weights


def _compute_adj_r2(data_x, data_y, weights):
    n_samples = len(data_x)
    n_features = len(data_x[0])

    if data_x.shape[1] == 0:
        # just intercept
        y_pred = numpy.array([numpy.average(data_y, weights=weights)] * data_y.size)
        raw_r2 = r2_score(data_y, y_pred, sample_weight=weights)
    else:
        base = LinearRegression(copy_X=True)
        base.fit(data_x, data_y, sample_weight=weights)
        raw_r2 = base.score(data_x, data_y, sample_weight=weights)
    adj_r2 = 1 - (1 - raw_r2) * ((n_samples - 1) / (n_samples - n_features - 1))
    return adj_r2


def _incremental_r2(args):
    features, y_name, inc_feature = args
    data_x = _setup_inc_r2.data[features].values
    data_x_add = _setup_inc_r2.data[features + [inc_feature]].values
    data_y = _setup_inc_r2.data[y_name].values
    weights = _setup_inc_r2.weights

    base_adj_r2 = _compute_adj_r2(data_x, data_y, weights)
    inc_adj_r2 = _compute_adj_r2(data_x_add, data_y, weights)
    return inc_adj_r2 - base_adj_r2


def estimate_average_partial_importance(  # pylint: disable=too-many-arguments
        data_frame, features, y_name, weights, features_inspect,
        submodel_search=100, processes=None):
    """Estimate partial importance using sub-model sampling method."""
    partial = {}
    with Pool(processes=processes, initializer=_setup_inc_r2,
              initargs=(data_frame, weights)) as pool:
        futures = []
        for feature_inspect in tqdm(features_inspect, desc='Submitting jobs'):
            complement_model = [
                feature for feature in features if feature != feature_inspect]
            submodels = [
                [feature for i, feature in enumerate(complement_model) if y[i] > 0]
                for y in numpy.random.randint(2, size=(submodel_search, len(complement_model)))
            ]

            future = pool.map_async(_incremental_r2,
                                    [(submodel, y_name, feature_inspect) for submodel in submodels],
                                    chunksize=len(submodels) // cpu_count())
            futures.append((feature_inspect, future))

        for feature_inspect, future in tqdm(futures, desc='Waiting for results'):
            partial[feature_inspect] = numpy.mean(future.get())
    return partial


def linear_feature_selection(  # pylint: disable=too-many-arguments
        data_frame: DataFrame, candidate_features: List[str], y_name: str, *,
        do_wls: bool = True, weights: Optional[numpy.ndarray] = None,
        submodel_search: int = 100, alpha: float = 0.05, processes: int = None):
    """
    Perform feature selection of linear models.

    1. All features with low F-statistic are dropped.
    2. Starting with full model, features with low t-statistic are dropped.

    :param data_frame: DataFrame containing training data.
    :param candidate_features: Candidate features.
    :param y_name: Name of y column.
    :param alpha: Significance level.
    :param do_wls: Perform WLS instead of OLS.
    :param weights: Weights of each sample.
    :param submodel_search: Size of submodel delta R estimation.
    :param processes: Number of processes to use.
    :return:
    """

    features = []
    for fea in candidate_features:
        if data_frame[fea].unique().size > 1:
            features.append(fea)

    print('Reducing features...', end='', flush=True)
    while len(features) > 0:
        print(f'\rReducing features {len(features)}/{len(candidate_features)}... ', end='',
              flush=True)
        result = do_linear_model(data_frame, features, y_name, weights, do_wls)
        coefs = [(val, name) for name, val in result.pvalues.items() if
                 name != 'Intercept' and val >= alpha]
        if not coefs:
            # all coefficients have sufficiently small p-value
            print('Done')
            return features, result

        features_inspect = [x[1].split('[')[0] for x in coefs]
        # remove the one that minimizes delta r
        partial_avg_importance = estimate_average_partial_importance(
            data_frame=data_frame, features=features, y_name=y_name, weights=weights,
            features_inspect=features_inspect, submodel_search=submodel_search,
            processes=processes
        )

        logging.debug('Current delta R2 stats: %r', partial_avg_importance)
        least_informative = min(partial_avg_importance, key=lambda x: partial_avg_importance[x])
        logging.debug('Removing %r because it has delta adj R2 of %r',
                      least_informative, partial_avg_importance[least_informative])
        features.remove(least_informative)

    print('')
    raise ValueError('After feature selection, no feature was significant.')


def feature_importance(  # pylint: disable=too-many-arguments
        data_frame: DataFrame, features: List[str], y_name: str, *,
        do_wls: bool = True, weights=None, submodel_search=100, processes=None):
    """
    Estimate average partial importance.

    The feature importance is calculated by sampling N=submodel_search submodels,
    since the total size of the submodel space may be huge to sample and calculate.

    Args:
        data_frame: Training data frame.
        features: List of features to use.
        y_name: Prediction target name.
        do_wls: Whether or not to perform WLS.
        weights: Sample weights.
        submodel_search: Number of sub-models when computing importance.
        processes: Number of processes.
    """
    # pylint: disable=too-many-locals
    full_model = do_linear_model(data_frame, features, y_name, weights, do_wls)

    interactional = []
    # interactional
    for fea in features:
        other_features = [x for x in features if x != fea]
        sub_model = do_linear_model(data_frame, other_features, y_name, weights, do_wls)
        interactional.append(full_model.rsquared_adj - sub_model.rsquared_adj)

    # individual
    individual = []
    for fea in features:
        sub_model = do_linear_model(data_frame, [fea], y_name, weights, do_wls)
        individual.append(sub_model.rsquared_adj)

    # estimated average partial
    partial = estimate_average_partial_importance(
        data_frame=data_frame, features=features, y_name=y_name, weights=weights,
        features_inspect=features, submodel_search=submodel_search, processes=processes
    )
    partial = [partial[x] for x in features]

    total = numpy.array([(a + b + c) for a, b, c in zip(interactional, individual, partial)])
    pct = total / total.sum()

    result = pandas.DataFrame({
        'Interactional Dominance': interactional,
        'Individual Dominance': individual,
        'Estimated Average Partial Dominance': partial,
        'Total Dominance': total,
        'Percentage Relative Importance': pct,
    }, index=features)

    return result.sort_values(by='Total Dominance', ascending=False)


def linear_model_analysis(
        data_frame: DataFrame, features: List[str], y_name: str, y_desc: str, *,
        do_wls: bool = True, weights: Optional[numpy.ndarray] = None,
        do_feature_selection: bool = True, do_plot: bool = True, processes=None):
    """
    Perform feature importance analysis and feature selection on a dataset using OLS/WLS.

    Args:
        data_frame: Training data frame.
        features: List of features to use.
        y_name: Prediction target name.
        y_desc: Description of prediction target, will used in plots.
        do_wls: Whether or not to perform WLS.
        weights: Sample weights.
        do_feature_selection: Whether or not to perform feature selection.
        do_plot: Whether or not to plot data.
        processes: Number of processes to use.
    """

    # pylint: disable=too-many-locals
    if do_wls and weights is None:
        raise ValueError('Must have weights when performing WLS')

    if do_feature_selection:
        subset_features, result = linear_feature_selection(
            data_frame=data_frame, candidate_features=features, y_name=y_name, alpha=0.05,
            do_wls=do_wls, weights=weights, processes=processes)
    else:
        subset_features = features
        result = do_linear_model(data_frame=data_frame, features=features, y_name=y_name,
                                 weights=weights,
                                 do_wls=do_wls)
    print(result.summary())

    importance = feature_importance(
        data_frame, features=subset_features, y_name=y_name,
        do_wls=do_wls, weights=weights, processes=processes)
    display(importance.style.background_gradient(cmap='coolwarm', axis=0))

    features_corr = data_frame[subset_features].corr()

    if do_plot:
        _, (ax1, ax2) = plt.subplots(ncols=2, figsize=(10, 4))
        ax1.scatter(data_frame[y_name], result.predict(data_frame), color=weights2colors(weights))
        ax1.set_xlabel('Observed ' + y_desc)
        ax1.set_ylabel('Predicted ' + y_desc)

        ax2.hist(result.resid, bins=50)
        ax2.set_xlabel('Residual')
        plt.show()

        display(features_corr.style.background_gradient(cmap='coolwarm'))

    return {
        'Selected features': subset_features,
        'Fit result': result,
        'Importance': importance,
        'Features correlation': features_corr
    }
