"""
Functions for creating nonlinear models for synthesis prediction.
"""
from itertools import product
from multiprocessing import cpu_count
from multiprocessing.pool import Pool
from typing import List, Union, Optional

import matplotlib.pyplot as plt
import numpy
from pandas import DataFrame
from sklearn import clone, preprocessing
from sklearn.base import RegressorMixin
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import r2_score
from sklearn.model_selection import LeaveOneOut, KFold
from tqdm import tqdm

from s4.ml.utils import weights2colors


def _setup_cv(data_x, data_y, weights, estimator, do_normalization):
    """Setup cross-validation calculation in a subprocess."""
    _setup_cv.data_x = data_x
    _setup_cv.data_y = data_y
    _setup_cv.weights = weights
    _setup_cv.estimator = estimator
    _setup_cv.do_normalization = do_normalization

    if do_normalization:
        _setup_cv.data_x = preprocessing.normalize(data_x, norm='l2')


def _run_cv(args):
    """Run cross validation for a pair of (train_indices, test_indices)"""
    train_ind, test_ind = args

    train_x, train_y = _setup_cv.data_x[train_ind], _setup_cv.data_y[train_ind]
    test_x, test_y = _setup_cv.data_x[test_ind], _setup_cv.data_y[test_ind]
    if _setup_cv.weights is not None:
        train_weights = _setup_cv.weights[train_ind]
    else:
        train_weights = None

    regressor = clone(_setup_cv.estimator)
    regressor.fit(train_x, train_y, sample_weight=train_weights)

    return {
        # 'OOB score': r2_score(train_y, regressor.oob_prediction_, sample_weight=train_weights),
        'Train residual': train_y - regressor.predict(train_x),
        'Test residual': test_y - regressor.predict(test_x),
        'Train indices': train_ind,
        'Test indices': test_ind,
    }


def _fast_cv(estimator, data_x, data_y, weights, cv_method):
    """Perform all cross validation calculations in one single process"""
    if cv_method == 'loo':
        cv_split = LeaveOneOut()
    elif isinstance(cv_method, int):
        cv_split = KFold(cv_method)
    else:
        raise ValueError(f'Invalid CV {cv_method}')

    test_val = numpy.empty_like(data_y)
    for train_ind, test_ind in cv_split.split(data_x):
        model = clone(estimator)
        model.fit(data_x[train_ind], data_y[train_ind], sample_weight=weights[train_ind])
        test_val[test_ind] = model.predict(data_x[test_ind])

    return r2_score(data_y, test_val, sample_weight=weights)


def estimate_average_partial_importance(  # pylint: disable=too-many-arguments
        estimator: RegressorMixin, data_frame: DataFrame, features: List[str], y_name: str, *,
        cv_method: Union[int, str] = 'loo', weights: Optional[numpy.ndarray] = None,
        features_inspect: Optional[List[str]] = None, submodel_search: int = 100,
        processes: Optional[int] = None):
    """
    Estimate average partial importance using cross-validation method specified.

    Since the total size of the submodel space may be huge to sample and calculate,
    the feature importance is calculated by sampling `N=submodel_search` submodels,

    Args:
        estimator: scikit-learn estimator to use.
        data_frame: Training data frame.
        features: List of features to use.
        y_name: Prediction target name.
        cv_method: What cross-validation to use, can be 'loo', int (k-fold).
        weights: Sample weights.
        features_inspect: What features to use when computing importance.
        submodel_search: Number of sub-models when computing importance.
        processes: Number of multiprocessing processes.
    """
    # pylint: disable=too-many-locals

    partial = {}
    features_inspect = features_inspect or features

    data_y = data_frame[y_name].values.astype(float)
    if weights is not None:
        weights = weights.astype(float)

    with Pool(processes=processes or cpu_count()) as pool:
        submodel_to_metric = {}
        fea2submodels = {}
        for fea in tqdm(features_inspect, desc='Dispatch jobs'):
            complement_model = [f for f in features if f != fea]
            submodels = numpy.random.randint(2, size=(submodel_search, len(complement_model)))

            for submodel_coded in submodels:
                submodel = tuple(f for i, f in enumerate(complement_model) if submodel_coded[i] > 0)
                if submodel not in submodel_to_metric:
                    submodel_to_metric[submodel] = pool.apply_async(
                        _fast_cv, (
                            estimator, data_frame[list(submodel)].values.astype(float), data_y,
                            weights, cv_method
                        ))
                inc_submodel = tuple(f for f in features if f in submodel or f == fea)
                if inc_submodel not in submodel_to_metric:
                    submodel_to_metric[inc_submodel] = pool.apply_async(
                        _fast_cv, (
                            estimator, data_frame[list(inc_submodel)].values.astype(float), data_y,
                            weights, cv_method
                        ))
            fea2submodels[fea] = submodels

        for submodel in tqdm(submodel_to_metric, desc='Collecting results'):
            submodel_to_metric[submodel] = submodel_to_metric[submodel].get()

        for fea in tqdm(features_inspect, desc='Computing metrics'):
            complement_model = [f for f in features if f != fea]
            metric_inc = []
            for submodel_coded in fea2submodels[fea]:
                submodel = tuple(f for i, f in enumerate(complement_model) if submodel_coded[i] > 0)
                inc_submodel = tuple(f for f in features if f in submodel or f == fea)
                metric_inc.append(submodel_to_metric[inc_submodel] - submodel_to_metric[submodel])
            partial[fea] = numpy.mean(metric_inc)
    return partial


def model_cv_analysis(  # pylint: disable=too-many-arguments
        data_frame: DataFrame, features: List[str], y_name: str, y_desc: str, *,
        weights: Optional[numpy.ndarray] = None, estimator: RegressorMixin = None,
        do_normalization: bool = False, cv_method: Union[str, int] = 'loo',
        do_plot: bool = True, processes: Optional[int] = None, display_pbar: bool = True):
    """
    Perform cross validation on a dataset using multiprocessing.

    Args:
        estimator: scikit-learn estimator to use.
        data_frame: Training data frame.
        features: List of features to use.
        y_name: Prediction target name.
        y_desc: Description of prediction target, will used in plots.
        do_normalization: Do normalization before running model.
        cv_method: What cross-validation to use, can be 'loo', int (k-fold).
        weights: Sample weights.
        processes: Number of multiprocessing processes.
        do_plot: Whether or not to plot data.
        display_pbar: Whether or not to display progress bar.
    """
    # pylint: disable=too-many-locals

    if estimator is None:
        estimator = RandomForestRegressor({
            'n_estimators': 500,
            'max_depth': 10,
            'max_features': 'sqrt',
            'bootstrap': True,
            'oob_score': True,
        })

    if cv_method == 'loo':
        cv_split = LeaveOneOut()
    elif isinstance(cv_method, int):
        cv_split = KFold(cv_method)
    else:
        raise ValueError(f'Invalid cv_method {cv_method}')

    data_x = data_frame[features].values.astype(float)
    data_y = data_frame[y_name].values.astype(float)

    with Pool(processes=processes or cpu_count(), initializer=_setup_cv,
              initargs=(data_x, data_y, weights, estimator, do_normalization)) as pool:
        # oob_scores = []
        num_splits = cv_split.get_n_splits(data_x)
        train_residual = numpy.full((data_y.size, num_splits), fill_value=float('nan'))
        test_residual = numpy.empty_like(data_y)

        results = pool.imap_unordered(_run_cv, cv_split.split(data_x), chunksize=16)
        if display_pbar:
            results = tqdm(results, total=num_splits, desc='Computing CV')

        for i, result in enumerate(results):
            test_ind = result['Test indices']
            test_residual[test_ind] = result['Test residual']

            train_ind = result['Train indices']
            train_residual[train_ind, i] = result['Train residual']
            # oob_scores.append(result['OOB score'])

        train_predictions = data_y - numpy.nanmean(train_residual, axis=1)
        cv_predictions = data_y - test_residual
        train_score = r2_score(data_y, train_predictions, sample_weight=weights)
        cv_score = r2_score(data_y, cv_predictions, sample_weight=weights)

        # print('Average OOB score:', numpy.mean(oob_scores))
        # print('Cross-validated score:', cv_score)

    if do_plot:
        _, (ax1, ax2) = plt.subplots(ncols=2, figsize=(10, 4))
        ax1.scatter(data_frame[y_name], cv_predictions, color=weights2colors(weights))
        ax1.set_xlabel('Observed ' + y_desc)
        ax1.set_ylabel('Predicted ' + y_desc)

        ax2.hist(test_residual, bins=50)
        ax2.set_xlabel('Residual')
        plt.show()

    return {
        # 'Average OOB score': numpy.mean(oob_scores),
        'Train score': train_score,
        'CV score': cv_score,
        'CV predictions': cv_predictions,
    }


def custom_grid_search(  # pylint: disable=too-many-arguments
        estimator_cls, params, data_frame, features, y_name, weights=None,
        do_normalization=False, cv_method='loo', processes=None):
    """
    Perform Grid search CV and find the best parameters.

    :param estimator_cls: scikit-learn learner class.
    :param params: Parameters to use for `estimator_cls`.
    :param data_frame: Training data frame.
    :param features: Features to use.
    :param y_name: Prediction target variable name.
    :param weights: Sample weights.
    :param do_normalization: Whether to do normalization for training features.
    :param cv_method: Cross-validation method.
    :param processes: Number of multiprocessing processes to use.
    """
    # pylint: disable=too-many-locals
    param_names = list(params)
    param_values = list(product(*[params[x] for x in param_names]))
    param_cv_results = {}

    for param_value in param_values:
        params = dict(zip(param_names, param_value))
        estimator = estimator_cls(**params)
        print('Estimator params:', params, end=' ==> ')
        cv_result = model_cv_analysis(
            data_frame=data_frame,
            features=features,
            y_name=y_name, y_desc='',
            estimator=estimator,
            do_normalization=do_normalization,
            cv_method=cv_method,
            weights=weights,
            do_plot=False,
            display_pbar=False,
            processes=processes,
        )
        print('CV score:', '%.3f' % cv_result['CV score'])
        param_cv_results[param_value] = cv_result['CV score']

    best_params = max(param_cv_results, key=lambda x: param_cv_results[x])

    return {
        'CV results': param_cv_results,
        'Best param': best_params,
        'Best CV result': param_cv_results[best_params]
    }
