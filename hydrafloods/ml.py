import ee
import copy
import numpy as np
import pandas as pd
from hydrafloods import geeutils, decorators


def fcnn(img, out_band_bames=None, probabilities=True, **kwargs):
    @decorators.carry_metadata
    def _apply(image):
        # run the predictions
        prediction = model.predictImage(image.toFloat().toArray())

        if probabilities is False:
            # find highest probability class
            precition = prediction.arrayArgmax().arrayFlatten([outBandNames])
        else:
            # flatten probability array to image with bands
            prediction = prediction.arrayFlatten([outBandNames]).toFloat()

        return prediction

    # Load the trained model and use it for prediction.
    model = ee.Model.fromAiPlatformPredictor(**kwargs)

    # run predictions over image collection and
    return _apply(img)


def logistic_regression():

    return


def fc_to_scaler(fc):
    return

def random_forest_ee(
    n_trees,
    feature_collection,
    feature_names,
    label,
    normalize_features=True,
    mode="classification",
):
    def feature_scaling(feature):
        def iter_cols(i):
            i = ee.String(i)
            v = ee.Number(feature.get(i))
            minv = ee.Number(min_max_dict.get(i.cat("_min")))
            maxv = ee.Number(min_max_dict.get(i.cat("_max")))
            return v.subtract(minv).divide(maxv.subtract(minv))

        scaled = ee_feature_names.map(iter_cols)
        newVals = ee.Dictionary.fromLists(ee_feature_names, scaled)
        return feature.set(newVals)

    fc = ee.FeatureCollection(feature_collection)
    ee_feature_names = ee.List(feature_names)

    if normalize_features:
        feature_min_max = fc.reduceColumns(
            ee.Reducer.minMax().repeat(ee_feature_names.length()), ee_feature_names
        )

        # min/max feature names
        names = ee_feature_names.map(
            lambda x: ee.List([ee.String(x).cat("_min"), ee.String(x).cat("_max")])
        ).flatten()

        # get the min/max values for each feature
        # used to scale values from 0-1
        min_max_dict = ee.Dictionary.fromLists(
            names,
            ee.List(feature_min_max.get("min"))
            .zip(feature_min_max.get("max"))
            .flatten(),
        )

        # normalize the features in the entire featureCollection
        fc_norm = fc.map(feature_scaling)

    else:
        fc_norm = fc

    classifier = (
        ee.Classifier.smileRandomForest(n_trees)
        .setOutputMode(mode.upper())
        .train(fc_norm, label, ee_feature_names)
    )

    return classifier, min_max_dict


def random_forest_from_str(path):
    if sample_path.startswith("gs://"):
        estimator_files = ee.List(utils.list_gcs_objs(path, pattern="estimator*.txt"))

        forest_strings = estimator_files.map(lambda x: ee.Blob(x).string())

        classifier = ee.Classifier.decisionTreeEnsemble(forest_strings)

    else:
        raise NotImplementedError("currently only reading files from GCS is supported")

    return classifier


def optimize(model, features, labels, opt_params=None):
    from bayes_opt import BayesianOptimization

    def crossval(n_estimators, min_samples_split, max_features):
        """Wrapper of RandomForest cross validation.
        Notice how we ensure n_estimators and min_samples_split are casted
        to integer before we pass them along. Moreover, to avoid max_features
        taking values outside the (0, 1) range, we also ensure it is capped
        accordingly.
        """
        return cv(
            n_estimators=int(n_estimators),
            min_samples_split=int(min_samples_split),
            max_features=max(min(max_features, 0.999), 1e-3),
            data=data,
            targets=targets,
        )

    optimizer = BayesianOptimization(
        f=crossval,
        pbounds={
            "n_estimators": (10, 250),
            "min_samples_split": (2, 25),
            "max_features": (0.1, 0.999),
        },
        random_state=1234,
        verbose=2,
    )
    optimizer.maximize(n_iter=10)

    print("Final result:", optimizer.max)

    return


def cv(model, features, labels, n_cv=5, classifier=False, **kwargs):
    """Random Forest cross validation.
    This function will instantiate a random forest classifier with parameters
    n_estimators, min_samples_split, and max_features. Combined with data and
    targets this will in turn be used to perform cross validation. The result
    of cross validation is returned.
    Our goal is to find combinations of n_estimators, min_samples_split, and
    max_features that minimzes the log loss.
    """

    score_metric = "neg_log_loss" if classifier else "mean_squared_error"
    estimator = model(**kwargs)
    cval = cross_val_score(estimator, features, labels, scoring=score_metric, cv=n_cv)
    return cval.mean()


def sklearn_tree_to_string(estimator, feature_names):

    # extract out the information need to build the tree string
    n_nodes = estimator.tree_.node_count
    children_left = estimator.tree_.children_left
    children_right = estimator.tree_.children_right
    feature_idx = estimator.tree_.feature
    values = estimator.tree_.value
    impurities = estimator.tree_.impurity
    n_samples = estimator.tree_.n_node_samples
    thresholds = estimator.tree_.threshold
    features = [feature_names[i] for i in feature_idx]

    # use iterative pre-order search to extract node depth and leaf information
    node_ids = np.zeros(shape=n_nodes, dtype=np.int64)
    node_depth = np.zeros(shape=n_nodes, dtype=np.int64)
    is_leaves = np.zeros(shape=n_nodes, dtype=bool)
    stack = [(0, -1)]  # seed is the root node id and its parent depth
    while len(stack) > 0:
        node_id, parent_depth = stack.pop()
        node_depth[node_id] = parent_depth + 1
        node_ids[node_id] = node_id

        # If we have a test node
        if children_left[node_id] != children_right[node_id]:
            stack.append((children_left[node_id], parent_depth + 1))
            stack.append((children_right[node_id], parent_depth + 1))
        else:
            is_leaves[node_id] = True

    # create a table of the initial structure
    # each row is a node or leaf
    df = pd.DataFrame(
        {
            "node_id": node_ids,
            "node_depth": node_depth,
            "is_leaf": is_leaves,
            "children_left": children_left,
            "children_right": children_right,
            "value": np.squeeze(values),
            "criterion": impurities,
            "n_samples": n_samples,
            "threshold": thresholds,
            "feature_name": features,
            "sign": ["<="] * n_nodes,
        }
    )

    # the table representation does not have lef vs right node structure
    # so we need to add in right nodes in the correct location
    # we do this by first calculating which nodes are right and then insert them at the correct index

    # get a dict of right node rows and assign key based on index where to insert
    inserts = {}
    for row in df.itertuples():
        child_r = row.children_right
        if child_r > row.Index:
            ordered_row = np.array(row)
            ordered_row[-1] = ">"
            inserts[child_r] = ordered_row[1:]  # drop index value
    # sort the inserts as to keep track of the additive indexing
    inserts_sorted = {k: inserts[k] for k in sorted(inserts.keys())}

    # loop through the row inserts and add to table (array)
    table_values = df.values
    for i, k in enumerate(inserts_sorted.keys()):
        table_values = np.insert(table_values, (k + i), inserts_sorted[k], axis=0)

    # make the ordered table array into a dataframe
    # note: df is dtype "object", need to cast later on
    ordered_df = pd.DataFrame(table_values, columns=df.columns)

    max_depth = np.max(ordered_df.node_depth.astype(int))
    tree_str = f"1) root 0 0 0 (0)\n"
    previous_depth = -1
    cnts = []
    # loop through the nodes and calculate the node number and values per node
    for row in ordered_df.itertuples():
        node_depth = int(row.node_depth)
        left = int(row.children_left)
        right = int(row.children_right)
        if left != right:
            if row.Index == 0:
                cnt = 2
            elif previous_depth > node_depth:
                depths = ordered_df.node_depth.values[: row.Index]
                idx = np.where(depths == node_depth)[0][-1]
                # cnt = (cnts[row.Index-1] // 2) + 1
                cnt = cnts[idx] + 1
            elif previous_depth < node_depth:
                cnt = cnts[row.Index - 1] * 2
            elif previous_depth == node_depth:
                cnt = cnts[row.Index - 1] + 1

            if node_depth == (max_depth - 1):
                value = float(ordered_df.iloc[row.Index + 1].value)
                samps = int(ordered_df.iloc[row.Index + 1].n_samples)
                criterion = float(ordered_df.iloc[row.Index + 1].criterion)
                tail = " *\n"
            else:
                if (
                    (bool(ordered_df.loc[ordered_df.node_id == left].iloc[0].is_leaf))
                    and (
                        bool(
                            int(row.Index)
                            < int(ordered_df.loc[ordered_df.node_id == left].index[0])
                        )
                    )
                    and (str(row.sign) == "<=")
                ):
                    rowx = ordered_df.loc[ordered_df.node_id == left].iloc[0]
                    tail = " *\n"
                    value = float(rowx.value)
                    samps = int(rowx.n_samples)
                    criterion = float(rowx.criterion)

                elif (
                    (bool(ordered_df.loc[ordered_df.node_id == right].iloc[0].is_leaf))
                    and (
                        bool(
                            int(row.Index)
                            < int(ordered_df.loc[ordered_df.node_id == right].index[0])
                        )
                    )
                    and (str(row.sign) == ">")
                ):
                    rowx = ordered_df.loc[ordered_df.node_id == right].iloc[0]
                    tail = " *\n"
                    value = float(rowx.value)
                    samps = int(rowx.n_samples)
                    criterion = float(rowx.criterion)

                else:
                    value = float(row.value)
                    samps = 0  # int(row.n_samples)
                    criterion = 0  # float(row.criterion)
                    tail = "\n"

            # extract out the information needed in each line
            # spacing = (node_depth + 1) * "  "  # for pretty printing
            fname = str(row.feature_name)  # name of the feature (i.e. band name)
            tresh = float(row.threshold)  # threshold
            sign = str(row.sign)

            # tree_str += f"{cnt}) {fname} {sign} {tresh:.6f} {samps} {criterion:.4f} {value:.6f}{tail}"
            tree_str += f"{cnt}) {fname} {sign} {tresh:.4f} 0 0 {value:.4f}{tail}"
            previous_depth = node_depth
        cnts.append(cnt)

    return tree_str

