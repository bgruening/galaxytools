"""
Predict tool usage to weigh the predicted tools
"""

import collections
import csv
import os
import warnings

import numpy as np
import utils
from sklearn.model_selection import GridSearchCV
from sklearn.pipeline import Pipeline
from sklearn.svm import SVR

warnings.filterwarnings("ignore")

main_path = os.getcwd()


class ToolPopularity:
    def __init__(self):
        """ Init method. """

    def extract_tool_usage(self, tool_usage_file, cutoff_date, dictionary):
        """
        Extract the tool usage over time for each tool
        """
        tool_usage_dict = dict()
        all_dates = list()
        all_tool_list = list(dictionary.keys())
        with open(tool_usage_file, "rt") as usage_file:
            tool_usage = csv.reader(usage_file, delimiter="\t")
            for index, row in enumerate(tool_usage):
                row = [item.strip() for item in row]
                if (str(row[1]).strip() > cutoff_date) is True:
                    tool_id = utils.format_tool_id(row[0])
                    if tool_id in all_tool_list:
                        all_dates.append(row[1])
                        if tool_id not in tool_usage_dict:
                            tool_usage_dict[tool_id] = dict()
                            tool_usage_dict[tool_id][row[1]] = int(row[2])
                        else:
                            curr_date = row[1]
                            # merge the usage of different version of tools into one
                            if curr_date in tool_usage_dict[tool_id]:
                                tool_usage_dict[tool_id][curr_date] += int(row[2])
                            else:
                                tool_usage_dict[tool_id][curr_date] = int(row[2])
        # get unique dates
        unique_dates = list(set(all_dates))
        for tool in tool_usage_dict:
            usage = tool_usage_dict[tool]
            # extract those dates for which tool's usage is not present in raw data
            dates_not_present = list(set(unique_dates) ^ set(usage.keys()))
            # impute the missing values by 0
            for dt in dates_not_present:
                tool_usage_dict[tool][dt] = 0
            # sort the usage list by date
            tool_usage_dict[tool] = collections.OrderedDict(sorted(usage.items()))
        return tool_usage_dict

    def learn_tool_popularity(self, x_reshaped, y_reshaped):
        """
        Fit a curve for the tool usage over time to predict future tool usage
        """
        epsilon = 0.0
        cv = 5
        s_typ = "neg_mean_absolute_error"
        n_jobs = 4
        s_error = 1
        tr_score = False
        try:
            pipe = Pipeline(steps=[("regressor", SVR(gamma="scale"))])
            param_grid = {
                "regressor__kernel": ["rbf", "poly", "linear"],
                "regressor__degree": [2, 3],
            }
            search = GridSearchCV(
                pipe,
                param_grid,
                cv=cv,
                scoring=s_typ,
                n_jobs=n_jobs,
                error_score=s_error,
                return_train_score=tr_score,
            )
            search.fit(x_reshaped, y_reshaped.ravel())
            model = search.best_estimator_
            # set the next time point to get prediction for
            prediction_point = np.reshape([x_reshaped[-1][0] + 1], (1, 1))
            prediction = model.predict(prediction_point)
            if prediction < epsilon:
                prediction = [epsilon]
            return prediction[0]
        except Exception as e:
            print(e)
            return epsilon

    def get_pupularity_prediction(self, tools_usage):
        """
        Get the popularity prediction for each tool
        """
        usage_prediction = dict()
        for tool_name, usage in tools_usage.items():
            y_val = list()
            x_val = list()
            for x, y in usage.items():
                x_val.append(x)
                y_val.append(y)
            x_pos = np.arange(len(x_val))
            x_reshaped = x_pos.reshape(len(x_pos), 1)
            y_reshaped = np.reshape(y_val, (len(x_pos), 1))
            prediction = np.round(self.learn_tool_popularity(x_reshaped, y_reshaped), 8)
            usage_prediction[tool_name] = prediction
        return usage_prediction
