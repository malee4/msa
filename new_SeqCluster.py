# AUTHOR:Jimmy Ka Ho Chiu & Rick Twee-Hee Ong
# CODE VERSION: 2022, last updated locally June 21, 2022
# LANGUAGE: Python
# SOURCE: "Clustering biological sequences with dynamic sequence similarity threshold"
# URL: https://doi.org/10.1186/s12859-022-04643-9


class SeqCluster:
    # default parameters (unfilled)
    _res_param_start = None
    _res_param_end = None
    _res_param_step_size = None
    _precision = None
    _seed = None
    _is_verbose = None
    _init = False

    @classmethod
    def init(cls, user_params, is_verbose=True):
        cls._res_param_start = user_params.res_param_start
        cls._res_param_end = user_params.res_param_end
        cls._res_param_step_size = user_params.res_param_step_size
        cls._precision = user_params.precision
        cls._seed = user_params.seed
        cls._is_verbose = is_verbose
        cls._init = True
