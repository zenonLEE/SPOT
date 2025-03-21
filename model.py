from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor, StackingRegressor
from sklearn.svm import SVR
from sklearn.neural_network import MLPRegressor
from sklearn.linear_model import Ridge, Lasso, ElasticNet, BayesianRidge
from catboost import CatBoostRegressor
from sklearn.neighbors import KNeighborsRegressor
from lightgbm import LGBMRegressor
from xgboost import XGBRegressor
from tensorflow import keras
from tensorflow.keras.layers import BatchNormalization
from tensorflow.keras.regularizers import l2
from sklearn.metrics import mean_squared_error
import optuna
import numpy as np
from sklearn.model_selection import KFold

from sklearn.model_selection import cross_val_score

def get_model(model_name, random_state=42, X_train=None, y_train=None, sample_weight=None, **model_params):
    """
    Returns a configured regression model instance, supporting custom parameters via **model_params.

    Parameters:
    - model_name : str - Model name (case-insensitive)
    - random_state : int - Random seed
    - model_params : dict - Model-specific parameters (overriding defaults)

    Supported models:
    'stacking', 'catboost', 'deep_neural', 'randomforest', 'gradientboosting',
    'xgboost', 'xgboost_optuna', 'lightgbm', 'svr', 'mlp', 'ridge', 'lasso', 'elasticnet',
    'bayesianridge', 'knn'
    """
    model_name = model_name.lower()

    # Function to apply parameters selectively
    def apply_params(model_class, default_params, **extra_params):
        params = {**default_params, **extra_params}
        return model_class(**params)

    if model_name == 'stacking':
        estimators = [
            ('rf', RandomForestRegressor(
                n_estimators=200,
                max_depth=12,
                max_samples=0.8,
                random_state=random_state,
                n_jobs=-1
            )),
            ('gb', GradientBoostingRegressor(
                n_estimators=300,
                learning_rate=0.05,
                max_depth=5,
                subsample=0.8,
                random_state=random_state
            )),
            ('cat', CatBoostRegressor(
                iterations=800,
                learning_rate=0.05,
                depth=8,
                l2_leaf_reg=5,
                verbose=False,
                random_state=random_state
            )),
            ('svr', SVR(
                C=10,
                epsilon=0.1,
                kernel='rbf'
            ))
        ]

        final_estimator = XGBRegressor(
            n_estimators=200,
            learning_rate=0.01,
            max_depth=4,
            subsample=0.8,
            colsample_bytree=0.8,
            reg_alpha=0.1,
            reg_lambda=0.1,
            random_state=random_state,
            n_jobs=-1
        )

        return StackingRegressor(
            estimators=estimators,
            final_estimator=final_estimator,
            cv=5,
            passthrough=False,
            **model_params
        )

    elif model_name == 'catboost':
        default_params = {
            'iterations': 1000,
            'learning_rate': 0.05,
            'depth': 8,
            'l2_leaf_reg': 5,
            'early_stopping_rounds': 30,
            'verbose': False,
            'random_state': random_state
        }
        return apply_params(CatBoostRegressor, default_params, **model_params)

    elif model_name == 'deep_neural':
        # Deep neural networks typically have different interfaces; consider handling separately
        model = keras.Sequential([
            keras.layers.Dense(256, activation='relu', kernel_initializer='he_normal',
                               input_shape=(6,), kernel_regularizer=l2(0.01)),
            BatchNormalization(),
            keras.layers.Dropout(0.4),

            keras.layers.Dense(128, activation='relu', kernel_regularizer=l2(0.005)),
            BatchNormalization(),
            keras.layers.Dropout(0.3),

            keras.layers.Dense(64, activation='relu'),
            BatchNormalization(),

            keras.layers.Dense(1, kernel_regularizer=l2(0.001))
        ])

        # Optimizer configuration
        optimizer = keras.optimizers.Adam(
            learning_rate=0.001,
            beta_1=0.9,
            beta_2=0.999,
            epsilon=1e-07
        )

        model.compile(
            optimizer=optimizer,
            loss='mse',
            metrics=['mae', keras.metrics.RootMeanSquaredError()]
        )
        return model

    elif model_name == 'randomforest':
        default_params = {
            'n_estimators': 300,
            'max_depth': 12,
            'max_features': 'sqrt',
            'min_samples_split': 5,
            'random_state': random_state,
            'n_jobs': -1
        }
        return apply_params(RandomForestRegressor, default_params, **model_params)

    elif model_name == 'gradientboosting':
        default_params = {
            'n_estimators': 500,
            'learning_rate': 0.05,
            'max_depth': 5,
            'subsample': 0.8,
            'min_samples_split': 5,
            'random_state': random_state
        }
        return apply_params(GradientBoostingRegressor, default_params, **model_params)

    elif model_name == 'xgboost':
        default_params = {
            'n_estimators': 500,
            'learning_rate': 0.05,
            'max_depth': 6,
            'subsample': 0.8,
            'colsample_bytree': 0.8,
            'reg_alpha': 0.1,
            'reg_lambda': 0.1,
            'random_state': random_state,
            'n_jobs': -1
        }
        return apply_params(XGBRegressor, default_params, **model_params)



    elif model_name == 'mlp':
        default_params = {
            'hidden_layer_sizes': (128, 64),
            'activation': 'tanh',
            'solver': 'adam',
            'alpha': 0.0001,
            'batch_size': 128,
            'learning_rate': 'adaptive',
            'learning_rate_init': 0.01,
            'early_stopping': True,
            'validation_fraction': 0.2,
            'max_iter': 2000,
            'random_state': random_state,
            'momentum': 0.9,
        }
        return apply_params(MLPRegressor, default_params, **model_params)

    elif model_name == 'svr':
        default_params = {
            'C': 10,
            'epsilon': 0.1,
            'kernel': 'rbf',
            'gamma': 'scale'
        }
        return apply_params(SVR, default_params, **model_params)

    elif model_name == 'ridge':
        default_params = {
            'alpha': 0.5,
            'random_state': random_state
        }
        return apply_params(Ridge, default_params, **model_params)

    elif model_name == 'lasso':
        default_params = {
            'alpha': 0.01,
            'random_state': random_state
        }
        return apply_params(Lasso, default_params, **model_params)

    elif model_name == 'elasticnet':
        default_params = {
            'alpha': 0.01,
            'l1_ratio': 0.5,
            'random_state': random_state
        }
        return apply_params(ElasticNet, default_params, **model_params)

    elif model_name == 'bayesianridge':
        default_params = {
            'random_state': random_state
        }
        return apply_params(BayesianRidge, default_params, **model_params)

    elif model_name == 'knn':
        default_params = {
            'n_neighbors': 7,
            'weights': 'distance',
            'algorithm': 'auto',
            'n_jobs': -1
        }
        default_params.pop('random_state', None)
        return apply_params(KNeighborsRegressor, default_params, **model_params)

    elif model_name == 'xgboost_optuna':
        if X_train is None or y_train is None:
            raise ValueError("X_train and y_train must be provided for Optuna optimization")

        def objective(trial):
            params = {
                'n_estimators': trial.suggest_int('n_estimators', 100, 1000),
                'learning_rate': trial.suggest_float('learning_rate', 0.001, 0.1, log=True),
                'max_depth': trial.suggest_int('max_depth', 3, 10),
                'subsample': trial.suggest_float('subsample', 0.5, 1.0),
                'colsample_bytree': trial.suggest_float('colsample_bytree', 0.5, 1.0),
                'min_child_weight': trial.suggest_int('min_child_weight', 1, 7),
                'reg_alpha': trial.suggest_float('reg_alpha', 1e-8, 1.0, log=True),
                'reg_lambda': trial.suggest_float('reg_lambda', 1e-8, 1.0, log=True),
                'random_state': random_state,
                'n_jobs': -1
            }

            # Apply any user-provided parameters
            trial_params = {**params, **model_params}
            model = XGBRegressor(**trial_params)

            # 使用传入的 train_model 中的 KFold 结果
            kf = KFold(n_splits=5, shuffle=True, random_state=random_state)
            scores = []

            for train_idx, val_idx in kf.split(X_train):
                X_fold_train, X_fold_val = X_train[train_idx], X_train[val_idx]
                y_fold_train, y_fold_val = y_train[train_idx], y_train[val_idx]

                # 使用样本权重（如果提供）
                if sample_weight is not None:
                    sw_fold = sample_weight[train_idx]
                    model.fit(
                        X_fold_train, y_fold_train,
                        sample_weight=sw_fold,
                        eval_set=[(X_fold_val, y_fold_val)],
                        early_stopping_rounds=50,
                        verbose=False
                    )
                else:
                    model.fit(
                        X_fold_train, y_fold_train,
                        eval_set=[(X_fold_val, y_fold_val)],
                        early_stopping_rounds=50,
                        verbose=False
                    )

                y_pred = model.predict(X_fold_val)
                rmse = np.sqrt(mean_squared_error(y_fold_val, y_pred))
                scores.append(rmse)

            return np.mean(scores)

        study = optuna.create_study(
            direction='minimize',
            sampler=optuna.samplers.TPESampler(seed=random_state)
        )

        study.optimize(
            objective,
            n_trials=50,
            show_progress_bar=True
        )

        best_params = study.best_trial.params
        best_params.update({
            'random_state': random_state,
            'n_jobs': -1,
         #   'early_stopping_rounds': 50,
            'verbose': False
        })

        final_params = {**best_params, **model_params}
        return XGBRegressor(**final_params)

    elif model_name == 'xgboost_optuna':

        def objective(trial):
            params = {
                'n_estimators': trial.suggest_int('n_estimators', 100, 1000),
                'learning_rate': trial.suggest_float('learning_rate', 0.001, 0.1, log=True),
                'max_depth': trial.suggest_int('max_depth', 3, 10),
                'subsample': trial.suggest_float('subsample', 0.5, 1.0),
                'colsample_bytree': trial.suggest_float('colsample_bytree', 0.5, 1.0),
                'min_child_weight': trial.suggest_int('min_child_weight', 1, 7),
                'reg_alpha': trial.suggest_float('reg_alpha', 1e-8, 1.0, log=True),
                'reg_lambda': trial.suggest_float('reg_lambda', 1e-8, 1.0, log=True),
                'random_state': random_state,
                'n_jobs': -1
            }
            # 合并任何用户传入的额外参数
            trial_params = {**params, **model_params}
            model = XGBRegressor(**trial_params)
            # 使用 3 折交叉验证评估 RMSE
            kf = KFold(n_splits=5, shuffle=True, random_state=random_state)
            scores = []
            for train_idx, val_idx in kf.split(X_train):
                X_fold_train, X_fold_val = X_train[train_idx], X_train[val_idx]
                y_fold_train, y_fold_val = y_train[train_idx], y_train[val_idx]
                if 'sample_weight' in model_params and model_params['sample_weight'] is not None:
                    sw_fold = model_params['sample_weight'][train_idx]
                    model.fit(
                        X_fold_train, y_fold_train,
                        sample_weight=sw_fold,
                        eval_set=[(X_fold_val, y_fold_val)],
                        early_stopping_rounds=50,
                        verbose=False
                    )
                else:
                    model.fit(
                        X_fold_train, y_fold_train,
                        eval_set=[(X_fold_val, y_fold_val)],
                        early_stopping_rounds=50,
                        verbose=False
                    )
                y_pred = model.predict(X_fold_val)
                rmse = np.sqrt(mean_squared_error(y_fold_val, y_pred))
                scores.append(rmse)
            return np.mean(scores)

        study = optuna.create_study(
            direction='minimize',
            sampler=optuna.samplers.TPESampler(seed=random_state)
        )
        study.optimize(objective, n_trials=50, show_progress_bar=True)
        best_params = study.best_trial.params
        best_params.update({
            'random_state': random_state,
            'n_jobs': -1,
            'verbose': False
        })
        # 与 xgboost 分支对齐，使用 apply_params 构造模型实例
        return apply_params(XGBRegressor, best_params, **model_params)
