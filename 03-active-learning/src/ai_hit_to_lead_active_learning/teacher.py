"""
Teacher ensemble made of Random Forest and XGBoost regressors.
"""

from __future__ import annotations

from typing import Dict, Tuple

import numpy as np
from sklearn.ensemble import RandomForestRegressor
from xgboost import XGBRegressor


class TeacherEnsemble:
    def __init__(self) -> None:
        self.models: Dict[str, object] = {
            "rf_1": RandomForestRegressor(n_estimators=200, max_depth=10, random_state=42, n_jobs=-1),
            "rf_2": RandomForestRegressor(n_estimators=200, max_depth=15, random_state=43, n_jobs=-1),
            "xgb_1": XGBRegressor(
                n_estimators=200,
                max_depth=6,
                learning_rate=0.1,
                random_state=42,
                n_jobs=-1,
                verbosity=0,
            ),
            "xgb_2": XGBRegressor(
                n_estimators=200,
                max_depth=8,
                learning_rate=0.05,
                random_state=43,
                n_jobs=-1,
                verbosity=0,
            ),
        }
        self.is_fitted = False

    def fit(self, X: np.ndarray, y: np.ndarray) -> "TeacherEnsemble":
        for model in self.models.values():
            model.fit(X, y)
        self.is_fitted = True
        return self

    def predict(self, X: np.ndarray) -> np.ndarray:
        if not self.is_fitted:
            raise RuntimeError("TeacherEnsemble must be fitted before predict()")
        predictions = np.array([model.predict(X) for model in self.models.values()])
        return predictions.mean(axis=0)

    def predict_with_std(self, X: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        if not self.is_fitted:
            raise RuntimeError("TeacherEnsemble must be fitted before predict_with_std()")
        predictions = np.array([model.predict(X) for model in self.models.values()])
        return predictions.mean(axis=0), predictions.std(axis=0)
