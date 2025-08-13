#!/usr/bin/env python3
"""
Two learning pipelines in one file:

1. Simple tabular RF/GBDT (scikit‑learn)
2. 1‑D CNN that takes a one‑hot 20×4 spacer as input

The model outputs a prediction of the editing efficiency (0‑1).
"""

import numpy as np
import pandas as pd
from pathlib import Path
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score, mean_absolute_error
import tensorflow as tf
from tensorflow import layers, models, callbacks

# ------------------------------
# Utility to load data
# ------------------------------
def load_data():
    feat_path = Path("data/features/SlAREB1_features.csv")
    spacer_path = Path("data/guides/SlAREB1_guides.csv")  # contains the original spacer
    feat_df = pd.read_csv(feat_path)
    spacer_df = pd.read_csv(spacer_path)
    # one‑hot encode
    onehot = np.array([np.eye(4)[list("ACGT".index(c) for c in s)] for s in spacer_df["spacer"]])
    # target – use column “score” from guide generator
    target = spacer_df["score"].values
    return feat_df, onehot, target

# ------------------------------
# 1. Tabular model
# ------------------------------
def train_tabular():
    X, onehot, y = load_data()
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=.2, random_state=42)

    rf = RandomForestRegressor(n_estimators=300, random_state=42, n_jobs=-1)
    rf.fit(X_train, y_train)
    preds = rf.predict(X_test)
    print("Random Forest R2:", r2_score(y_test, preds))
    print("MAE:", mean_absolute_error(y_test, preds))
    # Save
    Path("models/rf.pkl").write_bytes(pickle.dumps(rf))

# ------------------------------
# 2. CNN model
# ------------------------------
def build_cnn():
    inputs = layers.Input(shape=(20,4))          # 20‑mer, 4‑nt bases
    x = layers.Conv1D(64, 3, padding='same', activation='relu')(inputs)
    x = layers.BatchNormalization()(x)
    x = layers.Conv1D(128, 3, padding='same', activation='relu')(x)
    x = layers.GlobalMaxPooling1D()(x)
    x = layers.Dense(64, activation='relu')(x)
    outputs = layers.Dense(1, activation='linear')(x)  # regression
    model = models.Model(inputs, outputs)
    model.compile(optimizer='adam', loss='mse')
    return model

def train_cnn():
    X_tab, onehot, y = load_data()
    X_train, X_test, y_train, y_test = train_test_split(onehot, y, test_size=.2, random_state=42)

    model = build_cnn()
    es = callbacks.EarlyStopping(monitor='val_loss', patience=4, restore_best_weights=True)
    history = model.fit(X_train, y_train, validation_split=.2,
                        batch_size=128, epochs=50, callbacks=[es], verbose=1)
    # Evaluation
    preds = model.predict(X_test).flatten()
    print("CNN R2:", r2_score(y_test, preds))
    print("CNN MAE:", mean_absolute_error(y_test, preds))
    # Save
    model.save("models/cnn_slareb1.h5")

# ------------------------------
# Main
# ------------------------------
if __name__ == "__main__":
    import argparse, pickle
    parser = argparse.ArgumentParser()
    parser.add_argument("--train", choices=["rf","cnn","both"], default="both")
    args = parser.parse_args()
    Path("models").mkdir(parents=True, exist_ok=True)
    if args.train in ("rf","both"):
        train_tabular()
    if args.train in ("cnn","both"):
        train_cnn()
