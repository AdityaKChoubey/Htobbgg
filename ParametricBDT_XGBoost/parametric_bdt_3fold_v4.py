import os
import re
import glob
import uproot
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score, roc_curve
from xgboost import XGBClassifier

#configurations

SIGNAL_DIR = "."
SIGNAL_PATTERN = "curated_signal_m*.root"
BACKGROUND_FILE = "curated_ttbar_semileptonic.root"
TREE_NAME = "Events"
OUTDIR = "xgb_outputs_parametric_3fold"

os.makedirs(OUTDIR, exist_ok=True)

plt.rcParams.update({
    "figure.figsize": (7, 6),
    "font.size": 12,
    "axes.grid": True
})

#Models

MODEL_CONFIGS = {
    "model1":  dict(n_estimators=500,  max_depth=4, learning_rate=0.05, subsample=0.8, colsample_bytree=0.8),
    "model2":  dict(n_estimators=800,  max_depth=6, learning_rate=0.05, subsample=0.8, colsample_bytree=0.8),
    "model3":  dict(n_estimators=600,  max_depth=4, learning_rate=0.03, subsample=0.7, colsample_bytree=0.7, reg_alpha=1.0, reg_lambda=2.0),
    "model4":  dict(n_estimators=400,  max_depth=2, learning_rate=0.1,  subsample=0.8, colsample_bytree=0.8),
    "model5":  dict(n_estimators=1200, max_depth=3, learning_rate=0.03, subsample=0.8, colsample_bytree=0.8),
    "model6":  dict(n_estimators=1200, max_depth=8, learning_rate=0.1,  subsample=1.0, colsample_bytree=1.0),
    "model7":  dict(n_estimators=700,  max_depth=4, learning_rate=0.03, subsample=0.6, colsample_bytree=0.6, reg_alpha=5.0, reg_lambda=5.0),
    "model8":  dict(n_estimators=600,  max_depth=5, learning_rate=0.05, subsample=0.8, colsample_bytree=0.8, min_child_weight=5),
    "model9":  dict(n_estimators=600,  max_depth=4, learning_rate=0.05, subsample=0.8, colsample_bytree=0.8, gamma=2.0),
    "model10": dict(n_estimators=1500, max_depth=3, learning_rate=0.02, subsample=0.7, colsample_bytree=0.7, reg_alpha=1.0, reg_lambda=3.0),
}

FEATURES = [
    "dphi_gg", "dphi_bb", "dphi_bbgg", "dphi_gg_met",
    "b1_pt", "b2_pt", "g1_pt", "g2_pt",
    "b1_eta", "b2_eta", "lep_pt",
    "theta"
]

RANDOM_STATE = 42

BASE_PARAMS = dict(
    objective="binary:logistic",
    eval_metric="auc",
    tree_method="hist",
    random_state=RANDOM_STATE,
    n_jobs=8,
)

#helper functions

def extract_mass(path):
    m = re.search(r"m(\d+)", path)
    return int(m.group(1)) if m else None

def load_df(file):
    with uproot.open(file) as f:
        return f[TREE_NAME].arrays(FEATURES[:-1], library="pd")

def clean(df):
    return df.replace([np.inf, -np.inf], np.nan).dropna()

def avg_roc(roc_list):
    """Average a list of (fpr, tpr) curves onto a common FPR grid."""
    mean_fpr = np.linspace(0, 1, 500)
    tprs = []
    for fpr, tpr in roc_list:
        f = interp1d(fpr, tpr, bounds_error=False, fill_value=(0.0, 1.0))
        tprs.append(f(mean_fpr))
    return mean_fpr, np.mean(tprs, axis=0)

#Data Loading

signal_files = glob.glob(os.path.join(SIGNAL_DIR, SIGNAL_PATTERN))
signal_map = {extract_mass(f): f for f in signal_files}

mass_points = sorted(signal_map.keys())
print("Mass points:", mass_points)

sig_all = pd.concat([
    clean(load_df(signal_map[m])).assign(label=1, theta=m)
    for m in mass_points
], ignore_index=True)

bkg_df = clean(load_df(BACKGROUND_FILE))
bkg_df["label"] = 0
bkg_df["theta"] = np.random.choice(mass_points, size=len(bkg_df))

n = min(len(sig_all), len(bkg_df))
sig_all = sig_all.sample(n, random_state=RANDOM_STATE)
bkg_df  = bkg_df.sample(n,  random_state=RANDOM_STATE)

df_all = pd.concat([sig_all, bkg_df], ignore_index=True)

X = df_all[FEATURES].astype(np.float32)
y = df_all["label"].astype(int)

kf = StratifiedKFold(n_splits=3, shuffle=True, random_state=RANDOM_STATE)

#Training 

mass_auc_records   = []
all_model_mass_auc = {}

for name, params in MODEL_CONFIGS.items():

    print(f"\nTraining {name}")

    roc_train_folds  = []
    roc_test_folds   = []
    roc_mass_folds   = {m: [] for m in mass_points}
    test_aucs_mass   = {m: [] for m in mass_points}
    train_aucs       = []
    test_aucs_global = []

    # Keep last fold for BDT distribution plot
    last_train_scores = last_y_train = None
    last_test_scores  = last_y_test  = None

    for train_idx, test_idx in kf.split(X, y):

        X_train, X_test = X.iloc[train_idx], X.iloc[test_idx]
        y_train, y_test = y.iloc[train_idx], y.iloc[test_idx]

        model = XGBClassifier(**{**BASE_PARAMS, **params})
        model.fit(X_train, y_train)

        # TRAIN (global)
        train_scores = model.predict_proba(X_train)[:, 1]
        train_aucs.append(roc_auc_score(y_train, train_scores))
        fpr_tr, tpr_tr, _ = roc_curve(y_train, train_scores)
        roc_train_folds.append((fpr_tr, tpr_tr))

        # TEST (global)
        scores = model.predict_proba(X_test)[:, 1]
        test_aucs_global.append(roc_auc_score(y_test, scores))
        fpr_te, tpr_te, _ = roc_curve(y_test, scores)
        roc_test_folds.append((fpr_te, tpr_te))

        last_train_scores = train_scores
        last_y_train      = y_train.values
        last_test_scores  = scores
        last_y_test       = y_test.values

        # PER-MASS
        for m in mass_points:

            sig = X_test[(y_test == 1) & (X_test["theta"] == m)]
            if len(sig) == 0:
                continue

            bkg    = X_test[y_test == 0]
            X_eval = pd.concat([sig, bkg]).copy()
            y_eval = np.concatenate([np.ones(len(sig)), np.zeros(len(bkg))])
            X_eval["theta"] = m

            scores_m = model.predict_proba(X_eval)[:, 1]
            auc = roc_auc_score(y_eval, scores_m)
            test_aucs_mass[m].append(auc)

            fpr_m, tpr_m, _ = roc_curve(y_eval, scores_m)
            roc_mass_folds[m].append((fpr_m, tpr_m))

            # Train mass AUC
            sig_train = X_train[(y_train == 1) & (X_train["theta"] == m)]
            if len(sig_train) > 0:
                bkg_train  = X_train[y_train == 0]
                X_tr_eval  = pd.concat([sig_train, bkg_train]).copy()
                y_tr_eval  = np.concatenate([np.ones(len(sig_train)), np.zeros(len(bkg_train))])
                X_tr_eval["theta"] = m
                train_auc_m = roc_auc_score(y_tr_eval, model.predict_proba(X_tr_eval)[:, 1])
            else:
                train_auc_m = np.nan

            mass_auc_records.append({
                "model": name, "mass": m,
                "train_auc": train_auc_m, "test_auc": auc,
            })

    # SUMMARY
    train_auc_mean = np.mean(train_aucs)
    test_auc_mean  = np.mean(test_aucs_global)
    delta_auc      = train_auc_mean - test_auc_mean

    mass_mean_auc = {m: np.mean(test_aucs_mass[m]) for m in mass_points if test_aucs_mass[m]}
    all_model_mass_auc[name] = mass_mean_auc

    #
    # PLOT 1 — All mass point ROC
    #
    fig, ax = plt.subplots()

    for fpr, tpr in roc_train_folds:
        ax.plot(fpr, tpr, color="tab:red",  alpha=0.55, lw=1.5)
    for fpr, tpr in roc_test_folds:
        ax.plot(fpr, tpr, color="tab:blue", alpha=0.55, lw=1.5)

    ax.plot([0, 1], [0, 1], "k--", lw=1)
    ax.plot([], [], color="tab:red",  lw=1.5, label="Train")
    ax.plot([], [], color="tab:blue", lw=1.5, label="Test")

    ax.text(
        0.55, 0.18,
        f"Train AUC = {train_auc_mean:.3f}\n"
        f"Test AUC  = {test_auc_mean:.3f}\n"
        f"\u0394AUC      = {delta_auc:.3f}",
        transform=ax.transAxes,
        bbox=dict(facecolor="white", edgecolor="black", alpha=0.85),
        fontsize=11, family="monospace",
    )

    ax.set_xlabel("False Positive Rate")
    ax.set_ylabel("True Positive Rate")
    ax.legend(loc="lower right")
    ax.set_title(f"Global ROC ({name})")
    fig.tight_layout()
    fig.savefig(f"{OUTDIR}/{name}_roc_global.png", dpi=300,bbox_inches="tight")
    plt.close(fig)

    
    # PLOT 2 — BDT Dist
    
    fig, ax = plt.subplots()

    bins = 50
    ax.hist(last_train_scores[last_y_train == 1], bins=bins, density=True,
            histtype="step", lw=1.5, color="tab:blue",   label="Train Signal")
    ax.hist(last_train_scores[last_y_train == 0], bins=bins, density=True,
            histtype="step", lw=1.5, color="tab:orange", linestyle="--", label="Train Bkg")
    ax.hist(last_test_scores[last_y_test == 1],   bins=bins, density=True,
            alpha=0.45, color="tab:green", label="Test Signal")
    ax.hist(last_test_scores[last_y_test == 0],   bins=bins, density=True,
            alpha=0.45, color="tab:red",   label="Test Bkg")

    ax.text(
        0.63, 0.72,
        f"Train AUC = {train_auc_mean:.3f}\n"
        f"Test AUC  = {test_auc_mean:.3f}",
        transform=ax.transAxes,
        bbox=dict(facecolor="white", edgecolor="black", alpha=0.85),
        fontsize=11, family="monospace",
    )

    ax.set_xlabel("BDT Score")
    ax.set_ylabel("Density")
    ax.legend(loc="upper right", fontsize=10)
    ax.set_title(f"BDT Output ({name})")
    fig.tight_layout()
    fig.savefig(f"{OUTDIR}/{name}_bdt.png", dpi=300,bbox_inches="tight")
    plt.close(fig)

    
    # PLOT 3 — MASS point ROC (fold-averaged, one curve per mass)
    
    cmap   = plt.get_cmap("tab10")
    colors = [cmap(i % 10) for i in range(len(mass_points))]

    fig, ax = plt.subplots()

    for i, m in enumerate(mass_points):
        if not roc_mass_folds[m]:
            continue
        mean_fpr, mean_tpr = avg_roc(roc_mass_folds[m])
        auc_m = mass_mean_auc.get(m, float("nan"))
        ax.plot(mean_fpr, mean_tpr, color=colors[i], lw=1.6,
                label=f"m={m} (AUC={auc_m:.3f})")

    ax.plot([0, 1], [0, 1], "k--", lw=1)
    ax.set_xlabel("False Positive Rate")
    ax.set_ylabel("True Positive Rate")
    ax.legend(loc="lower right", fontsize=9)
    ax.set_title(f"Mass ROC ({name})")
    fig.tight_layout()
    fig.savefig(f"{OUTDIR}/{name}_roc_mass.png", dpi=300,bbox_inches="tight")
    plt.close(fig)

    
    # PLOT 4 — AUC vs MASS
    
    masses = sorted(mass_mean_auc.keys())
    aucs   = [mass_mean_auc[m] for m in masses]

    fig, ax = plt.subplots()

    ax.plot(masses, aucs, "o-", color="tab:blue", lw=2, ms=6, label="Test AUC(m)")
    ax.axhline(test_auc_mean, color="tab:blue", lw=1.5, linestyle="--", label="Global Test AUC")

    auc_range = max(aucs) - min(aucs) if len(aucs) > 1 else 0.01
    for m, a in zip(masses, aucs):
        ax.annotate(
            f"{a:.3f}", xy=(m, a),
            xytext=(m, a + auc_range * 0.12),
            ha="center", fontsize=9,
        )

    ax.set_xlabel("Mass (GeV)")
    ax.set_ylabel("Test AUC")
    ax.legend(loc="upper right", fontsize=10)
    ax.set_title(f"AUC vs Mass ({name})")
    fig.tight_layout()
    fig.savefig(f"{OUTDIR}/{name}_auc_vs_mass.png", dpi=300,bbox_inches="tight")
    plt.close(fig)



df_mass_auc = pd.DataFrame(mass_auc_records)
df_mass_auc = df_mass_auc.groupby(["model", "mass"]).mean().reset_index()
df_mass_auc["delta_auc"] = df_mass_auc["train_auc"] - df_mass_auc["test_auc"]

df_mass_auc.to_csv(os.path.join(OUTDIR, "mass_auc_train_test.csv"), index=False)
print("Saved CSV with train/test AUC per mass.")
