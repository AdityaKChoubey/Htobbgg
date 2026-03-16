import uproot
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_curve, roc_auc_score
import xgboost as xgb


# -------------------------------------------------
# Load ROOT file
# -------------------------------------------------

file = "training_dataset.root"
tree = "Events"

features = [
"dr_bb","dr_gg",
"dr_b1g1","dr_b1g2",
"dr_b2g1","dr_b2g2",
"deta_bb","deta_gg",
"dphi_bb","dphi_gg"
]

branches = features + ["label"]

df = uproot.open(file)[tree].arrays(branches, library="pd")

print("Dataset size:", len(df))


# -------------------------------------------------
# Train/Test split
# -------------------------------------------------

X = df[features]
y = df["label"]

X_train, X_test, y_train, y_test = train_test_split(
X, y,
test_size=0.3,
random_state=42,
stratify=y
)


# -------------------------------------------------
# Hyperparameters
# -------------------------------------------------

param_sets = [

{"max_depth":3,"learning_rate":0.1,"n_estimators":200},
{"max_depth":5,"learning_rate":0.1,"n_estimators":300},
{"max_depth":6,"learning_rate":0.05,"n_estimators":400},
{"max_depth":8,"learning_rate":0.03,"n_estimators":500}

]


# -------------------------------------------------
# ROC figure
# -------------------------------------------------

fig_roc, ax_roc = plt.subplots(figsize=(8,6))


# -------------------------------------------------
# Train models
# -------------------------------------------------

for i,params in enumerate(param_sets):

    print("\n==============================")
    print("Model",i+1)
    print("Params:",params)

    model = xgb.XGBClassifier(
        objective="binary:logistic",
        eval_metric="logloss",
        use_label_encoder=False,
        n_jobs = 8,
        tree_method = "hist",
        **params
    )

    model.fit(X_train,y_train)

    train_prob = model.predict_proba(X_train)[:,1]
    test_prob = model.predict_proba(X_test)[:,1]

    train_auc = roc_auc_score(y_train,train_prob)
    test_auc = roc_auc_score(y_test,test_prob)

    print("Train AUC:",train_auc)
    print("Test AUC:",test_auc)


    # -------------------------
    # ROC curve
    # -------------------------

    fpr,tpr,_ = roc_curve(y_test,test_prob)

    ax_roc.plot(
        fpr,
        tpr,
        label=f"Model {i+1} (AUC={test_auc:.3f})"
    )


    # -------------------------------------------------
    # BDT Score distribution
    # -------------------------------------------------

    fig, ax = plt.subplots()

    ax.hist(train_prob[y_train==1],bins=40,density=True,
            histtype="step",label="Signal (train)")

    ax.hist(train_prob[y_train==0],bins=40,density=True,
            histtype="step",label="Background (train)")

    ax.hist(test_prob[y_test==1],bins=40,density=True,
            alpha=0.4,label="Signal (test)")

    ax.hist(test_prob[y_test==0],bins=40,density=True,
            alpha=0.4,label="Background (test)")

    ax.set_xlabel("BDT score")
    ax.set_ylabel("Normalized events")
    ax.set_title(f"BDT Score Distribution (Model {i+1})")

    ax.legend()
    ax.grid()

    fig.savefig(f"bdt_score_model{i+1}.png")
    plt.close(fig)


    # -------------------------------------------------
    # Feature importance
    # -------------------------------------------------

    fig, ax = plt.subplots()

    xgb.plot_importance(
        model,
        importance_type="gain",
        ax=ax
    )

    ax.set_title(f"Feature Importance (Model {i+1})")

    fig.savefig(f"feature_importance_model{i+1}.png")
    plt.close(fig)



# -------------------------------------------------
# Final ROC plot
# -------------------------------------------------

ax_roc.plot([0,1],[0,1],'k--')

ax_roc.set_xlabel("False Positive Rate")
ax_roc.set_ylabel("True Positive Rate")
ax_roc.set_title("ROC Curve Comparison")

ax_roc.legend()
ax_roc.grid()

fig_roc.savefig("roc_models.png")
plt.show()
