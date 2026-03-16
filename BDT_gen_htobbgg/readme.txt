README



Workflow



1. Generate background features from ttbar GEN-level information
2. Inspect kinematic distributions and investigate unusual peaks
3. Prepare the training dataset
4. Train multiple BDT models with different hyperparameters
5. Compare model performance using ROC and AUC


Files Description


1. bkg_gen.py

extracts required gen info using GEN-level information from the ttbar sample.

Tasks performed:
- Select b and bbar quarks from top decay
- Save pT, eta and phi into a ROOT file

Variables include:
dr_bb
dr_gg
dr_b1g1
dr_b1g2
dr_b2g1
dr_b2g2
deta_bb
deta_gg
dphi_bb
dphi_gg


2. input_feature_hcobbgg.py
Variables include:
dr_bb
dr_gg
dr_b1g1
dr_b1g2
dr_b2g1
dr_b2g2
deta_bb
deta_gg
dphi_bb
dphi_gg

Generates the input features used for BDT training.

This script:
- Reads signal and bkg_gen samples
- pt orders signal photons and b quarks
- pt orders bkg_gen b quarks
- computes input variables and saves in root file



3. dr_gg_check.py

Investigates the peak near zero in the ΔR_gg distribution.

This script:
- Finds events with very small photon separation
- Checks the mother and grandmother particles of photons
- Identifies the origin of the photons

Observation:
Many photons originate from neutral meson decays such as
pi0 -> gamma gamma
eta -> gamma gamma
inside jets.


4. zeroth_bin_check.py

Further investigates events in the first bin of dR_gg ,dR_g1b2 , dR_g2,b1 ,etc.

Tasks performed:
- Inspect photon ancestry
- Check overlap between photons and b quarks
- checks the origin of  photons

Conclusion:
The peak near zero arises from neutral meson decays inside b-jets.


5. xgboost_bdt_v2.py

Main script for training the BDT classifier using XGBoost.

The script:
- Loads the training dataset
- Splits data into train and test samples
- Trains four BDT models with different hyperparameters
- Compares model performance

Metrics computed:
- ROC curve
- AUC
- BDT score distributions
- Feature importance


Model Configurations Tested

Four BDT models were trained with different hyperparameters:

Model 1: max_depth = 3, learning_rate = 0.1, n_estimators = 200
Model 2: max_depth = 5, learning_rate = 0.1, n_estimators = 300
Model 3: max_depth = 6, learning_rate = 0.05, n_estimators = 400
Model 4: max_depth = 8, learning_rate = 0.03, n_estimators = 500

The models were compared using ROC curves and AUC values.



Author

Aditya Kumar Choubey

