import uproot
import awkward as ak
import numpy as np

# ----------------------------
# Utility functions
# ----------------------------

def dphi(phi1, phi2):
    d = phi1 - phi2
    d = np.where(d > np.pi, d - 2*np.pi, d)
    d = np.where(d < -np.pi, d + 2*np.pi, d)
    return d

def dr(eta1,phi1,eta2,phi2):
    dp = dphi(phi1,phi2)
    de = eta1-eta2
    return np.sqrt(dp**2 + de**2)

# ----------------------------
# Background
# ----------------------------

bkg_file = "root://eosuser.cern.ch//eos/user/a/achoubey/htobbgg/ttbar_training_gen.root"

bkg = uproot.open(bkg_file)["Events"].arrays(library="ak")


#pt order b
b_eta  = ak.to_numpy(ak.firsts(bkg["b_eta"]))
b_phi  = ak.to_numpy(ak.firsts(bkg["b_phi"]))
b_pt   = ak.to_numpy(ak.firsts(bkg["b_pt"]))

bbar_eta = ak.to_numpy(ak.firsts(bkg["bbar_eta"]))
bbar_phi = ak.to_numpy(ak.firsts(bkg["bbar_phi"]))
bbar_pt  = ak.to_numpy(ak.firsts(bkg["bbar_pt"]))

eta = np.stack([b_eta, bbar_eta], axis=1)
phi = np.stack([b_phi, bbar_phi], axis=1)
pt  = np.stack([b_pt,  bbar_pt],  axis=1)

order = np.argsort(-pt, axis=1)
eta_sorted = np.take_along_axis(eta, order, axis=1)
phi_sorted = np.take_along_axis(phi, order, axis=1)
pt_sorted  = np.take_along_axis(pt,  order, axis=1)

b1_eta = eta_sorted[:,0] #lead
b1_phi = phi_sorted[:,0]
b1_pt  = pt_sorted[:,0]

b2_eta = eta_sorted[:,1]  #sublead
b2_phi = phi_sorted[:,1]
b2_pt  = pt_sorted[:,1]


g1_eta = bkg["pho1_eta"]  #already pt ordered pho1 is lead photon
g1_phi = bkg["pho1_phi"]  #already pt ordered ph02 is sublead

g2_eta = bkg["pho2_eta"]
g2_phi = bkg["pho2_phi"]

bkg_dict = {

"dr_bb":dr(b1_eta,b1_phi,b2_eta,b2_phi),
"dr_gg":dr(g1_eta,g1_phi,g2_eta,g2_phi),

"dr_b1g1":dr(b1_eta,b1_phi,g1_eta,g1_phi),
"dr_b1g2":dr(b1_eta,b1_phi,g2_eta,g2_phi),

"dr_b2g1":dr(b2_eta,b2_phi,g1_eta,g1_phi),
"dr_b2g2":dr(b2_eta,b2_phi,g2_eta,g2_phi),

"deta_bb":b1_eta-b2_eta,
"deta_gg":g1_eta-g2_eta,

"dphi_bb":dphi(b1_phi,b2_phi),
"dphi_gg":dphi(g1_phi,g2_phi),

"label":np.zeros(len(b1_eta))
}

n_bkg = len(b1_eta)

print("Background events:",n_bkg)

# ----------------------------
# Signal
# ----------------------------

sig_file =  "root://eosuser.cern.ch//eos/user/a/achoubey/htobbgg/WH_M40_RunIISummer20UL18NanoAODv9_genlevel.root"

sig = uproot.open(sig_file)["Events"].arrays(library="np")


# ----------------------------
# Photon cuts
# ----------------------------

g1_pt  = sig["pho1_pt"]
g2_pt  = sig["pho2_pt"]

g1_eta = sig["pho1_eta"]
g2_eta = sig["pho2_eta"]

def photon_eta_cut(eta):
    a = np.abs(eta)
    return (a < 1.442) | ((a > 1.556) & (a < 2.5))

mask = (
    (g1_pt > 10) &
    (g2_pt > 10) &
    photon_eta_cut(g1_eta) &
    photon_eta_cut(g2_eta)
)

sig = {k: v[mask] for k,v in sig.items()}

print("Signal after photon cuts:", len(sig["pho1_eta"]))

# ==================================================
# Photon pt ordering
# ==================================================

pho_eta = np.stack([sig["pho1_eta"], sig["pho2_eta"]], axis=1)
pho_phi = np.stack([sig["pho1_phi"], sig["pho2_phi"]], axis=1)
pho_pt  = np.stack([sig["pho1_pt"],  sig["pho2_pt"]],  axis=1)

order = np.argsort(-pho_pt, axis=1)

pho_eta = np.take_along_axis(pho_eta, order, axis=1)
pho_phi = np.take_along_axis(pho_phi, order, axis=1)
pho_pt  = np.take_along_axis(pho_pt,  order, axis=1)

g1_eta = pho_eta[:,0]
g1_phi = pho_phi[:,0]

g2_eta = pho_eta[:,1]
g2_phi = pho_phi[:,1]

# ==================================================
# b-quark pt ordering
# ==================================================

b_eta = np.stack([sig["b1_eta"], sig["b2_eta"]], axis=1)
b_phi = np.stack([sig["b1_phi"], sig["b2_phi"]], axis=1)
b_pt  = np.stack([sig["b1_pt"],  sig["b2_pt"]],  axis=1)

order = np.argsort(-b_pt, axis=1)

b_eta = np.take_along_axis(b_eta, order, axis=1)
b_phi = np.take_along_axis(b_phi, order, axis=1)
b_pt  = np.take_along_axis(b_pt,  order, axis=1)

b1_eta = b_eta[:,0]
b1_phi = b_phi[:,0]

b2_eta = b_eta[:,1]
b2_phi = b_phi[:,1]
sig_dict = {

"dr_bb":dr(b1_eta,b1_phi,b2_eta,b2_phi),
"dr_gg":dr(g1_eta,g1_phi,g2_eta,g2_phi),

"dr_b1g1":dr(b1_eta,b1_phi,g1_eta,g1_phi),
"dr_b1g2":dr(b1_eta,b1_phi,g2_eta,g2_phi),

"dr_b2g1":dr(b2_eta,b2_phi,g1_eta,g1_phi),
"dr_b2g2":dr(b2_eta,b2_phi,g2_eta,g2_phi),

"deta_bb":b1_eta-b2_eta,
"deta_gg":g1_eta-g2_eta,

"dphi_bb":dphi(b1_phi,b2_phi),
"dphi_gg":dphi(g1_phi,g2_phi),

"label":np.ones(len(b1_eta))
}
# Randomly pick signal


rng = np.random.default_rng()

n = min(len(sig_dict["label"]), n_bkg)

idx = rng.choice(len(sig_dict["label"]), size=n, replace=False)

sig_small = {k:v[idx] for k,v in sig_dict.items()}


# Merge datasets


merged = {}

for k in bkg_dict:
    merged[k] = np.concatenate([sig_small[k], bkg_dict[k]])


# Shuffle dataset

perm = rng.permutation(len(merged["label"]))

for k in merged:
    merged[k] = merged[k][perm]

# ----------------------------
# Save to ROOT
# ----------------------------

with uproot.recreate("training_dataset.root") as f:
    f["Events"] = merged

print("Saved training_dataset.root")