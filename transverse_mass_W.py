import ROOT

# Enable multithreading for faster RDataFrame execution
ROOT.EnableImplicitMT()

# ============================================================
# Load snapshot produced after GEN–RECO matching
# This file already contains:
#  - matchedGenIdx
#  - matched GEN muon kinematics
#  - leading RECO muon kinematics
# ============================================================
df = ROOT.RDataFrame("Events", "Matched_muon.root")

# Keep only events where the leading RECO muon
# was successfully matched to a GEN muon
df_matched = df.Filter("matchedGenIdx >= 0")

# ============================================================
# Basic sanity plots: GEN vs RECO muon pT
# ============================================================

# Leading RECO muon pT (GEN-matched)
h_leadReco_pt = df_matched.Histo1D(
    ("h_leadReco_pt",
     "Leading RECO muon p_{T} (GEN-matched);p_{T}^{RECO} [GeV];Events",
     100, 0, 400),
    "leadReco_pt"
)

# Matched GEN muon pT
h_matchedGen_pt = df_matched.Histo1D(
    ("h_matchedGen_pt",
     "Matched GEN muon p_{T};p_{T}^{GEN} [GeV];Events",
     100, 0, 400),
    "matchedGen_pt"
)

# 2D correlation: RECO vs GEN muon pT
h2_pt = df_matched.Histo2D(
    ("h2_reco_vs_gen_pt",
     "Leading RECO p_{T} vs Matched GEN p_{T};p_{T}^{GEN} [GeV];p_{T}^{RECO} [GeV]",
     100, 0, 400,
     100, 0, 400),
    "matchedGen_pt",
    "leadReco_pt"
)

# Write these sanity-check histograms
out = ROOT.TFile("output.root", "RECREATE")
h_leadReco_pt.GetValue().Write()
h_matchedGen_pt.GetValue().Write()
h2_pt.GetValue().Write()
out.Close()

# ============================================================
# Define transverse mass function (shared by GEN and RECO)
# ============================================================

ROOT.gInterpreter.Declare("""
#include <cmath>

float mtw(float lep_pt, float lep_phi,
          float met_pt, float met_phi)
{
    float dphi = std::fabs(lep_phi - met_phi);
    if (dphi > M_PI) dphi = 2.f * M_PI - dphi;

    return std::sqrt(
        2.f * lep_pt * met_pt * (1.f - std::cos(dphi))
    );
}
""")

# ============================================================
# Compute MTW at RECO and GEN level (no cuts)
# ============================================================

df_matched = df_matched.Define(
    "mtw_reco",
    "mtw(leadReco_pt, leadReco_phi, PuppiMET_pt, PuppiMET_phi)"
)

df_matched = df_matched.Define(
    "mtw_gen",
    "mtw(matchedGen_pt, matchedGen_phi, GenMET_pt, GenMET_phi)"
)

# MTW distributions
h_mtw_reco = df_matched.Histo1D(
    ("h_mtw_reco",
     "W transverse mass (RECO);M_{T}^{W} [GeV];Events",
     120, 0, 240),
    "mtw_reco"
)

h_mtw_gen = df_matched.Histo1D(
    ("h_mtw_gen",
     "W transverse mass (GEN);M_{T}^{W} [GeV];Events",
     120, 0, 240),
    "mtw_gen"
)

# GEN vs RECO MTW correlation
h2_mtw = df_matched.Histo2D(
    ("h2_mtw_reco_vs_gen",
     "RECO vs GEN W transverse mass;M_{T}^{GEN} [GeV];M_{T}^{RECO} [GeV]",
     120, 0, 240,
     120, 0, 240),
    "mtw_gen",
    "mtw_reco"
)

out = ROOT.TFile("output.root", "UPDATE")
h_mtw_reco.GetValue().Write()
h_mtw_gen.GetValue().Write()
h2_mtw.GetValue().Write()
out.Close()

# ============================================================
# Apply analysis-level kinematic cuts
# ============================================================

df_matched_cut = (
    df_matched
    .Filter("leadReco_pt > 25")
    .Filter("abs(leadReco_eta) < 2.4")
    .Filter("matchedGen_pt > 25")
    .Filter("abs(matchedGen_eta) < 2.4")
    .Filter("PuppiMET_pt > 25")
    .Filter("GenMET_pt > 25")
)

# Muon pT after cuts
h_leadReco_pt_cut = df_matched_cut.Histo1D(
    ("h_leadReco_pt_cut",
     "Leading RECO muon p_{T} (after cuts);p_{T}^{RECO} [GeV];Events",
     100, 0, 400),
    "leadReco_pt"
)

h_matchedGen_pt_cut = df_matched_cut.Histo1D(
    ("h_matchedGen_pt",
     "Matched GEN muon p_{T} (after cuts);p_{T}^{GEN} [GeV];Events",
     100, 0, 400),
    "matchedGen_pt"
)

# MTW after cuts
df_matched_cut = df_matched_cut.Define(
    "mtw_reco_cut",
    "mtw(leadReco_pt, leadReco_phi, PuppiMET_pt, PuppiMET_phi)"
)

df_matched_cut = df_matched_cut.Define(
    "mtw_gen_cut",
    "mtw(matchedGen_pt, matchedGen_phi, GenMET_pt, GenMET_phi)"
)

h_mtw_reco_cut = df_matched_cut.Histo1D(
    ("h_mtw_reco_cut",
     "W transverse mass (RECO, cuts);M_{T}^{W} [GeV];Events",
     180, 0, 240),
    "mtw_reco_cut"
)

h_mtw_gen_cut = df_matched_cut.Histo1D(
    ("h_mtw_gen_cut",
     "W transverse mass (GEN, cuts);M_{T}^{W} [GeV];Events",
     180, 0, 240),
    "mtw_gen_cut"
)

# GEN vs RECO MTW after cuts
h2_mtw_cut = df_matched_cut.Histo2D(
    ("h2_mtw_reco_vs_gen_cut",
     "RECO vs GEN W transverse mass (cuts);"
     "M_{T}^{GEN} [GeV];M_{T}^{RECO} [GeV]",
     180, 0, 240,
     180, 0, 240),
    "mtw_gen_cut",
    "mtw_reco_cut"
)

out = ROOT.TFile("mtw_after_cuts.root", "RECREATE")
h_mtw_reco_cut.GetValue().Write()
h_mtw_gen_cut.GetValue().Write()
h2_mtw_cut.GetValue().Write()
out.Close()

# ============================================================
# MET validation: GEN MET vs PuppiMET
# ============================================================

h_gen_met_pt = df_matched_cut.Histo1D(
    ("h_gen_met_pt",
     "GEN MET p_{T};p_{T}^{GEN MET} [GeV];Events",
     100, 0, 400),
    "GenMET_pt"
)

h_puppi_met_pt = df_matched_cut.Histo1D(
    ("h_puppi_met_pt",
     "Puppi MET p_{T};p_{T}^{Puppi MET} [GeV];Events",
     100, 0, 400),
    "PuppiMET_pt"
)

h2_gen_vs_puppi_met = df_matched_cut.Histo2D(
    ("h2_gen_vs_puppi_met",
     "GEN MET vs Puppi MET (after cuts);"
     "p_{T}^{GEN MET} [GeV];p_{T}^{Puppi MET} [GeV]",
     100, 0, 400,
     100, 0, 400),
    "GenMET_pt",
    "PuppiMET_pt"
)

# ============================================================
# Compute W transverse momentum (boost) at GEN and RECO level
# ============================================================

# GEN W pT from mu + nu
df_matched_cut = df_matched_cut.Define(
    "W_gen_px",
    "matchedGen_pt * cos(matchedGen_phi) + GenMET_pt * cos(GenMET_phi)"
)

df_matched_cut = df_matched_cut.Define(
    "W_gen_py",
    "matchedGen_pt * sin(matchedGen_phi) + GenMET_pt * sin(GenMET_phi)"
)

df_matched_cut = df_matched_cut.Define(
    "boostW_gen",
    "sqrt(W_gen_px*W_gen_px + W_gen_py*W_gen_py)"
)

# RECO W pT from mu + MET
df_matched_cut = df_matched_cut.Define(
    "W_reco_px",
    "leadReco_pt * cos(leadReco_phi) + PuppiMET_pt * cos(PuppiMET_phi)"
)

df_matched_cut = df_matched_cut.Define(
    "W_reco_py",
    "leadReco_pt * sin(leadReco_phi) + PuppiMET_pt * sin(PuppiMET_phi)"
)

df_matched_cut = df_matched_cut.Define(
    "boostW_reco",
    "sqrt(W_reco_px*W_reco_px + W_reco_py*W_reco_py)"
)

# ============================================================
# MTW in W-boost bins (GEN and RECO definitions)
# ============================================================

boost_bins = [(0,20),(20,40),(40,60),(60,100),(100,200)]

h_mtw_gen_in_genBoost   = []
h_mtw_reco_in_genBoost  = []
h_mtw_gen_in_recoBoost  = []
h_mtw_reco_in_recoBoost = []

# Bin by GEN W boost
for lo, hi in boost_bins:
    df_bin = df_matched_cut.Filter(f"boostW_gen >= {lo} && boostW_gen < {hi}")

    h_mtw_gen_in_genBoost.append(
        df_bin.Histo1D(
            (f"h_mtw_gen_genBoost_{lo}_{hi}",
             f"GEN M_T^W, GEN boost {lo}-{hi} GeV;M_T [GeV];Events",
             120, 0, 240),
            "mtw_gen_cut"
        )
    )

    h_mtw_reco_in_genBoost.append(
        df_bin.Histo1D(
            (f"h_mtw_reco_genBoost_{lo}_{hi}",
             f"RECO M_T^W, GEN boost {lo}-{hi} GeV;M_T [GeV];Events",
             120, 0, 240),
            "mtw_reco_cut"
        )
    )

# Bin by RECO W boost
for lo, hi in boost_bins:
    df_bin = df_matched_cut.Filter(f"boostW_reco >= {lo} && boostW_reco < {hi}")

    h_mtw_gen_in_recoBoost.append(
        df_bin.Histo1D(
            (f"h_mtw_gen_recoBoost_{lo}_{hi}",
             f"GEN M_T^W, RECO boost {lo}-{hi} GeV;M_T [GeV];Events",
             120, 0, 240),
            "mtw_gen_cut"
        )
    )

    h_mtw_reco_in_recoBoost.append(
        df_bin.Histo1D(
            (f"h_mtw_reco_recoBoost_{lo}_{hi}",
             f"RECO M_T^W, RECO boost {lo}-{hi} GeV;M_T [GeV];Events",
             120, 0, 240),
            "mtw_reco_cut"
        )
    )

# Write boost-binned MTW histograms
out = ROOT.TFile("mtw_in_boost_bins.root", "RECREATE")
for hlist in [h_mtw_gen_in_genBoost,
              h_mtw_reco_in_genBoost,
              h_mtw_gen_in_recoBoost,
              h_mtw_reco_in_recoBoost]:
    for h in hlist:
        h.GetValue().Write()
out.Close()
