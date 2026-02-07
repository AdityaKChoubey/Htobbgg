
import ROOT
import os

Fname ="/home/aditya/Desktop/ttbar_run3/ANALYSIS 1/TTtoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8_2024.root"
df = ROOT.RDataFrame("Events",Fname)

df0 = df

df0 = df0.Define(
    "nGenMuon",
    "Sum(ROOT::VecOps::abs(GenPart_pdgId) == 13)"
);


# Final-state GEN muons (status == 1) and counting  multiple copies of same muon as
#new column added which counts total final states muons in each event
df = df.Define(
    "nGenMuon_final",
    "Sum( (ROOT::VecOps::abs(GenPart_pdgId) == 13) && "
    "     (GenPart_status == 1) && "
    "     (GenPart_statusFlags & (1 << 13)) )"
);



totalGenMuons = df0.Sum("nGenMuon").GetValue()
print("Total gen muon entries =", totalGenMuons)
totalGenMuonsFinal = df.Sum("nGenMuon_final").GetValue()
print("Total final gen muons =", totalGenMuonsFinal)



#nBranches = len(df.GetColumnNames())
#print("Number of branches =", nBranches)
#nEvents = df.Count().GetValue()
#print("Number of events =", nEvents)


# Histogram of number of final-state GEN muons per event

h_nGenMuon_final = df.Histo1D(
    ("h_nGenMuon_final",
     "Final-state GEN muons per event;N_{#mu}^{GEN, final};Events",
     6, 0, 6),
    "nGenMuon_final"
)

# Number of events win each bin

#n_events_1mu = df.Filter("nGenMuon_final == 1").Count()
#n_events_2mu = df.Filter("nGenMuon_final == 2").Count()
#n_events_0mu = df.Filter("nGenMuon_final == 0").Count()
#n_events_3mu = df.Filter("nGenMuon_final == 3").Count()
#n_events_4mu = df.Filter("nGenMuon_final == 4").Count()
#n_events_5mu = df.Filter("nGenMuon_final == 5").Count()

#print("Number of events with exactly 0 final-state GEN muon:", n_events_0mu.GetValue())
#print("Number of events with exactly 1 final-state GEN muon:", n_events_1mu.GetValue())
#print("Number of events with exactly 2 final-state GEN muon:", n_events_2mu.GetValue())
#print("Number of events with exactly 3 final-state GEN muon:", n_events_3mu.GetValue())
#print("Number of events with exactly 4 final-state GEN muon:", n_events_4mu.GetValue())
#print("Number of events with exactly 5 final-state GEN muon:", n_events_5mu.GetValue())









#n_events_1mu = df.Filter("nMuon == 1").Count()
#n_events_2mu = df.Filter("nMuon == 2").Count()
#n_events_0mu = df.Filter("nMuon == 0").Count()
#n_events_3mu = df.Filter("nMuon == 3").Count()
#n_events_4mu = df.Filter("nMuon == 4").Count()
#n_events_5mu = df.Filter("nMuon == 5").Count()

#print("Number of events with exactly 0 final-state Reco muon:", n_events_0mu.GetValue())
#print("Number of events with exactly 1 final-state Reco muon:", n_events_1mu.GetValue())
#print("Number of events with exactly 2 final-state Reco muon:", n_events_2mu.GetValue())
#print("Number of events with exactly 3 final-state Reco muon:", n_events_3mu.GetValue())
#print("Number of events with exactly 4 final-state Reco muon:", n_events_4mu.GetValue())
#print("Number of events with exactly 5 final-state Reco muon:", n_events_5mu.GetValue())



#DEL R Matching ========================================

# Select events with ≥1 final-state GEN muon and >=1 Reco Muons
df_genFinal = df.Filter("nGenMuon_final >= 1 && nMuon >= 1")

#df_genFinal.Display(["nMuon"]).Print()

ROOT.gInterpreter.Declare("""
#include <cmath>
#include <ROOT/RVec.hxx>

using ROOT::VecOps::RVec;

float deltaR(float eta1, float phi1, float eta2, float phi2) {
    const float PI = 3.14159265358979323846f;
    float dphi = std::fabs(phi1 - phi2);
    if (dphi > PI) dphi = 2.f * PI - dphi;
    float deta = eta1 - eta2;
    return std::sqrt(deta*deta + dphi*dphi);
}

int matchLeadRecoToGen(
    float reco_eta,
    float reco_phi,
    const RVec<float>& gen_eta,
    const RVec<float>& gen_phi,
    const RVec<int>&   gen_pdgId,
    const RVec<int>&   gen_status,
    const RVec<int>&   gen_statusFlags
) {
    for (size_t i = 0; i < gen_eta.size(); ++i) {

        // Only physical final-state GEN muons
        if (std::abs(gen_pdgId[i]) != 13) continue;
        if (gen_status[i] != 1) continue;
        if (!(gen_statusFlags[i] & (1 << 13))) continue; // isLastCopy

        float dR = deltaR(
            reco_eta, reco_phi,
            gen_eta[i], gen_phi[i]
        );

        if (dR < 0.3f) {
            return i;
        }
    }
    return -1;
}
""")


# ======================================
# Leading RECO muon (pt-leading)
# ======================================
df_genFinal = df_genFinal.Define(
    "leadRecoIdx", "ROOT::VecOps::ArgMax(Muon_pt)"
)

df_genFinal = df_genFinal.Define(
    "leadReco_eta", "Muon_eta[leadRecoIdx]"
)

df_genFinal = df_genFinal.Define(
    "leadReco_phi", "Muon_phi[leadRecoIdx]"
)
df_genFinal = df_genFinal.Define(
    "leadReco_pt", "Muon_pt[leadRecoIdx]"
)

df_genFinal = df_genFinal.Define(
    "matchedGenIdx",
    "matchLeadRecoToGen("
    "leadReco_eta, leadReco_phi, "
    "GenPart_eta, GenPart_phi, "
    "GenPart_pdgId, GenPart_status, GenPart_statusFlags)"
)



# ======================================
# Matched GEN muon kinematics

df_genFinal = df_genFinal.Define(
    "matchedGen_pt",
    "matchedGenIdx >= 0 ? GenPart_pt[matchedGenIdx] : -1.f"
)

df_genFinal = df_genFinal.Define(
    "matchedGen_eta",
    "matchedGenIdx >= 0 ? GenPart_eta[matchedGenIdx] : -99.f"
)

df_genFinal = df_genFinal.Define(
    "matchedGen_phi",
    "matchedGenIdx >= 0 ? GenPart_phi[matchedGenIdx] : -99.f"
)


# ======================================
# Final-state GEN muons (before matching)

df_genFinal = df_genFinal.Define(
    "GenMuonFinal_pt",
    "GenPart_pt[(ROOT::VecOps::abs(GenPart_pdgId) == 13) && (GenPart_status == 1)]"
)

df_genFinal = df_genFinal.Define(
    "GenMuonFinal_eta",
    "GenPart_eta[(ROOT::VecOps::abs(GenPart_pdgId) == 13) && (GenPart_status == 1)]"
)

df_genFinal = df_genFinal.Define(
    "GenMuonFinal_phi",
    "GenPart_phi[(ROOT::VecOps::abs(GenPart_pdgId) == 13) && (GenPart_status == 1)]"
)

# ======================================
# Histograms: GEN muons before matching

h_gen_pt_all = df_genFinal.Histo1D(
    ("h_genFinal_pt", "Final-state GEN muon pT;pT [GeV];Entries", 100, 0, 400),
    "GenMuonFinal_pt"
)

h_gen_eta_all = df_genFinal.Histo1D(
    ("h_genFinal_eta", "Final-state GEN muon #eta;#eta;Entries", 80, -2.5, 2.5),
    "GenMuonFinal_eta"
)

h_gen_phi_all = df_genFinal.Histo1D(
    ("h_genFinal_phi", "Final-state GEN muon #phi;#phi;Entries", 64, -3.2, 3.2),
    "GenMuonFinal_phi"
)


# ======================================
# Histograms: matched GEN muon

h_matched_pt = df_genFinal.Filter("matchedGenIdx >= 0").Histo1D(
    ("h_matchedGen_pt", "Matched GEN muon pT;pT [GeV];Entries", 100, 0, 400),
    "matchedGen_pt"
)

h_matched_eta = df_genFinal.Filter("matchedGenIdx >= 0").Histo1D(
    ("h_matchedGen_eta", "Matched GEN muon #eta;#eta;Entries", 80, -2.5, 2.5),
    "matchedGen_eta"
)

h_matched_phi = df_genFinal.Filter("matchedGenIdx >= 0").Histo1D(
    ("h_matchedGen_phi", "Matched GEN muon #phi;#phi;Entries", 64, -3.2, 3.2),
    "matchedGen_phi"
)



# ======================================
# Write output ROOT file
# ======================================
out = ROOT.TFile("gen_reco_matching_histsv_2.root", "RECREATE")

out.mkdir("GEN_before_matching")
out.cd("GEN_before_matching")
h_gen_pt_all.GetValue().Write()
h_gen_eta_all.GetValue().Write()
h_gen_phi_all.GetValue().Write()

out.mkdir("GEN_matched")
out.cd("GEN_matched")
h_matched_pt.GetValue().Write()
h_matched_eta.GetValue().Write()
h_matched_phi.GetValue().Write()

out.Close()

df_genFinal.Snapshot("Events","Matched_muon.root")
# Retrieve histograms from RDataFrame
h_final   = h_gen_pt_all.GetValue()      # denominator
h_matched = h_matched_pt.GetValue()      # numerator

# --------------------------------------------------
# Constructing TEfficiency

h_eff = ROOT.TEfficiency(h_matched, h_final)

h_eff.SetTitle(
    "Gen and Reco muon matching efficiency;"
    "p_{T} [GeV];Efficiency"
)


# Default is Clopper-Pearson (exact)
# h_eff.SetStatisticOption(ROOT.TEfficiency.kClopperPearson)
# h_eff.SetStatisticOption(ROOT.TEfficiency.kWilson)   # alternative


out = ROOT.TFile("gen_reco_matching_histsv_2.root", "RECREATE")
h_eff.Write("h_eff_matched_gen")
out.Close()


ROOT.gInterpreter.Declare(r"""
#include <ROOT/RVec.hxx>
using ROOT::VecOps::RVec;

RVec<int> selectPromptWGenMuons(
    const RVec<int>& pdgId,
    const RVec<int>& status,
    const RVec<int>& statusFlags,
    const RVec<int>& motherIdx
) {
    RVec<int> mask(pdgId.size(), 0);

    for (size_t i = 0; i < pdgId.size(); ++i) {
        if (std::abs(pdgId[i]) != 13) continue;
        if (status[i] != 1) continue;
        if ((statusFlags[i] & (1 << 13)) == 0) continue;

        int midx = motherIdx[i];
        if (midx < 0) continue;
        if (std::abs(pdgId[midx]) != 24) continue;

        mask[i] = 1;
    }
    return mask;
}
""")


df = df.Define(
    "isPromptWGenMuon",
    "selectPromptWGenMuons("
    "GenPart_pdgId, "
    "GenPart_status, "
    "GenPart_statusFlags, "
    "GenPart_genPartIdxMother)"
)


df = df.Define(
    "nPromptWGenMuon",
    "Sum(isPromptWGenMuon)"
)
df = df.Define(
    "PromptWGenMuon_pt",
    "GenPart_pt[isPromptWGenMuon]"
)

df = df.Define(
    "PromptWGenMuon_eta",
    "GenPart_eta[isPromptWGenMuon]"
)

df = df.Define(
    "PromptWGenMuon_phi",
    "GenPart_phi[isPromptWGenMuon]"
)


# ======================================
# Histograms: GEN truth (prompt W muons)
# ======================================
h_truth_pt = df.Histo1D(
    ("h_truth_promptW_pt",
     "Prompt W-->mu GEN p_{T};p_{T} [GeV];Entries",
     100, 0, 400),
    "PromptWGenMuon_pt"
)

h_truth_eta = df.Histo1D(
    ("h_truth_promptW_eta",
     "Prompt W-->mu GEN #eta;#eta;Entries",
     150, -2.5, 2.5),
    "PromptWGenMuon_eta"
)

h_truth_phi = df.Histo1D(
    ("h_truth_promptW_phi",
     "Prompt W-->mu GEN #phi;#phi;Entries",
     64, -3.2, 3.2),
    "PromptWGenMuon_phi"
)



#First i matched the leading RECO muon to any GEN muon geometrically. Then i classify whether that matched GEN muon is a prompt W muon. Separately, we keep the full prompt W GEN muon truth for efficiency.


df_genFinal = df_genFinal.Define(
    "matchedIsPromptW",
    "matchedGenIdx >= 0 && "
    "std::abs(GenPart_pdgId[matchedGenIdx]) == 13 && "
    "GenPart_status[matchedGenIdx] == 1 && "
    "((GenPart_statusFlags[matchedGenIdx] & (1 << 13)) != 0) && "
    "GenPart_genPartIdxMother[matchedGenIdx] >= 0 && "
    "std::abs(GenPart_pdgId[GenPart_genPartIdxMother[matchedGenIdx]]) == 24"
)


h_matched_promptW_pt = df_genFinal.Filter("matchedIsPromptW").Histo1D(
    ("h_matchedPromptW_pt",
     "Matched prompt W--> mu GEN p_{T};p_{T} [GeV];Entries",
     100, 0, 400),
    "matchedGen_pt"
)

h_matched_promptW_eta = df_genFinal.Filter("matchedIsPromptW").Histo1D(
    ("h_matchedPromptW_eta",
     "Matched prompt W-->mu GEN #eta;#eta;Entries",
     150, -2.5, 2.5),
    "matchedGen_eta"
)

h_matched_promptW_phi = df_genFinal.Filter("matchedIsPromptW").Histo1D(
    ("h_matchedPromptW_phi",
     "Matched prompt W-->mu GEN #phi;#phi;Entries",
     64, -3.2, 3.2),
    "matchedGen_phi"
)


eff_pt = ROOT.TEfficiency(
    h_matched_promptW_pt.GetValue(),
    h_truth_pt.GetValue()
)
eff_pt.SetName("eff_promptW_pt")
eff_pt.SetTitle("Prompt (W --> mu) Gen mu vs Reco Matched Prompt Gen mu efficiency;p_{T} [GeV];Efficiency")

eff_eta = ROOT.TEfficiency(
    h_matched_promptW_eta.GetValue(),
    h_truth_eta.GetValue()
)
eff_eta.SetName("eff_promptW_eta")
eff_eta.SetTitle("Prompt (W --> mu) Gen mu vs Reco Matched Prompt Gen mu efficiency;#eta;Efficiency")

eff_phi = ROOT.TEfficiency(
    h_matched_promptW_phi.GetValue(),
    h_truth_phi.GetValue()
)
eff_phi.SetName("eff_promptW_phi")
eff_phi.SetTitle("Prompt W --> mu efficiency;#phi;Efficiency")


out = ROOT.TFile.Open("gen_reco_matching_histsv_2.root", "UPDATE")
if out.GetDirectory("PromptWMuon"):
    prompt_dir = out.GetDirectory("PromptWMuon")
else:
    prompt_dir = out.mkdir("PromptWMuon")

prompt_dir.cd()


h_truth_pt.GetValue().Write("h_truth_promptW_pt", ROOT.TObject.kOverwrite)
h_truth_eta.GetValue().Write("h_truth_promptW_eta", ROOT.TObject.kOverwrite)
h_truth_phi.GetValue().Write("h_truth_promptW_phi", ROOT.TObject.kOverwrite)
h_matched_promptW_pt.GetValue().Write("h_matchedPromptW_pt", ROOT.TObject.kOverwrite)
h_matched_promptW_eta.GetValue().Write("h_matchedPromptW_eta", ROOT.TObject.kOverwrite)
h_matched_promptW_phi.GetValue().Write("h_matchedPromptW_phi", ROOT.TObject.kOverwrite)
eff_pt.Write("eff_promptW_pt", ROOT.TObject.kOverwrite)
eff_eta.Write("eff_promptW_eta", ROOT.TObject.kOverwrite)
eff_phi.Write("eff_promptW_phi", ROOT.TObject.kOverwrite)
out.Close()
print("Prompt W muon histograms and TEfficiency objects stored successfully.")


#Transverse mass reconstruction


#Gen Prompt Mu and MET (Transvere component of Neutrino)

ROOT.gInterpreter.Declare(r"""
#include <ROOT/RVec.hxx>
using ROOT::VecOps::RVec;

RVec<int> selectPromptWGenNeutrinos(
    const RVec<int>& pdgId,
    const RVec<int>& status,
    const RVec<int>& statusFlags,
    const RVec<int>& motherIdx
) {
    RVec<int> mask(pdgId.size(), 0);

    for (size_t i = 0; i < pdgId.size(); ++i) {
        if (std::abs(pdgId[i]) != 14) continue;   // ν_μ
        if (status[i] != 1) continue;
        if ((statusFlags[i] & (1 << 13)) == 0) continue; // isLastCopy

        int midx = motherIdx[i];
        if (midx < 0) continue;
        if (std::abs(pdgId[midx]) != 24) continue; // W

        mask[i] = 1;
    }
    return mask;
}
""")


df_genFinal = df_genFinal.Define(
    "isPromptWGenNu",
    "selectPromptWGenNeutrinos("
    "GenPart_pdgId, "
    "GenPart_status, "
    "GenPart_statusFlags, "
    "GenPart_genPartIdxMother)"
)

df_genFinal = df_genFinal.Define(
    "PromptWGenNu_pt",
    "GenPart_pt[isPromptWGenNu][0]"
)

df_genFinal = df_genFinal.Define(
    "PromptWGenNu_phi",
    "GenPart_phi[isPromptWGenNu][0]"
)

ROOT.gInterpreter.Declare(r"""
float deltaPhi(float phi1, float phi2) {
    const float PI = 3.14159265358979323846f;
    float dphi = std::fabs(phi1 - phi2);
    if (dphi > PI) dphi = 2.f * PI - dphi;
    return dphi;
}
""")

df_genFinal = df_genFinal.Define(
    "MTW_genMu_genNu",
    "std::sqrt( 2.f * matchedGen_pt * PromptWGenNu_pt * "
    "(1.f - std::cos(deltaPhi(matchedGen_phi, PromptWGenNu_phi))) )"
)
h_mtw_gen = df_genFinal.Filter(
    "matchedIsPromptW"
).Histo1D(
    ("h_MTW_genMu_genNu",
     "M_{T}^{W} (GEN #mu matched to RECO + GEN #nu);M_{T} [GeV];Events",
     80, 0, 160),
    "MTW_genMu_genNu"
)
out = ROOT.TFile.Open("gen_reco_matching_histsv_2.root", "UPDATE")
out.cd()
h_mtw_gen.GetValue().Write("h_MTW_genMu_genNu", ROOT.TObject.kOverwrite)
out.Close()


df_genFinal = df_genFinal.Define(
    "isGenNu",
    "( (ROOT::VecOps::abs(GenPart_pdgId) == 12) || "
    "  (ROOT::VecOps::abs(GenPart_pdgId) == 14) || "
    "  (ROOT::VecOps::abs(GenPart_pdgId) == 16) ) && "
    "(GenPart_status == 1) && "
    "(GenPart_statusFlags & (1 << 13))"
)
df_genFinal = df_genFinal.Define(
    "GenMET_px",
    "Sum( GenPart_pt[isGenNu] * cos(GenPart_phi[isGenNu]) )"
)

df_genFinal = df_genFinal.Define(
    "GenMET_py",
    "Sum( GenPart_pt[isGenNu] * sin(GenPart_phi[isGenNu]) )"
)
df_genFinal = df_genFinal.Define(
    "genMET_pt_custom",
    "sqrt(GenMET_px*GenMET_px + GenMET_py*GenMET_py)"
)

df_genFinal = df_genFinal.Define(
    "genMET_phi_custom",
    "atan2(GenMET_py, GenMET_px)"
)


df_genFinal = df_genFinal.Define(
    "MTW_genMu_genMET_custom",
    "sqrt( 2.f * matchedGen_pt * genMET_pt_custom * "
    "(1.f - cos(deltaPhi(matchedGen_phi, genMET_phi_custom))) )"
)

h_mtw_genMET_custom = df_genFinal.Filter(
    "matchedIsPromptW"
).Histo1D(
    ("h_MTW_genMu_genMET_custom",
     "M_{T}^{W} (GEN #mu matched to RECO + custom GEN MET);"
     "M_{T} [GeV];Events",
     80, 0, 200),
    "MTW_genMu_genMET_custom"
)

out = ROOT.TFile.Open("gen_reco_matching_histsv_2.root", "UPDATE")
out.cd()
h_mtw_genMET_custom.GetValue().Write(
    "h_MTW_genMu_genMET_custom",
    ROOT.TObject.kOverwrite
)
out.Close()


df_genFinal = df_genFinal.Redefine(
    "leadReco_pt",
    "Muon_pt[leadRecoIdx]"
)

df_genFinal = df_genFinal.Redefine(
    "leadReco_phi",
    "Muon_phi[leadRecoIdx]"
)


df_genFinal = df_genFinal.Define(
    "MTW_recoMu_PuppiMET",
    "sqrt( 2.f * leadReco_pt * PuppiMET_pt * "
    "(1.f - cos(deltaPhi(leadReco_phi, PuppiMET_phi))) )"
)

h_mtw_reco_puppi = df_genFinal.Filter(
    "matchedIsPromptW"
).Histo1D(
    ("h_MTW_recoMu_PuppiMET",
     "M_{T}^{W} (Lead RECO #mu + PuppiMET);M_{T} [GeV];Events",
     80, 0, 200),
    "MTW_recoMu_PuppiMET"
)

out = ROOT.TFile.Open("gen_reco_matching_histsv_2.root", "UPDATE")
out.cd()
h_mtw_reco_puppi.GetValue().Write(
    "h_MTW_recoMu_PuppiMET",
    ROOT.TObject.kOverwrite
)
out.Close()





df_genAcc = df_genFinal.Filter(
    "matchedIsPromptW && "
    "matchedGen_pt > 26.0 && "
    "abs(matchedGen_eta) < 2.4"
)
df_recoAcc = df_genFinal.Filter(
    "matchedIsPromptW && "
    "leadReco_pt > 26.0 && "
    "abs(leadReco_eta) < 2.4"
)
h_mtw_gen_acc = df_genAcc.Histo1D(
    ("h_MTW_genMu_genNu_acc",
     "M_{T}^{W} (GEN #mu + GEN #nu, RECO acceptance);"
     "M_{T} [GeV];Events",
     80, 0, 160),
    "MTW_genMu_genNu"
)
h_mtw_genMET_acc = df_genAcc.Histo1D(
    ("h_MTW_genMu_genMET_acc",
     "M_{T}^{W} (GEN #mu + GEN MET, RECO acceptance);"
     "M_{T} [GeV];Events",
     80, 0, 200),
    "MTW_genMu_genMET_custom"
)
h_mtw_reco_acc = df_recoAcc.Histo1D(
    ("h_MTW_recoMu_PuppiMET_acc",
     "M_{T}^{W} (RECO #mu + PuppiMET, acceptance);"
     "M_{T} [GeV];Events",
     80, 0, 200),
    "MTW_recoMu_PuppiMET"
)
out = ROOT.TFile.Open("gen_reco_matching_histsv_2.root", "UPDATE")
out.cd()

h_mtw_gen_acc.GetValue().Write("h_MTW_genMu_genNu_acc", ROOT.TObject.kOverwrite)
h_mtw_genMET_acc.GetValue().Write("h_MTW_genMu_genMET_acc", ROOT.TObject.kOverwrite)
h_mtw_reco_acc.GetValue().Write("h_MTW_recoMu_PuppiMET_acc", ROOT.TObject.kOverwrite)

out.Close()



c = ROOT.TCanvas("c_MTW_acc", "W Transverse Mass (RECO-like acceptance)", 900, 700)
ROOT.gStyle.SetOptStat(0)

h1 = h_mtw_gen_acc.GetValue()
h2 = h_mtw_genMET_acc.GetValue()
h3 = h_mtw_reco_acc.GetValue()

# Normalized for shape comparison
for h in [h1, h2, h3]:
    if h.Integral() > 0:
        h.Scale(1.0 / h.Integral())
    h.SetLineWidth(2)

h1.SetLineColor(ROOT.kBlack)
h2.SetLineColor(ROOT.kBlue + 1)
h3.SetLineColor(ROOT.kRed + 1)

h1.SetTitle(
    "W Transverse Mass with RECO-like Acceptance;"
    "M_{T} [GeV];Normalized Events"
)

h1.Draw("HIST")
h2.Draw("HIST SAME")
h3.Draw("HIST SAME")

leg = ROOT.TLegend(0.55, 0.65, 0.88, 0.88)
leg.SetBorderSize(0)
leg.SetFillStyle(0)
leg.AddEntry(h1, "GEN #mu + GEN #nu", "l")
leg.AddEntry(h2, "GEN #mu + GEN MET", "l")
leg.AddEntry(h3, "RECO #mu + PuppiMET", "l")
leg.Draw()

c.SaveAs("MTW_acceptance_matched.png")
c.SaveAs("MTW_acceptance_matched.pdf")



# Muon Fake Rate (Leading RECO muon)
# Fake = RECO muon with NO matched GEN muon

# Denominator: all events with a leading RECO muon
h_reco_all_pt = df_genFinal.Histo1D(
    ("h_reco_all_pt",
     "All leading RECO muons;p_{T} [GeV];Entries",
     100, 0, 400),
    "leadReco_pt"
)

h_reco_all_eta = df_genFinal.Histo1D(
    ("h_reco_all_eta",
     "All leading RECO muons;#eta;Entries",
     80, -2.5, 2.5),
    "leadReco_eta"
)

# Numerator: fake RECO muons (no GEN match)
df_fake = df_genFinal.Filter("matchedGenIdx < 0")

h_fake_pt = df_fake.Histo1D(
    ("h_fake_pt",
     "Fake leading RECO muons;p_{T} [GeV];Entries",
     100, 0, 400),
    "leadReco_pt"
)

h_fake_eta = df_fake.Histo1D(
    ("h_fake_eta",
     "Fake leading RECO muons;#eta;Entries",
     80, -2.5, 2.5),
    "leadReco_eta"
)

# Fake rate using TEfficiency (binomial errors)
fakeRate_pt = ROOT.TEfficiency(
    h_fake_pt.GetValue(),
    h_reco_all_pt.GetValue()
)
fakeRate_pt.SetName("fakeRate_pt")
fakeRate_pt.SetTitle("Muon fake rate vs p_{T};p_{T} [GeV];Fake rate")

fakeRate_eta = ROOT.TEfficiency(
    h_fake_eta.GetValue(),
    h_reco_all_eta.GetValue()
)
fakeRate_eta.SetName("fakeRate_eta")
fakeRate_eta.SetTitle("Muon fake rate vs #eta;#eta;Fake rate")

# Store
out = ROOT.TFile.Open("gen_reco_matching_histsv_2.root", "UPDATE")
out.cd()
fakeRate_pt.Write("fakeRate_pt", ROOT.TObject.kOverwrite)
fakeRate_eta.Write("fakeRate_eta", ROOT.TObject.kOverwrite)
out.Close()

# Inclusive
n_all  = df_genFinal.Count().GetValue()
n_fake = df_genFinal.Filter("matchedGenIdx < 0").Count().GetValue()
print("Inclusive fake fraction =", n_fake / n_all)
print(n_fake)


















