import ROOT
import sys

ROOT.ROOT.EnableImplicitMT(8)

input_file  = sys.argv[1]
output_file = sys.argv[2]

tree_name = "Events"

ROOT.gInterpreter.Declare(r"""
#include <vector>
#include <cmath>
#include <algorithm>
#include "ROOT/RVec.hxx"
#include "TLorentzVector.h"

using namespace ROOT::VecOps;

// ---------- BASIC ----------
float deltaPhi(float phi1, float phi2){
    float dphi = phi1 - phi2;
    while (dphi > M_PI) dphi -= 2*M_PI;
    while (dphi <= -M_PI) dphi += 2*M_PI;
    return std::abs(dphi);
}

float deltaR(float eta1, float phi1, float eta2, float phi2){
    float deta = eta1 - eta2;
    float dphi = deltaPhi(phi1, phi2);
    return std::sqrt(deta*deta + dphi*dphi);
}

// ---------- ELECTRONS ----------
RVec<int> SelectElectrons(const RVec<float>& pt,
                          const RVec<float>& eta,
                          const RVec<float>& iso,
                          const RVec<bool>& mva80)
{
    RVec<int> idx;
    for (size_t i = 0; i < pt.size(); ++i){
        float aeta = std::abs(eta[i]);
        bool pass_eta = (aeta < 2.5) && !(aeta > 1.44 && aeta < 1.57);
        if (pt[i] > 30 && pass_eta && iso[i] < 0.15 && mva80[i])
            idx.push_back(i);
    }
    return idx;
}

// ---------- MUONS ----------
RVec<int> SelectMuons(const RVec<float>& pt,
                      const RVec<float>& eta,
                      const RVec<float>& iso)
{
    RVec<int> idx;
    for (size_t i = 0; i < pt.size(); ++i){
        if (pt[i] > 24 && std::abs(eta[i]) < 2.4 && iso[i] < 0.15)
            idx.push_back(i);
    }
    return idx;
}

// ---------- PHOTONS ----------
RVec<int> SelectPhotons(const RVec<float>& pt,
                        const RVec<float>& eta,
                        const RVec<float>& mvaID,
                        const RVec<bool>& pixelSeed)
{
    RVec<int> idx;
    for (size_t i = 0; i < pt.size(); ++i){
        float aeta = std::abs(eta[i]);

        bool valid_eta = (aeta <= 2.5) && !(aeta >= 1.442 && aeta <= 1.566);

        bool barrel = aeta < 1.442;
        bool endcap = (aeta > 1.566 && aeta < 2.5);

        bool mva_pass = (barrel && mvaID[i] > 0.0439603) ||
                        (endcap && mvaID[i] > -0.249526);

        if (pt[i] > 10 && valid_eta && mva_pass && !pixelSeed[i])
            idx.push_back(i);
    }
    return idx;
}

// ---------- JETS ----------
RVec<int> SelectJets(const RVec<float>& pt,
                     const RVec<float>& eta,
                     const RVec<float>& btag)
{
    RVec<int> idx;
    for (size_t i = 0; i < pt.size(); ++i){
        if (pt[i] > 20 && std::abs(eta[i]) < 2.4 && btag[i] > 0.0246)
            idx.push_back(i);
    }
    return idx;
}

// ---------- DR CLEANING ----------
bool PassDR(const RVec<int>& lepIdx,
            const RVec<float>& lep_eta,
            const RVec<float>& lep_phi,
            const RVec<int>& phoIdx,
            const RVec<float>& pho_eta,
            const RVec<float>& pho_phi)
{
    for (auto l : lepIdx){
        for (auto p : phoIdx){
            if (deltaR(lep_eta[l], lep_phi[l], pho_eta[p], pho_phi[p]) < 0.4)
                return false;
        }
    }
    return true;
}

// ---------- SORT ----------
RVec<int> SortByPt(const RVec<int>& idx, const RVec<float>& pt)
{
    RVec<int> out = idx;
    std::sort(out.begin(), out.end(),
              [&](int a, int b){ return pt[a] > pt[b]; });
    return out;
}

// ---------- 4-VECTOR ----------
TLorentzVector MakeP4(float pt, float eta, float phi, float m){
    TLorentzVector v;
    v.SetPtEtaPhiM(pt, eta, phi, m);
    return v;
}

float DiMass(float pt1, float eta1, float phi1, float m1,
             float pt2, float eta2, float phi2, float m2)
{
    return (MakeP4(pt1,eta1,phi1,m1) + MakeP4(pt2,eta2,phi2,m2)).M();
}

float DiPhi(float pt1, float eta1, float phi1, float m1,
            float pt2, float eta2, float phi2, float m2)
{
    return (MakeP4(pt1,eta1,phi1,m1) + MakeP4(pt2,eta2,phi2,m2)).Phi();
}
""")

df = ROOT.RDataFrame(tree_name, input_file)

# ---------- SELECTION ----------
df = (
    df
    .Define("selEleIdx", "SelectElectrons(Electron_pt, Electron_eta, Electron_pfRelIso03_all, Electron_mvaIso_WP80)")
    .Define("selMuIdx",  "SelectMuons(Muon_pt, Muon_eta, Muon_pfRelIso03_all)")
    .Define("selPhoIdx", "SelectPhotons(Photon_pt, Photon_eta, Photon_mvaID, Photon_pixelSeed)")
    .Define("selJetIdx", "SelectJets(Jet_pt, Jet_eta, Jet_btagUParTAK4B)")
)

df = (
    df
    .Define("oneEle", "selEleIdx.size() == 1")
    .Define("oneMu",  "selMuIdx.size() == 1")
    .Define("lepMask", "oneEle || oneMu")

    .Filter("lepMask")
    .Filter("selPhoIdx.size() >= 2")
    .Filter("selJetIdx.size() >= 2")

    .Filter("PassDR(selEleIdx, Electron_eta, Electron_phi, selPhoIdx, Photon_eta, Photon_phi)")
    .Filter("PassDR(selMuIdx, Muon_eta, Muon_phi, selPhoIdx, Photon_eta, Photon_phi)")
)

# ---------- SORT ----------
df = (
    df
    .Define("pho_sorted", "SortByPt(selPhoIdx, Photon_pt)")
    .Define("jet_sorted", "SortByPt(selJetIdx, Jet_pt)")
)

# ---------- LEADING ----------
df = (
    df
    .Define("g1", "pho_sorted[0]")
    .Define("g2", "pho_sorted[1]")
    .Define("b1", "jet_sorted[0]")
    .Define("b2", "jet_sorted[1]")
)

# ---------- LEPTON ----------
df = (
    df
    .Define("lep_pt",
        "selEleIdx.size()==1 ? Electron_pt[selEleIdx[0]] : Muon_pt[selMuIdx[0]]")
)

# ---------- FEATURES ----------
df = (
    df
    .Define("g1_pt", "Photon_pt[g1]")
    .Define("g2_pt", "Photon_pt[g2]")
    .Define("b1_pt", "Jet_pt[b1]")
    .Define("b2_pt", "Jet_pt[b2]")

    .Define("b1_eta", "Jet_eta[b1]")
    .Define("b2_eta", "Jet_eta[b2]")

    .Define("mgg", "DiMass(g1_pt, Photon_eta[g1], Photon_phi[g1], 0.0, g2_pt, Photon_eta[g2], Photon_phi[g2], 0.0)")
    .Define("mbb", "DiMass(b1_pt, Jet_eta[b1], Jet_phi[b1], Jet_mass[b1], b2_pt, Jet_eta[b2], Jet_phi[b2], Jet_mass[b2])")

    .Define("dphi_gg", "deltaPhi(Photon_phi[g1], Photon_phi[g2])")
    .Define("dphi_bb", "deltaPhi(Jet_phi[b1], Jet_phi[b2])")

    .Define("gg_phi", "DiPhi(g1_pt, Photon_eta[g1], Photon_phi[g1], 0.0, g2_pt, Photon_eta[g2], Photon_phi[g2], 0.0)")
    .Define("bb_phi", "DiPhi(b1_pt, Jet_eta[b1], Jet_phi[b1], Jet_mass[b1], b2_pt, Jet_eta[b2], Jet_phi[b2], Jet_mass[b2])")

    .Define("dphi_bbgg", "deltaPhi(bb_phi, gg_phi)")
    .Define("dphi_gg_met", "deltaPhi(gg_phi, PuppiMET_phi)")

    .Define("label", "0")
    .Define("mass_point", "-1.0")
)

columns = [
    "g1_pt", "g2_pt",
    "b1_pt", "b2_pt",
    "b1_eta", "b2_eta",
    "mgg", "mbb",
    "dphi_gg", "dphi_bb",
    "dphi_bbgg", "dphi_gg_met",
    "lep_pt",
    "label", "mass_point"
]

df.Snapshot("Events", output_file, columns)
print("Saved:", output_file)