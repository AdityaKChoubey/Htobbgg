import ROOT

#ROOT.EnableImplicitMT()

file = "/home/aditya/Downloads/TTtoLNu2Q_preselected.root"

df = ROOT.RDataFrame("Events", file)

df1 = (

    df

    # --------------------------------------------------
    # mother pdgId
    # --------------------------------------------------
    .Define(
        "mother_pdgId",
        "Take(GenPart_pdgId, GenPart_genPartIdxMother)"
    )

    # --------------------------------------------------
    # photon selection
    # --------------------------------------------------
    .Define(
        "ph_mask",
        "GenPart_pdgId == 22 && GenPart_status == 1 && GenPart_pt > 10 && "
        "(abs(GenPart_eta) < 1.4442 || "
        "(abs(GenPart_eta) > 1.566 && abs(GenPart_eta) < 2.5))"
    )

    # photon indices
    .Define(
        "ph_idx",
        "ROOT::VecOps::Nonzero(ph_mask)"
    )

    # photon kinematics
    .Define("ph_pt",  "Take(GenPart_pt, ph_idx)")
    .Define("ph_eta", "Take(GenPart_eta, ph_idx)")
    .Define("ph_phi", "Take(GenPart_phi, ph_idx)")

    .Filter("ph_pt.size() >= 2")

    # --------------------------------------------------
    # sort photons by pt
    # --------------------------------------------------
    .Define(
        "ph_order",
        "Reverse(Argsort(ph_pt))"
    )

    .Define(
        "ph_idx_sorted",
        "Take(ph_idx, ph_order)"
    )

    # leading photon
    .Define("pho1_idx", "ph_idx_sorted[0]")

    # subleading photon
    .Define("pho2_idx", "ph_idx_sorted[1]")

    # --------------------------------------------------
    # photon kinematics
    # --------------------------------------------------
    .Define("pho1_pt",  "GenPart_pt[pho1_idx]")
    .Define("pho1_eta", "GenPart_eta[pho1_idx]")
    .Define("pho1_phi", "GenPart_phi[pho1_idx]")

    .Define("pho2_pt",  "GenPart_pt[pho2_idx]")
    .Define("pho2_eta", "GenPart_eta[pho2_idx]")
    .Define("pho2_phi", "GenPart_phi[pho2_idx]")

    # --------------------------------------------------
    # photon mothers
    # --------------------------------------------------
    .Define(
        "pho1_mother_idx",
        "GenPart_genPartIdxMother[pho1_idx]"
    )

    .Define(
        "pho2_mother_idx",
        "GenPart_genPartIdxMother[pho2_idx]"
    )

    .Define(
        "pho1_mother_pdg",
        "pho1_mother_idx >= 0 ? GenPart_pdgId[pho1_mother_idx] : -999"
    )

    .Define(
        "pho2_mother_pdg",
        "pho2_mother_idx >= 0 ? GenPart_pdgId[pho2_mother_idx] : -999"
    )

    # --------------------------------------------------
    # photon grandmothers
    # --------------------------------------------------
    .Define(
        "pho1_grandmother_idx",
        "pho1_mother_idx >= 0 ? GenPart_genPartIdxMother[pho1_mother_idx] : -999"
    )

    .Define(
        "pho2_grandmother_idx",
        "pho2_mother_idx >= 0 ? GenPart_genPartIdxMother[pho2_mother_idx] : -999"
    )

    .Define(
        "pho1_grandmother_pdg",
        "pho1_grandmother_idx >= 0 ? GenPart_pdgId[pho1_grandmother_idx] : -999"
    )

    .Define(
        "pho2_grandmother_pdg",
        "pho2_grandmother_idx >= 0 ? GenPart_pdgId[pho2_grandmother_idx] : -999"
    )

    # --------------------------------------------------
    # photon separation
    # --------------------------------------------------
    .Define(
        "dphi",
        "ROOT::VecOps::DeltaPhi(pho1_phi,pho2_phi)"
    )

    .Define(
        "deta",
        "pho1_eta - pho2_eta"
    )

    .Define(
        "dr_gg",
        "sqrt(dphi*dphi + deta*deta)"
    )

)

# --------------------------------------------------
# Print events in zero ΔR bin
# --------------------------------------------------



print("\nEvents with very small photon separation (dr_gg < 0.05)\n")

df1.Filter("dr_gg < 0.05").Display([
"dr_gg",
"pho1_mother_pdg","pho2_mother_pdg",
"pho1_grandmother_pdg","pho2_grandmother_pdg"
],20).Print()

# --------------------------------------------------
# Count how many such events exist
# --------------------------------------------------

print("\nNumber of dr_gg < 0.05 events:")
print(df1.Filter("dr_gg < 0.05").Count().GetValue())
