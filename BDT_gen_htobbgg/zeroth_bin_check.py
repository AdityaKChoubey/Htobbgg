import ROOT

# ROOT.EnableImplicitMT()

file = "/home/aditya/Downloads/TTtoLNu2Q_preselected.root"

df = ROOT.RDataFrame("Events", file)

df1 = (

    df

    # --------------------------------------------------
    # Safe mother pdgId
    # --------------------------------------------------
    .Define(
        "mother_pdgId",
        "Where(GenPart_genPartIdxMother >= 0, "
        "Take(GenPart_pdgId, GenPart_genPartIdxMother), -999)"
    )

    # --------------------------------------------------
    # b and bbar from top decay
    # --------------------------------------------------
    .Define(
        "b_mask",
        "GenPart_pdgId == 5 && abs(mother_pdgId) == 6"
    )

    .Define(
        "bbar_mask",
        "GenPart_pdgId == -5 && abs(mother_pdgId) == 6"
    )

    .Define("b_pt","GenPart_pt[b_mask]")
    .Define("b_eta","GenPart_eta[b_mask]")
    .Define("b_phi","GenPart_phi[b_mask]")

    .Define("bbar_pt","GenPart_pt[bbar_mask]")
    .Define("bbar_eta","GenPart_eta[bbar_mask]")
    .Define("bbar_phi","GenPart_phi[bbar_mask]")

    # --------------------------------------------------
    # combine b and bbar
    # --------------------------------------------------
    .Define("b_all_pt","Concatenate(b_pt,bbar_pt)")
    .Define("b_all_eta","Concatenate(b_eta,bbar_eta)")
    .Define("b_all_phi","Concatenate(b_phi,bbar_phi)")

    .Filter("b_all_pt.size() >= 2")

    # --------------------------------------------------
    # sort b quarks by pt
    # --------------------------------------------------
    .Define("b_order","Reverse(Argsort(b_all_pt))")

    .Define("b_pt_sorted","Take(b_all_pt,b_order)")
    .Define("b_eta_sorted","Take(b_all_eta,b_order)")
    .Define("b_phi_sorted","Take(b_all_phi,b_order)")

    .Define("b1_pt","b_pt_sorted[0]")
    .Define("b1_eta","b_eta_sorted[0]")
    .Define("b1_phi","b_phi_sorted[0]")

    .Define("b2_pt","b_pt_sorted[1]")
    .Define("b2_eta","b_eta_sorted[1]")
    .Define("b2_phi","b_phi_sorted[1]")

    # --------------------------------------------------
    # photon selection
    # --------------------------------------------------
    .Define(
        "ph_mask",
        "GenPart_pdgId == 22 && GenPart_status == 1 && GenPart_pt > 10 && "
        "(abs(GenPart_eta) < 1.4442 || "
        "(abs(GenPart_eta) > 1.566 && abs(GenPart_eta) < 2.5))"
    )

    .Define("ph_idx","ROOT::VecOps::Nonzero(ph_mask)")

    .Define("ph_pt","Take(GenPart_pt, ph_idx)")
    .Define("ph_eta","Take(GenPart_eta, ph_idx)")
    .Define("ph_phi","Take(GenPart_phi, ph_idx)")

    .Filter("ph_pt.size() >= 2")

    # --------------------------------------------------
    # sort photons by pt
    # --------------------------------------------------
    .Define("ph_order","Reverse(Argsort(ph_pt))")

    .Define("ph_idx_sorted","Take(ph_idx, ph_order)")

    .Define("pho1_idx","ph_idx_sorted[0]")
    .Define("pho2_idx","ph_idx_sorted[1]")

    # --------------------------------------------------
    # photon kinematics
    # --------------------------------------------------
    .Define("pho1_pt","GenPart_pt[pho1_idx]")
    .Define("pho1_eta","GenPart_eta[pho1_idx]")
    .Define("pho1_phi","GenPart_phi[pho1_idx]")

    .Define("pho2_pt","GenPart_pt[pho2_idx]")
    .Define("pho2_eta","GenPart_eta[pho2_idx]")
    .Define("pho2_phi","GenPart_phi[pho2_idx]")

    # --------------------------------------------------
    # photon mothers
    # --------------------------------------------------
    .Define("pho1_mother_idx","GenPart_genPartIdxMother[pho1_idx]")
    .Define("pho2_mother_idx","GenPart_genPartIdxMother[pho2_idx]")

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
    # photon-photon separation
    # --------------------------------------------------
    .Define("dphi_gg","ROOT::VecOps::DeltaPhi(pho1_phi,pho2_phi)")
    .Define("deta_gg","pho1_eta - pho2_eta")
    .Define("dr_gg","sqrt(dphi_gg*dphi_gg + deta_gg*deta_gg)")

    # --------------------------------------------------
    # b-photon separations
    # --------------------------------------------------
    .Define("dphi_b1g1","ROOT::VecOps::DeltaPhi(b1_phi,pho1_phi)")
    .Define("dphi_b1g2","ROOT::VecOps::DeltaPhi(b1_phi,pho2_phi)")
    .Define("dphi_b2g1","ROOT::VecOps::DeltaPhi(b2_phi,pho1_phi)")
    .Define("dphi_b2g2","ROOT::VecOps::DeltaPhi(b2_phi,pho2_phi)")

    .Define("deta_b1g1","b1_eta - pho1_eta")
    .Define("deta_b1g2","b1_eta - pho2_eta")
    .Define("deta_b2g1","b2_eta - pho1_eta")
    .Define("deta_b2g2","b2_eta - pho2_eta")

    .Define("dr_b1g1","sqrt(dphi_b1g1*dphi_b1g1 + deta_b1g1*deta_b1g1)")
    .Define("dr_b1g2","sqrt(dphi_b1g2*dphi_b1g2 + deta_b1g2*deta_b1g2)")
    .Define("dr_b2g1","sqrt(dphi_b2g1*dphi_b2g1 + deta_b2g1*deta_b2g1)")
    .Define("dr_b2g2","sqrt(dphi_b2g2*dphi_b2g2 + deta_b2g2*deta_b2g2)")
)

# --------------------------------------------------
# Print gamma-gamma near zero
# --------------------------------------------------

print("\nEvents with dr_gg < 0.05\n")

df1.Filter("dr_gg < 0.05").Display([
"dr_gg",
"pho1_mother_pdg","pho2_mother_pdg",
"pho1_grandmother_pdg","pho2_grandmother_pdg"
],20).Print()

# --------------------------------------------------
# Check photon overlapping with b
# --------------------------------------------------

print("\nEvents with photon overlapping b (dr_b1g2 < 0.02)\n")

df1.Filter("dr_b1g2 < 0.02").Display([
"dr_b1g2",
"pho2_mother_pdg",
"pho2_grandmother_pdg"
],20).Print()


print("\nEvents with photon1 overlapping b1 (dr_b1g1 < 0.02)\n")

df1.Filter("dr_b1g1 < 0.02").Display([
"dr_b1g1",
"pho1_mother_pdg",
"pho1_grandmother_pdg"
],20).Print()


print("\nEvents with photon2 overlapping b1 (dr_b1g2 < 0.02)\n")

df1.Filter("dr_b1g2 < 0.02").Display([
"dr_b1g2",
"pho2_mother_pdg",
"pho2_grandmother_pdg"
],20).Print()


print("\nEvents with photon1 overlapping b2 (dr_b2g1 < 0.02)\n")

df1.Filter("dr_b2g1 < 0.02").Display([
"dr_b2g1",
"pho1_mother_pdg",
"pho1_grandmother_pdg"
],20).Print()

print("\nEvents with photon2 overlapping b2 (dr_b2g2 < 0.02)\n")

df1.Filter("dr_b2g2 < 0.02").Display([
"dr_b2g2",
"pho2_mother_pdg",
"pho2_grandmother_pdg"
],20).Print()
print("\nCounts of overlaps\n")

print("b1-g1:", df1.Filter("dr_b1g1 < 0.02").Count().GetValue())
print("b1-g2:", df1.Filter("dr_b1g2 < 0.02").Count().GetValue())
print("b2-g1:", df1.Filter("dr_b2g1 < 0.02").Count().GetValue())
print("b2-g2:", df1.Filter("dr_b2g2 < 0.02").Count().GetValue())
