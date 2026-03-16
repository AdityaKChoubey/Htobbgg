import ROOT

ROOT.EnableImplicitMT()

file = "root://eosuser.cern.ch//eos/user/a/achoubey/htobbgg/TTtoLNu2Q_preselected.root"

df = ROOT.RDataFrame("Events", file)

df1 = (
    df

    # mother pdgId
    .Define("mother_pdgId",
        "Take(GenPart_pdgId, GenPart_genPartIdxMother)")

    # -------------------------
    # b and bbar from top decay
    # -------------------------
    .Define("b_mask",
        "GenPart_pdgId == 5 && abs(mother_pdgId) == 6")

    .Define("bbar_mask",
        "GenPart_pdgId == -5 && abs(mother_pdgId) == 6")

    .Define("b_pt",  "GenPart_pt[b_mask]")
    .Define("b_eta", "GenPart_eta[b_mask]")
    .Define("b_phi", "GenPart_phi[b_mask]")

    .Define("bbar_pt",  "GenPart_pt[bbar_mask]")
    .Define("bbar_eta", "GenPart_eta[bbar_mask]")
    .Define("bbar_phi", "GenPart_phi[bbar_mask]")

    # leading b and bbar
    #.Define("b_pt_lead",  "b_pt.size()>0 ? b_pt[0] : -1")
    #.Define("b_eta_lead", "b_eta.size()>0 ? b_eta[0] : 0")
    #.Define("b_phi_lead", "b_phi.size()>0 ? b_phi[0] : 0")

    #.Define("bbar_pt_lead",  "bbar_pt.size()>0 ? bbar_pt[0] : -1")
    #.Define("bbar_eta_lead", "bbar_eta.size()>0 ? bbar_eta[0] : 0")
    #.Define("bbar_phi_lead", "bbar_phi.size()>0 ? bbar_phi[0] : 0")

    # -------------------------
    # photons
    # -------------------------
    .Define("ph_mask",
        "GenPart_pdgId == 22 && GenPart_status == 1 && GenPart_pt > 10 && "
        "(abs(GenPart_eta) < 1.4442 || "
        "(abs(GenPart_eta) > 1.566 && abs(GenPart_eta) < 2.5))")

    .Define("ph_pt",  "GenPart_pt[ph_mask]")
    .Define("ph_eta", "GenPart_eta[ph_mask]")
    .Define("ph_phi", "GenPart_phi[ph_mask]")

    .Filter("ph_pt.size() >= 2")

    # sort photons by pt
    .Define("ph_order",
        "Reverse(Argsort(ph_pt))")

    .Define("ph_pt_sorted",
        "Take(ph_pt, ph_order)")
    .Define("ph_eta_sorted",
        "Take(ph_eta, ph_order)")
    .Define("ph_phi_sorted",
        "Take(ph_phi, ph_order)")

    # leading photon
    .Define("pho1_pt",  "ph_pt_sorted[0]")
    .Define("pho1_eta", "ph_eta_sorted[0]")
    .Define("pho1_phi", "ph_phi_sorted[0]")

    # subleading photon
    .Define("pho2_pt",  "ph_pt_sorted[1]")
    .Define("pho2_eta", "ph_eta_sorted[1]")
    .Define("pho2_phi", "ph_phi_sorted[1]")
)

branches = [

"b_pt","b_eta","b_phi",
"bbar_pt","bbar_eta","bbar_phi",

"pho1_pt","pho1_eta","pho1_phi",
"pho2_pt","pho2_eta","pho2_phi"

]

df1.Snapshot(
    "Events",
    "root://eosuser.cern.ch//eos/user/a/achoubey/htobbgg/ttbar_training_gen.root",
    branches
)

