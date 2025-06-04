import glob
import math
import ROOT
from pyLCIO import IOIMPL

# --- Matching functions ---
def match_deltaR(tlv1, tlv2):
    return tlv1.DeltaR(tlv2) < 0.1

def match_relPt(tlv1, tlv2):
    return abs(tlv1.Perp() - tlv2.Perp()) / tlv2.Perp() < 0.2

# --- File setup ---
files = [
    "electronGun_pT_0_50",
    "electronGun_pT_50_250",
    "electronGun_pT_250_1000",
    "electronGun_pT_1000_5000"
]

samples = glob.glob("/data/fmeloni/DataMuC_MuColl10_v0A/v0/reco/electronGun*")

label_map = {
    "electronGun_pT_0_50": "pT 0-50 GeV",
    "electronGun_pT_50_250": "pT 50-250 GeV",
    "electronGun_pT_250_1000": "pT 250-1000 GeV",
    "electronGun_pT_1000_5000": "pT 1000-5000 GeV"
}

# --- Histograms ---
hists = {}
for s in files:
    hists[s] = {
        "mcp_el_eta": ROOT.TH1F(f"{s}_mcp_el_eta", "", 50, -2.5, 2.5),
        "clusters_eta": ROOT.TH1F(f"{s}_clusters_eta", "", 50, -2.5, 2.5),
        "clusters_match_dR_eta": ROOT.TH1F(f"{s}_clusters_match_dR_eta", "", 50, -2.5, 2.5),
        "clusters_match_relPt_eta": ROOT.TH1F(f"{s}_clusters_match_relPt_eta", "", 50, -2.5, 2.5),
        "clusters_match_both_eta": ROOT.TH1F(f"{s}_clusters_match_both_eta", "", 50, -2.5, 2.5)
    }

# --- Processing loop ---
for s in files:
    for f in samples:
        if s in f:
            reader = IOIMPL.LCFactory.getInstance().createLCReader()
            reader.open(f)

            while True:
                try:
                    event = reader.readNextEvent()
                    if event is None:
                        break

                    mcps = event.getCollection("MCParticle")
                    clusters = event.getCollection("ReconClusters")

                    my_mcp_el = None
                    for mcp in mcps:
                        if not mcp.getGeneratorStatus() == 1:
                            continue
                        px, py, pz = mcp.getMomentum()
                        energy = mcp.getEnergy()
                        mcp_tlv = ROOT.TLorentzVector(px, py, pz, energy)
                        if abs(mcp.getPDG()) == 11 and abs(mcp_tlv.Eta()) <= 2:
                            hists[s]["mcp_el_eta"].Fill(mcp_tlv.Eta())
                            my_mcp_el = mcp_tlv

                    for cluster in clusters:
                        cluster_pos = cluster.getPosition()
                        cluster_E = cluster.getEnergy()
                        if cluster_E <= 0 or not my_mcp_el:
                            continue

                        cluster_vec = ROOT.TVector3(cluster_pos[0], cluster_pos[1], cluster_pos[2])
                        direction = cluster_vec.Unit()
                        cluster_tlv = ROOT.TLorentzVector()
                        cluster_tlv.SetVect(direction * cluster_E)
                        cluster_tlv.SetE(cluster_E)

                        hists[s]["clusters_eta"].Fill(cluster_tlv.Eta())

                        if match_deltaR(cluster_tlv, my_mcp_el):
                            hists[s]["clusters_match_dR_eta"].Fill(cluster_tlv.Eta())

                        if match_relPt(cluster_tlv, my_mcp_el):
                            hists[s]["clusters_match_relPt_eta"].Fill(cluster_tlv.Eta())

                        if match_deltaR(cluster_tlv, my_mcp_el) and match_relPt(cluster_tlv, my_mcp_el):
                            hists[s]["clusters_match_both_eta"].Fill(cluster_tlv.Eta())

                except Exception:
                    break

            reader.close()
            break

# --- Output summary ---
print("\nCluster Matching Efficiency Summary:")
for s in hists:
    total = hists[s]["mcp_el_eta"].GetEntries()
    dR = hists[s]["clusters_match_dR_eta"].GetEntries()
    relPt = hists[s]["clusters_match_relPt_eta"].GetEntries()
    both = hists[s]["clusters_match_both_eta"].GetEntries()
    print(f"{label_map[s]}: total MC electrons = {int(total)}, matched dR = {int(dR)}, relPt = {int(relPt)}, both = {int(both)}")

# --- Save plots ---
if ROOT.gSystem.AccessPathName("plots/test_match"):
    pass
else:
    ROOT.gSystem.mkdir("plots/test_match", True)

for s in hists:
    for key in hists[s]:
        c = ROOT.TCanvas("c", "c", 800, 600)
        hists[s][key].Draw()
        outname = f"plots/test_match/{s}_{key}.png"
        c.SaveAs(outname)

