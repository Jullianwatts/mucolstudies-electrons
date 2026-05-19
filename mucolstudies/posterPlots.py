import math
import glob
import ROOT
import pyLCIO
import os

exec(open("./plotHelper.py").read())
ROOT.gROOT.SetBatch()

PLOT_DIR = "/scratch/jwatts/mucol/mucolstudies/plots2026"
os.makedirs(PLOT_DIR, exist_ok=True)

# Define two samples (using different files to represent the different Pandora versions)
files = {
    "sample_1": glob.glob("/scratch/jwatts/mucol/data/reco/electronGun_pT_0_50/electronGun_pT_0_50_reco_9.slcio"),
    "sample_2": glob.glob("/scratch/jwatts/mucol/data/reco/electronGun_pT_0_50/electronGun_pT_0_50_reco_7.slcio") 
}

label_map = {
    "sample_1": "no changes to Pandora",
    "sample_2": "updated electron id in Pandora"
}

hists = {}
for s in files:
    hists[s] = {}
    for obj in ["mcp_el", "trk_el_match", "pfo_el_match", "clusters_el_match"]:
        hists[s][obj + "_eta"] = ROOT.TH1F(f"{s}_{obj}_eta", ";#eta;Counts", 50, -2.5, 2.5)

reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()

for s in files:
    for f in files[s]:
        reader.open(f)
        for event in reader:
            mcps = event.getCollection("MCParticle")
            pfos = event.getCollection("PandoraPFOs")
            trks = event.getCollection("SiTracks")
            clusters = event.getCollection("PandoraClusters")

            # Truth selection (1 per event)
            mcp_el = None
            for mcp in mcps:
                # Updated truth control logic
                if abs(mcp.getPDG()) == 11 and mcp.getGeneratorStatus() == 1 and mcp.getEnergy() > 10:
                    tlv = getTLV(mcp)
                    if abs(tlv.Eta()) > 2.4: continue
                    mcp_el = tlv
                    hists[s]["mcp_el_eta"].Fill(tlv.Eta())
                    break 

            if not mcp_el: continue

            # Reconstruction object selection (2 GeV cut)
            trk_tlvs = [getTrackTLV(t) for t in trks if getTrackTLV(t).E() > 2]
            clu_tlvs = [getClusterTLV(c) for c in clusters if c.getEnergy() > 2]
            pfo_electrons = [getTLV(p) for p in pfos if abs(p.getType()) == 11 and getTLV(p).E() > 2]

            # Track Matching
            if any(mcp_el.DeltaR(t) < 0.2 for t in trk_tlvs):
                hists[s]["trk_el_match_eta"].Fill(mcp_el.Eta())
            
            # Cluster Matching
            if any(mcp_el.DeltaR(c) < 0.2 for c in clu_tlvs):
                hists[s]["clusters_el_match_eta"].Fill(mcp_el.Eta())

            # PFO Matching (ID as electron + DeltaR < 0.2)
            if any(mcp_el.DeltaR(p) < 0.2 for p in pfo_electrons):
                hists[s]["pfo_el_match_eta"].Fill(mcp_el.Eta())

        reader.close()

# Plot 1: Track and Cluster Efficiency on one plot (using Sample 1)
s1 = "sample_1"
eff_track_cluster = {
    "SiTracks": ROOT.TEfficiency(hists[s1]["trk_el_match_eta"], hists[s1]["mcp_el_eta"]).CreateGraph(),
    "PandoraClusters": ROOT.TEfficiency(hists[s1]["clusters_el_match_eta"], hists[s1]["mcp_el_eta"]).CreateGraph()
}
plotEfficiencies(eff_track_cluster, os.path.join(PLOT_DIR, "EFFICIENCY_TRACK_AND_CLUSTER.svg"), xlabel="#eta", ylabel="Efficiency")

# Plot 2: PFO Efficiency (ID as electron) comparison between Pandora versions
eff_pfo_comparison = {}
for s in files:
    eff_pfo_comparison[label_map[s]] = ROOT.TEfficiency(hists[s]["pfo_el_match_eta"], hists[s]["mcp_el_eta"]).CreateGraph()

plotEfficiencies(eff_pfo_comparison, os.path.join(PLOT_DIR, "EFFICIENCY_PFO_ELECTRON_ID.svg"), xlabel="#eta", ylabel="Efficiency")

for s in files:
    print(f"\nStats for {label_map[s]}:")
    for obj in ["trk_el_match", "clusters_el_match", "pfo_el_match"]:
        m = hists[s][obj + "_eta"].GetEntries()
        t = hists[s]["mcp_el_eta"].GetEntries()
        print(f"{obj}: matched = {int(m)}, total = {int(t)}, eff = {m/t if t>0 else 0:.3f}")
