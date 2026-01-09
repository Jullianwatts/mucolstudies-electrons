import math
import glob
import ROOT
import pyLCIO
import os

exec(open("./plotHelper.py").read())
ROOT.gROOT.SetBatch()

# =========================
# CONFIGURATION
# =========================
PLOT_DIR = "/scratch/jwatts/mucol/mucolstudies/plots2026"
os.makedirs(PLOT_DIR, exist_ok=True)

samples = glob.glob("/scratch/jwatts/mucol/data/reco/electronGun_pT_0_50/electronGun_pT_0_50_reco_2.slcio")
files = {"electronGun_pT_0_50": samples}

label_map = {"electronGun_pT_0_50": "0-50 GeV"}

# =========================
# MINIMAL HISTOGRAM SETUP (Only Eta)
# =========================
hists = {}
for s in files:
    hists[s] = {}
    # ONLY keeping Eta histograms as requested
    for obj in ["mcp_el", "trk_el_match", "pfo_el_match", "clusters_el_match"]:
        hists[s][obj + "_eta"] = ROOT.TH1F(f"{s}_{obj}_eta", ";#eta;Counts", 50, -2.5, 2.5)

reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()

# =========================
# EVENT LOOP (Your trusted logic)
# =========================
for s in files:
    for f in files[s]:
        reader.open(f)
        for event in reader:
            mcps = event.getCollection("MCParticle")
            pfos = event.getCollection("PandoraPFOs")
            trks = event.getCollection("SiTracks")
            clusters = event.getCollection("PandoraClusters")

            # Denominator
            mcp_electrons = []
            for mcp in mcps:
                if mcp.getGeneratorStatus() != 1 or abs(mcp.getPDG()) != 11: continue
                tlv = getTLV(mcp)
                if tlv.E() < 20 or abs(tlv.Eta()) > 2.4: continue
                hists[s]["mcp_el_eta"].Fill(tlv.Eta())
                mcp_electrons.append(tlv)

            # Reconstructed Candidates
            pfo_electrons = [getTLV(p) for p in pfos if getTLV(p).E() > 10 and abs(p.getType()) == 11]
            trk_electrons = [getTrackTLV(t) for t in trks if getTrackTLV(t).E() > 10]
            clu_electrons = [getClusterTLV(c) for c in clusters if c.getEnergy() > 10]

            # Matching (dR < 0.1)
            for mcp_el in mcp_electrons:
                if any(mcp_el.DeltaR(t) < 0.1 for t in trk_electrons):
                    hists[s]["trk_el_match_eta"].Fill(mcp_el.Eta())
                if any(mcp_el.DeltaR(p) < 0.1 for p in pfo_electrons):
                    hists[s]["pfo_el_match_eta"].Fill(mcp_el.Eta())
                if any(mcp_el.DeltaR(c) < 0.1 for c in clu_electrons):
                    hists[s]["clusters_el_match_eta"].Fill(mcp_el.Eta())
        reader.close()

# =========================
# SAVE PUBLICATION PNGs
# =========================
print(f"\nSaving publication plots to {PLOT_DIR}...")
for s in hists:
    denom = hists[s]["mcp_el_eta"]
    eff_map = {}
    
    # Map for clean legend labels
    display_names = {
        "trk_el_match": "SiTracks",
        "pfo_el_match": "PandoraPFOs",
        "clusters_el_match": "PandoraClusters"
    }

    for obj in display_names:
        num = hists[s][obj + "_eta"]
        
        # Correct statistical handling
        teff = ROOT.TEfficiency(num, denom)
        graph = teff.CreateGraph() # Necessary for plotHelper GetXaxis()
        eff_map[display_names[obj]] = graph
        
        # Save individual plot
        single_map = {display_names[obj]: graph}
        plotEfficiencies(single_map, os.path.join(PLOT_DIR, f"EFFICIENCY_{obj}.png"), xlabel="#eta", ylabel="Efficiency")

    # SAVE COMBINED PLOT WITH LEGEND
    plotEfficiencies(eff_map, os.path.join(PLOT_DIR, f"EFFICIENCY_COMBINED.png"), xlabel="#eta", ylabel="Efficiency")

# Final Summaries
for obj, label in [("trk_el_match", "TRACK"), ("pfo_el_match", "PFO"), ("clusters_el_match", "CLUSTER")]:
    m = hists[s][obj + "_eta"].GetEntries()
    t = hists[s]["mcp_el_eta"].GetEntries()
    print(f"{label}: matched = {int(m)}, total = {int(t)}, eff = {m/t if t>0 else 0:.3f}")
