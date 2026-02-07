import math
import glob
import ROOT
import pyLCIO
import os

exec(open("./plotHelper.py").read())
ROOT.gROOT.SetBatch()

PLOT_DIR = "/scratch/jwatts/mucol/mucolstudies/plots2026"
os.makedirs(PLOT_DIR, exist_ok=True)

samples = glob.glob("/scratch/jwatts/mucol/data/reco/pionGun_pT_0_50/pionGun_pT_0_50_reco_1.slcio")
files = {"pionGun_pT_0_50": samples}

label_map = {"pionGun_pT_0_50": "0-50 GeV"}

# HISTOGRAM SETUP (Only Eta)
hists = {}
for s in files:
    hists[s] = {}
    for obj in ["mcp_pi", "trk_pi_match", "pfo_pi_match", "clusters_pi_match"]:
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

            # Denominator
            mcp_pions = []
            for mcp in mcps:
                if mcp.getGeneratorStatus() != 1 or abs(mcp.getPDG()) != 211: continue
                tlv = getTLV(mcp)
                if tlv.E() < 10 or abs(tlv.Eta()) > 2.4: continue
                hists[s]["mcp_pi_eta"].Fill(tlv.Eta())
                mcp_pions.append(tlv)

            # Reconstructed Candidates with added Eta cuts
            pfo_pions = [getTLV(p) for p in pfos if getTLV(p).E() > 10 and abs(getTLV(p).Eta()) < 2.4 and abs(p.getType()) == 211]
            trk_pions = [getTrackTLV(t) for t in trks if getTrackTLV(t).E() > 10 and abs(getTrackTLV(t).Eta()) < 2.4]
            clu_pions = [getClusterTLV(c) for c in clusters if c.getEnergy() > 10 and abs(getClusterTLV(c).Eta()) < 2.4]

            # Matching 
            for mcp_pi in mcp_pions:
                # Tracks and Clusters still use Delta R matching
                if any(mcp_pi.DeltaR(t) < 0.2 for t in trk_pions):
                    hists[s]["trk_pi_match_eta"].Fill(mcp_pi.Eta())
                
                if any(mcp_pi.DeltaR(c) < 0.2 for c in clu_pions):
                    hists[s]["clusters_pi_match_eta"].Fill(mcp_pi.Eta())

                # PFO efficiency determined only by presence of PDGID 211 PFO in the event
                if len(pfo_pions) > 0:
                    hists[s]["pfo_pi_match_eta"].Fill(mcp_pi.Eta())

        reader.close()

print(f"\nSaving publication plots to {PLOT_DIR}...")
for s in hists:
    denom = hists[s]["mcp_pi_eta"]
    eff_map = {}
    
    display_names = {
        "trk_pi_match": "SiTracks",
        "pfo_pi_match": "PFOs with PDGID = 211",
        "clusters_pi_match": "PandoraClusters"
    }

    for obj in display_names:
        num = hists[s][obj + "_eta"]
        
        teff = ROOT.TEfficiency(num, denom)
        graph = teff.CreateGraph() # Necessary for plotHelper GetXaxis()
        eff_map[display_names[obj]] = graph
        
        # Save individual plot
        single_map = {display_names[obj]: graph}
        plotEfficiencies(single_map, os.path.join(PLOT_DIR, f"EFFICIENCY_{obj}.png"), xlabel="#eta", ylabel="Efficiency")

    plotEfficiencies(eff_map, os.path.join(PLOT_DIR, f"EFFICIENCY_COMBINED.png"), xlabel="#eta", ylabel="Efficiency")

for obj, label in [("trk_pi_match", "TRACK"), ("pfo_pi_match", "PFO"), ("clusters_pi_match", "CLUSTER")]:
    m = hists[s][obj + "_eta"].GetEntries()
    t = hists[s]["mcp_pi_eta"].GetEntries()
    print(f"{label}: matched = {int(m)}, total = {int(t)}, eff = {m/t if t>0 else 0:.3f}")
