import math
import glob
import ROOT
import pyLCIO
import os

exec(open("./plotHelper.py").read())
ROOT.gROOT.SetBatch()

PLOT_DIR = "/scratch/jwatts/mucol/mucolstudies/plots2026"
os.makedirs(PLOT_DIR, exist_ok=True)

samples = glob.glob("/scratch/jwatts/mucol/data/reco/electronGun_pT_0_50/electronGun_pT_0_50_reco_7.slcio")
files = {"electronGun_pT_0_50": samples}

hists = {}
for s in files:
    hists[s] = {}
    # pfos and clusters commented out per request
    for obj in ["mcp_el", "trk_el_match"]: 
        hists[s][obj + "_central"] = ROOT.TH1F(f"{s}_{obj}_central", ";#theta [rad];Counts", 50, 0.698, 2.443)
        hists[s][obj + "_forward"] = ROOT.TH1F(f"{s}_{obj}_forward", ";#theta [rad];Counts", 50, 0.0, 0.698)

reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()

for s in files:
    for f in files[s]:
        reader.open(f)
        for event in reader:
            mcps = event.getCollection("MCParticle")
            trks = event.getCollection("SelectedTracks")
            # pfos = event.getCollection("PandoraPFOs")
            # clusters = event.getCollection("PandoraClusters")

            mcp_electrons = []
            for mcp in mcps:
                if mcp.getGeneratorStatus() != 1 or abs(mcp.getPDG()) != 11: continue
                tlv = getTLV(mcp)
                theta = tlv.Theta()
                
                if 0.698 < theta < 2.443:
                    hists[s]["mcp_el_central"].Fill(theta)
                elif 0.0 <= theta <= 0.698:
                    hists[s]["mcp_el_forward"].Fill(theta)
                
                mcp_electrons.append(tlv)

            trk_electrons = [getTrackTLV(t) for t in trks]

            for mcp_el in mcp_electrons:
                m_theta = mcp_el.Theta()
                suffix = ""
                if 0.698 < m_theta < 2.443: suffix = "_central"
                elif 0.0 <= m_theta <= 0.698: suffix = "_forward"
                else: continue

                if any(mcp_el.DeltaR(t) < 0.2 for t in trk_electrons):
                    hists[s]["trk_el_match" + suffix].Fill(m_theta)

        reader.close()

print("\n--- Efficiency Results ---")
for s in hists:
    display_names = {"trk_el_match": "SelectedTracks"}
    
    for region in ["central", "forward"]:
        denom = hists[s]["mcp_el_" + region]
        eff_map = {}
        
        for obj in display_names:
            num = hists[s][obj + "_" + region]
            
            n_num = num.GetEntries()
            n_den = denom.GetEntries()
            efficiency = n_num / n_den if n_den > 0 else 0
            
            region_name = "40-140 deg" if region == "central" else "0-40 deg"
            print(f"Region: {region_name} | {display_names[obj]} Eff: {efficiency:.4f} ({int(n_num)}/{int(n_den)})")
            
            teff = ROOT.TEfficiency(num, denom)
            graph = teff.CreateGraph()
            eff_map[display_names[obj]] = graph
        
        suffix_label = "40-140deg" if region == "central" else "0-40deg"
        plotEfficiencies(eff_map, os.path.join(PLOT_DIR, f"EFFICIENCY_TRACKS_{suffix_label}.png"), 
                         xlabel="#theta [rad]", ylabel="Efficiency")

print(f"\nPlots saved to {PLOT_DIR}")
