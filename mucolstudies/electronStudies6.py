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
    # Defining pT histograms for the two theta regions
    # Range set to 0-100 GeV as requested
    for obj in ["mcp_el", "trk_el_match"]: 
        hists[s][obj + "_barrel"] = ROOT.TH1F(f"{s}_{obj}_barrel", ";p_{T} [GeV];Counts", 50, 0, 100)
        hists[s][obj + "_endcap"] = ROOT.TH1F(f"{s}_{obj}_endcap", ";p_{T} [GeV];Counts", 50, 0, 100)

reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()

for s in files:
    for f in files[s]:
        reader.open(f)
        for event in reader:
            mcps = event.getCollection("MCParticle")
            trks = event.getCollection("SelectedTracks")

            mcp_electrons = []
            for mcp in mcps:
                if mcp.getGeneratorStatus() != 1 or abs(mcp.getPDG()) != 11: continue
                tlv = getTLV(mcp)
                theta = tlv.Theta()
                pt = tlv.Pt()
                
                # Filter by theta regions and fill pT histograms
                # Barrel: 40-140 deg (0.698 - 2.443 rad)
                if 0.698 < theta < 2.443:
                    hists[s]["mcp_el_barrel"].Fill(pt)
                # Endcap: 0-40 deg (0.0 - 0.698 rad)
                elif 0.0 <= theta <= 0.698:
                    hists[s]["mcp_el_endcap"].Fill(pt)
                
                mcp_electrons.append(tlv)

            trk_electrons = [getTrackTLV(t) for t in trks]

            for mcp_el in mcp_electrons:
                m_theta = mcp_el.Theta()
                m_pt = mcp_el.Pt()
                
                suffix = ""
                if 0.698 < m_theta < 2.443: suffix = "_barrel"
                elif 0.0 <= m_theta <= 0.698: suffix = "_endcap"
                else: continue

                if any(mcp_el.DeltaR(t) < 0.2 for t in trk_electrons):
                    hists[s]["trk_el_match" + suffix].Fill(m_pt)

        reader.close()

print("\n--- Efficiency Results (vs pT) ---")
for s in hists:
    display_names = {"trk_el_match": "SelectedTracks"}
    
    for region in ["barrel", "endcap"]:
        denom = hists[s]["mcp_el_" + region]
        eff_map = {}
        
        for obj in display_names:
            num = hists[s][obj + "_" + region]
            
            n_num = num.GetEntries()
            n_den = denom.GetEntries()
            efficiency = n_num / n_den if n_den > 0 else 0
            
            region_name = "Barrel (40-140 deg)" if region == "barrel" else "Endcap (0-40 deg)"
            print(f"Region: {region_name} | {display_names[obj]} Eff: {efficiency:.4f} ({int(n_num)}/{int(n_den)})")
            
            teff = ROOT.TEfficiency(num, denom)
            graph = teff.CreateGraph()
            eff_map[display_names[obj]] = graph
        
        suffix_label = "BARREL_40-140deg" if region == "barrel" else "ENDCAP_0-40deg"
        plotEfficiencies(eff_map, os.path.join(PLOT_DIR, f"EFFICIENCY_pT_{suffix_label}.png"), 
                         xlabel="p_{T} [GeV]", ylabel="Efficiency")

print(f"\nPlots saved to {PLOT_DIR}")
