import math
import glob
import ROOT
import pyLCIO
import os

exec(open("./plotHelper.py").read())
ROOT.gROOT.SetBatch()

PLOT_DIR = "/scratch/jwatts/mucol/mucolstudies/plots2026"
os.makedirs(PLOT_DIR, exist_ok=True)

el_samples = glob.glob("/scratch/jwatts/mucol/data/reco/electronGun_pT_0_50/electronGun_pT_0_50_reco_7.slcio")
mu_path = "/ospool/uc-shared/project/futurecolliders/gpenn/v8_muons/no_BIB_reco"
mu_samples = glob.glob(os.path.join(mu_path, "0p5_10GeV/*.slcio")) + \
             glob.glob(os.path.join(mu_path, "10_50_GeV/*.slcio"))

files = {
    "electronGun_pT_0_50": el_samples,
    "muonGun": mu_samples
}

hists = {}
for s in files:
    hists[s] = {}
    prefix = "el" if "electron" in s else "mu"
    for obj in [f"mcp_{prefix}", f"trk_{prefix}_match"]: 
        hists[s][obj + "_barrel"] = ROOT.TH1F(f"{s}_{obj}_barrel", ";p_{T} [GeV];Counts", 50, 0, 100)
        hists[s][obj + "_endcap"] = ROOT.TH1F(f"{s}_{obj}_endcap", ";p_{T} [GeV];Counts", 50, 0, 100)

reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()

for s in files:
    current_pdg = 11 if "electron" in s else 13
    prefix = "el" if current_pdg == 11 else "mu"
    
    for f in files[s]:
        reader.open(f)
        for event in reader:
            mcps = event.getCollection("MCParticle")
            trks = event.getCollection("SelectedTracks")

            mcp_list = []
            for mcp in mcps:
                if mcp.getGeneratorStatus() != 1: continue
                if abs(mcp.getPDG()) != current_pdg: continue
                
                tlv = getTLV(mcp)
                theta = tlv.Theta()
                pt = tlv.Pt()

                if pt < 10: continue
                
                if 0.698 < theta < 2.443:
                    hists[s][f"mcp_{prefix}_barrel"].Fill(pt)
                    mcp_list.append(tlv)
                elif (0.0 <= theta <= 0.698) or (theta >= 2.443):
                    hists[s][f"mcp_{prefix}_endcap"].Fill(pt)
                    mcp_list.append(tlv)

            trk_tlvs = []
            for t in trks:
                tlv = getTrackTLV(t)
                if tlv.Pt() >= 2:
                    trk_tlvs.append(tlv)

            for mcp_tlv in mcp_list:
                m_theta = mcp_tlv.Theta()
                m_pt = mcp_tlv.Pt()
                
                suffix = "_barrel" if 0.698 < m_theta < 2.443 else "_endcap"

                if any(mcp_tlv.DeltaR(t) < 0.1 for t in trk_tlvs):
                    hists[s][f"trk_{prefix}_match" + suffix].Fill(m_pt)

        reader.close()

print("\n--- Efficiency Results (vs pT) ---")
for region in ["barrel", "endcap"]:
    eff_map = {}
    region_label = "BARREL_40-140deg" if region == "barrel" else "ENDCAP_0-40deg"
    
    for s in files:
        prefix = "el" if "electron" in s else "mu"
        label = "Electrons" if prefix == "el" else "Muons"
        
        num = hists[s][f"trk_{prefix}_match_{region}"]
        denom = hists[s][f"mcp_{prefix}_{region}"]
        
        n_num = num.GetEntries()
        n_den = denom.GetEntries()
        efficiency = n_num / n_den if n_den > 0 else 0
        
        print(f"Region: {region} | {label} Eff: {efficiency:.4f} ({int(n_num)}/{int(n_den)})")
        
        teff = ROOT.TEfficiency(num, denom)
        graph = teff.CreateGraph()
        eff_map[label] = graph

    plotEfficiencies(eff_map, os.path.join(PLOT_DIR, f"EFFICIENCY_COMPARISON_pT_{region_label}.png"), 
                     xlabel="p_{T} [GeV]", ylabel="Efficiency")

