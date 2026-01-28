import math
import glob
import ROOT
import pyLCIO
import os

exec(open("./plotHelper.py").read())
ROOT.gROOT.SetBatch()

PLOT_DIR = "/scratch/jwatts/mucol/mucolstudies/plots2026"
os.makedirs(PLOT_DIR, exist_ok=True)
file_path = "/scratch/jwatts/mucol/data/reco/electronGun_pT_0_50/electronGun_pT_0_50_reco_4.slcio"

h_sum_ratios = ROOT.TH1F("h_sum_ratios", "Sum of E/p Ratios", 48, -2.4, 2.4)
h_count = ROOT.TH1F("h_count", "Particle Count", 48, -2.4, 2.4)

reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()

total_mc_electrons = 0
matched_electrons = 0

try:
    reader.open(file_path)
    print(f"Reading file: {file_path}")

    for event in reader:
        mcps = event.getCollection("MCParticle")
        trks = event.getCollection("SiTracks")
        clusters = event.getCollection("PandoraClusters")

        for mcp in mcps:
            if mcp.getGeneratorStatus() != 1 or abs(mcp.getPDG()) != 11: 
                continue
            
            mcp_tlv = getTLV(mcp)
            mcp_eta = mcp_tlv.Eta()
            
            if mcp_tlv.E() < 10 or abs(mcp_eta) > 2.4: 
                continue

            total_mc_electrons += 1

            best_trk_p = -1
            min_trk_dr = 0.2 
            for t in trks:
                trk_tlv = getTrackTLV(t)
                dr = mcp_tlv.DeltaR(trk_tlv)
                if dr < min_trk_dr:
                    min_trk_dr = dr
                    best_trk_p = trk_tlv.P() 

            best_clu_e = -1
            min_clu_dr = 0.2
            for c in clusters:
                clu_tlv = getClusterTLV(c)
                dr = mcp_tlv.DeltaR(clu_tlv)
                if dr < min_clu_dr:
                    min_clu_dr = dr
                    best_clu_e = c.getEnergy() 

            if best_trk_p > 0 and best_clu_e > 0:
                matched_electrons += 1
                ep_ratio = best_clu_e / best_trk_p
                
                h_sum_ratios.Fill(mcp_eta, ep_ratio)
                h_count.Fill(mcp_eta, 1)

    reader.close()
except Exception as e: print(f"Error: {e}")

# 5. Calculate the final Average(E/p) as a function of Eta
h_avg_ep = h_sum_ratios.Clone("h_avg_ep")
h_avg_ep.Divide(h_count)

h_map = {"Average E/p per Particle": h_avg_ep}
plotHistograms(
    h_map, 
    os.path.join(PLOT_DIR, "ecal_avg_ep_vs_eta.png"), 
    xlabel="#eta", 
    ylabel="Mean(E_{cluster} / p_{track})", 
    atltext=["Muon Collider", "Simulation, no BIB", "|#eta| < 2.4", "MAIA Detector Concept"]
)

# PRINT STATISTICS
print("ECAL CALIBRATION STATISTICS")
print(f"Total MC Electrons:  {total_mc_electrons}")
print(f"Matched Electrons:   {matched_electrons}")
if total_mc_electrons > 0:
    print(f"Matching Efficiency: {100.0*matched_electrons/total_mc_electrons:.1f}%")
