import os
import pyLCIO
import ROOT
from pyLCIO import IOIMPL
import plotHelper

file_path = "/scratch/jwatts/mucol/data/reco/electronGun_pT_0_50/electronGun_pT_0_50_reco_9.slcio"

output_dir = "/scratch/jwatts/mucol/mucolstudies/plots2026"
output_file = os.path.join(output_dir, "ep_ratio_0_50_GeV.png")

hists = {
    11: ROOT.TH1F("Electrons", "Electrons", 60, 0.0, 3.0),
    211: ROOT.TH1F("Pions", "Pions", 60, 0.0, 3.0)
}

event_count = 0
reader = IOIMPL.LCFactory.getInstance().createLCReader()
reader.open(file_path)

for event in reader:
    event_count += 1
    try:
        pfo_coll = event.getCollection("PandoraPFOs")
    except:
        continue

    for pfo in pfo_coll:
        p_type = abs(pfo.getType())
        if p_type not in hists:
            continue
            
        if len(pfo.getTracks()) == 0 or len(pfo.getClusters()) == 0:
            continue
            
        track_tlv = plotHelper.getTrackTLV(pfo.getTracks()[0])
        cluster_tlv = plotHelper.getClusterTLV(pfo.getClusters()[0])
        
        if abs(track_tlv.Eta()) >= 2.4:
            continue
            
        p = track_tlv.P()
        if p > 0:
            hists[p_type].Fill(cluster_tlv.E() / p)

    if event_count == 1000:
        break

reader.close()

plotHelper.plotHistograms(
    h_map={"Electrons (0-50 GeV)": hists[11], "Pions (0-50 GeV)": hists[211]},
    save_name=output_file,
    xlabel="E/p",
    ylabel="Entries",
    interactive=False,
    logy=False,
    atltext=["Muon Collider", "Simulation, no BIB", "|#eta| < 2.4", "MAIA Detector Concept"]
)
