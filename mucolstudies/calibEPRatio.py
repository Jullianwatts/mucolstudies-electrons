import os
import pyLCIO
import ROOT
from pyLCIO import IOIMPL
import plotHelper

#file_path = "/scratch/jwatts/mucol/v2.11/reco/electronGun_pT_0_50/electronGun_pT_0_50_reco_0.slcio"
#file_path = "/scratch/jwatts/mucol/v2.11/reco/electronGun_pT_0_50/10kelectron0to50_recoCalib2.slcio"
file_path = "/scratch/jwatts/mucol/v11Container/reco/10kelectron0to50_reco.slcio"
output_dir = "/scratch/jwatts/mucol/mucolstudies/plots2026"
output_file = os.path.join(output_dir, "ep_ratio_0_50_GeV.png")

hists = {
    11: ROOT.TH1F("Electrons", "Electrons", 60, 0.0, 3.0)
}

h_calc_cluster = ROOT.TH1F("h_calc_cluster", "", 1, 0.0, 1000.0)
h_calc_track = ROOT.TH1F("h_calc_track", "", 1, 0.0, 1000.0)
h_calc_ep = ROOT.TH1F("h_calc_ep", "", 1, 0.0, 1000.0)

event_count = 0
signal_mcps = 0
matched_mcps = 0
n_high_ep = 0
n_within_02 = 0
n_within_04 = 0

reader = IOIMPL.LCFactory.getInstance().createLCReader()
reader.open(file_path)

for event in reader:
    event_count += 1
    try:
        mc_coll = event.getCollection("MCParticle")
        track_coll = event.getCollection("SelectedTracks")
        cluster_coll = event.getCollection("PandoraClusters")
    except:
        continue

    mc_tlv = None
    for mc_part in mc_coll:
        if abs(mc_part.getPDG()) == 11 and mc_part.getGeneratorStatus() == 1:
            if mc_part.getEnergy() > 10.0:
                signal_mcps += 1
                mc_tlv = plotHelper.getTLV(mc_part)
                break

    if not mc_tlv:
        continue

    best_trk_pt = -1.0
    matched_trk_tlv = None
    for trk in track_coll:
        trk_tlv = plotHelper.getTrackTLV(trk)
        if trk_tlv.Pt() > 2.0 and abs(trk_tlv.Eta()) < 2.4:
            if mc_tlv.DeltaR(trk_tlv) < 0.2 and trk_tlv.Pt() > best_trk_pt:
                best_trk_pt = trk_tlv.Pt()
                matched_trk_tlv = trk_tlv

    best_cls_pt = -1.0
    matched_cls_tlv = None
    for clus in cluster_coll:
        cls_tlv = plotHelper.getClusterTLV(clus)
        if cls_tlv.E() > 2.0:
            if mc_tlv.DeltaR(cls_tlv) < 0.2 and cls_tlv.Pt() > best_cls_pt:
                best_cls_pt = cls_tlv.Pt()
                matched_cls_tlv = cls_tlv

    if matched_trk_tlv and matched_cls_tlv:
        matched_mcps += 1
        p = matched_trk_tlv.P()
        if p > 0:
            ep_ratio = matched_cls_tlv.E() / p
            hists[11].Fill(ep_ratio)

            if abs(ep_ratio - 1.0) <= 0.2:
                n_within_02 += 1
            if abs(ep_ratio - 1.0) <= 0.4:
                n_within_04 += 1

            if ep_ratio > 1.0:
                n_high_ep += 1
                h_calc_cluster.Fill(matched_cls_tlv.E())
                h_calc_track.Fill(matched_trk_tlv.Pt())
                h_calc_ep.Fill(ep_ratio)
                print(f"Event {event_count}: E/p = {ep_ratio:.4f}, MCP Eta = {mc_tlv.Eta():.4f}, Cluster Eta = {matched_cls_tlv.Eta():.4f}, Cluster E = {matched_cls_tlv.E():.4f} GeV, Track pT = {matched_trk_tlv.Pt():.4f} GeV")

    if event_count == 10000:
        break

reader.close()
os.makedirs(output_dir, exist_ok=True)

plotHelper.plotHistograms(
    h_map={"Electrons (0-50 GeV)": hists[11]},
    save_name=output_file,
    xlabel="E/p",
    ylabel="Entries",
    interactive=False,
    logy=False,
    atltext=["Muon Collider", "Simulation, no BIB", "|#eta| < 2.4", "MAIA Detector Concept"]
)

avg_cluster_e = h_calc_cluster.GetMean()
err_cluster_e = h_calc_cluster.GetMeanError()
avg_track_pt = h_calc_track.GetMean()
err_track_pt = h_calc_track.GetMeanError()
avg_ep_ratio = h_calc_ep.GetMean()
err_ep_ratio = h_calc_ep.GetMeanError()

print("\n--- Summary Statement ---")
print(f"Total processed events: {event_count}")
print(f"Total primary signal electrons (>10 GeV): {signal_mcps}")
print(f"Total primary signal electrons matched to both a track and cluster: {matched_mcps}")
print(f"Total events with E/p > 1: {n_high_ep}")
print(f"Average Cluster Energy (for E/p > 1): {avg_cluster_e:.4f} +/- {err_cluster_e:.4f} GeV")
print(f"Average Track pT (for E/p > 1): {avg_track_pt:.4f} +/- {err_track_pt:.4f} GeV")
print(f"Average E/p (for E/p > 1): {avg_ep_ratio:.4f} +/- {err_ep_ratio:.4f}")
print(f"Total events with E/p within 1.0 +/- 0.2: {n_within_02}")
print(f"Total events with E/p within 1.0 +/- 0.4: {n_within_04}")
