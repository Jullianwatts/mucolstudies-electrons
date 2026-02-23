import pyLCIO
import ROOT
import os
import plotHelper

ROOT.gROOT.SetBatch()
fnames = ["/scratch/jwatts/mucol/data/reco/electronGun_pT_0_50/electronGun_pT_0_50_reco_9.slcio"]
out_dir = "plots2026"
if not os.path.exists(out_dir): os.makedirs(out_dir)

hists = {}
total_trk, total_clu, i_evt = 0, 0, 0
for obj in ["track", "cluster"]:
    v_n = plotHelper.variables["evt"]["n"]
    hists[f"{obj}_n"] = ROOT.TH1F(f"{obj}_n", f"{obj}_n", v_n["nbins"], v_n["xmin"], v_n["xmax"])
    for var, cfg in plotHelper.variables["obj"].items():
        hists[f"{obj}_{var}"] = ROOT.TH1F(f"{obj}_{var}", f"{obj}_{var}", cfg["nbins"], cfg["xmin"], cfg["xmax"])

reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()
for f in fnames:
    reader.open(f)
    for event in reader:
        trk_coll = event.getCollection("SelectedTracks")
        n_trk = 0
        for trk in trk_coll:
            tlv = plotHelper.getTrackTLV(trk, m=0.000511)
            if abs(tlv.Eta()) <= 2.4: # and tlv.Pt() > 10:
                plotHelper.fillObjHists(hists, "track", tlv)
                n_trk += 1
        hists["track_n"].Fill(n_trk)

        clu_coll = event.getCollection("PandoraClusters")
        n_clu = 0
        for clu in clu_coll:
            tlv = plotHelper.getClusterTLV(clu)
            if abs(tlv.Eta()) <= 2.4: # and tlv.E() > 10:
                plotHelper.fillObjHists(hists, "cluster", tlv)
                n_clu += 1
        hists["cluster_n"].Fill(n_clu)

        print(f"Evt {i_evt:3} | Tracks: {n_trk:2} | Clusters: {n_clu:2}")
        total_trk += n_trk; total_clu += n_clu; i_evt += 1
    reader.close()

if i_evt > 0:
    print(f"\n{'='*40}\nAVERAGES (Eta < 2.4):")
    print(f"Tracks/Evt: {total_trk/i_evt:.2f} | Clusters/Evt: {total_clu/i_evt:.2f}\n{'='*40}")

for obj in ["track", "cluster"]:
    plotHelper.plotHistograms({f"{obj}_n": hists[f"{obj}_n"]}, f"{out_dir}/{obj}_n.png", xlabel=plotHelper.variables["evt"]["n"]["label"], ylabel="Events", atltext=["Muon Collider", "Simulation, no BIB", "|#eta| < 2.4", "MAIA Concept"])
    for var, cfg in plotHelper.variables["obj"].items():
        plotHelper.plotHistograms({f"{obj}_{var}": hists[f"{obj}_{var}"]}, f"{out_dir}/{obj}_{var}.png", xlabel=cfg["label"], ylabel="Entries", atltext=["Muon Collider", "Simulation, no BIB", "|#eta| < 2.4", "MAIA Concept"])

print(f"Done. Results in {out_dir}/")
