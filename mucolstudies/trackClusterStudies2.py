import pyLCIO, ROOT, os, plotHelper, math
ROOT.gROOT.SetBatch()

config = {
    "v2.9.7": {"path": "/scratch/jwatts/mucol/data/reco/electronGun_pT_0_50/electronGun_pT_0_50_reco_5.slcio", "trk_coll": "SiTracks_Refitted"},
    "v11": {"path": "/scratch/jwatts/mucol/data/reco/electronGun_pT_0_50/electronGun_pT_0_50_reco_9.slcio", "trk_coll": "SelectedTracks"}
}

out_dir = "plots2026"
if not os.path.exists(out_dir): os.makedirs(out_dir)

def getDR(v1, v2): return math.sqrt((v1.Eta()-v2.Eta())**2 + v1.DeltaPhi(v2)**2)

hists, stats = {v: {} for v in config}, {v: {"evts": 0, "m_trk": 0, "m_clu": 0, "t_trk": 0, "t_clu": 0} for v in config}

for ver, cfg in config.items():
    print(f"\n--- {ver} ---"); reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()
    for obj in ["track", "cluster"]:
        hists[ver][f"{obj}_n"] = ROOT.TH1F(f"{obj}_n_{ver}", "", 20, 0, 20); hists[ver][f"{obj}_n"].Sumw2()
        for var, c in plotHelper.variables["obj"].items(): 
            hists[ver][f"{obj}_{var}"] = ROOT.TH1F(f"{obj}_{var}_{ver}", "", c["nbins"], c["xmin"], c["xmax"]); hists[ver][f"{obj}_{var}"].Sumw2()
    
    try: reader.open(cfg['path'])
    except: print(f"Skip {ver}"); continue

    for i_evt, event in enumerate(reader):
        mc_coll = event.getCollection("MCParticle"); truth_tlv = None
        for mcp in mc_coll:
            if abs(mcp.getPDG()) == 11 and mcp.getGeneratorStatus() == 1 and mcp.getEnergy() > 10:
                p = mcp.getMomentum(); truth_tlv = ROOT.TLorentzVector(); truth_tlv.SetPxPyPzE(p[0], p[1], p[2], mcp.getEnergy()); break
        if not truth_tlv: continue
        stats[ver]["evts"] += 1

        # TRACKS: No Pt cut, only Eta
        trk_coll = event.getCollection(cfg['trk_coll']); b_trk, min_dr_t, n_t = None, 0.1, 0
        for trk in trk_coll:
            tlv = plotHelper.getTrackTLV(trk)
            if abs(tlv.Eta()) < 2.4:
                n_t += 1; dr = getDR(truth_tlv, tlv)
                if dr < min_dr_t: min_dr_t, b_trk = dr, tlv
        
        # CLUSTERS: Small reco cut (E > 1 GeV)
        clu_coll = event.getCollection("PandoraClusters"); b_clu, min_dr_c, n_c = None, 0.2, 0
        for clu in clu_coll:
            tlv = plotHelper.getClusterTLV(clu)
            if abs(tlv.Eta()) < 2.4 and tlv.E() > 1:
                n_c += 1; dr = getDR(truth_tlv, tlv)
                if dr < min_dr_c: min_dr_c, b_clu = dr, tlv

        if b_trk: plotHelper.fillObjHists(hists[ver], "track", b_trk); stats[ver]["m_trk"] += 1
        if b_clu: plotHelper.fillObjHists(hists[ver], "cluster", b_clu); stats[ver]["m_clu"] += 1
        hists[ver]["track_n"].Fill(n_t); hists[ver]["cluster_n"].Fill(n_c)
        stats[ver]["t_trk"] += n_t; stats[ver]["t_clu"] += n_c

        if i_evt % 50 == 0:
            e_r, pt_r = (b_clu.E() if b_clu else 0), (b_trk.Pt() if b_trk else 0)
            print(f"Evt {i_evt:3} | Trks: {n_t} | Clus: {n_c} | E_miss: {truth_tlv.E()-e_r:.2f} | pT_miss: {truth_tlv.Pt()-pt_r:.2f}")
    reader.close()

print(f"\n{'='*50}\nFINAL AVERAGES (Truth E > 10 GeV)\n{'='*50}")
for v in config:
    n = stats[v]["evts"]
    if n > 0: print(f"{v}: Avg Trk {stats[v]['t_trk']/n:.2f}, Avg Clu {stats[v]['t_clu']/n:.2f}, Eff_T {stats[v]['m_trk']/n*100:.1f}%, Eff_C {stats[v]['m_clu']/n*100:.1f}%")

for obj in ["track", "cluster"]:
    m_map = {v: hists[v][f"{obj}_n"] for v in config}
    plotHelper.plotHistograms(m_map, f"{out_dir}/{obj}_n.png", xlabel="Number per Event", ylabel="Events", atltext=["Muon Collider", "Comparison", "|#eta| < 2.4", "MAIA"])
    for var, c in plotHelper.variables["obj"].items():
        v_map = {v: hists[v][f"{obj}_{var}"] for v in config}
        plotHelper.plotHistograms(v_map, f"{out_dir}/{obj}_{var}.png", xlabel=c["label"], ylabel="Entries", atltext=["Muon Collider", "Comparison", "|#eta| < 2.4", "MAIA"])

print(f"\nDone. Results in {out_dir}/")
