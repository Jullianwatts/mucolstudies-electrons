import pyLCIO, ROOT, os, plotHelper, math
ROOT.gROOT.SetBatch()

calib_factor = 1.08 

config = {
    "v2.9.7": {"path": "/scratch/jwatts/mucol/data/reco/electronGun_pT_0_50/electronGun_pT_0_50_reco_5.slcio", "trk_coll": "SiTracks_Refitted"},
    "v11": {"path": "/scratch/jwatts/mucol/data/reco/electronGun_pT_0_50/electronGun_pT_0_50_reco_9.slcio", "trk_coll": "SelectedTracks"}
}

out_dir = "plots2026"
if not os.path.exists(out_dir): os.makedirs(out_dir)

def getDR(v1, v2): return math.sqrt((v1.Eta()-v2.Eta())**2 + v1.DeltaPhi(v2)**2)

hists, stats = {v: {} for v in config}, {v: {"evts": 0, "m_trk": 0, "m_clu": 0, "t_trk": 0, "t_clu": 0, "sum_ptrue": 0, "sum_precon": 0, "sum_etrue": 0, "sum_erecon": 0} for v in config}

for ver, cfg in config.items():
    print(f"\n--- {ver} ---"); reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()
    
    hists[ver]["e_res"] = ROOT.TH1F(f"e_res_{ver}", "", 100, -1.2, 0.5); hists[ver]["e_res"].Sumw2()
    hists[ver]["fail_eta"] = ROOT.TH1F(f"fail_eta_{ver}", "Truth Eta of Failures", 50, -2.5, 2.5)
    hists[ver]["fail_e"] = ROOT.TH1F(f"fail_e_{ver}", "Truth Energy of Failures", 50, 0, 100)
    
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
        stats[ver]["sum_ptrue"] += truth_tlv.Pt()
        stats[ver]["sum_etrue"] += truth_tlv.E()

        trk_coll = event.getCollection(cfg['trk_coll']); b_trk, min_dr_t, n_t = None, 0.1, 0
        for trk in trk_coll:
            tlv = plotHelper.getTrackTLV(trk)
            if abs(tlv.Eta()) < 2.4:
                n_t += 1; dr = getDR(truth_tlv, tlv)
                if dr < min_dr_t: min_dr_t, b_trk = dr, tlv

        clu_coll = event.getCollection("PandoraClusters"); b_clu, min_dr_c, n_c = None, 0.2, 0
        for clu in clu_coll:
            tlv = plotHelper.getClusterTLV(clu)
            if abs(tlv.Eta()) < 2.4 and tlv.E() > 2:
                n_c += 1; dr = getDR(truth_tlv, tlv)
                if dr < min_dr_c:
                    min_dr_c, b_clu = dr, tlv

        e_reco = (b_clu.E() * calib_factor) if b_clu else 0
        hists[ver]["e_res"].Fill((e_reco - truth_tlv.E()) / truth_tlv.E())

        if not b_clu:
            hists[ver]["fail_eta"].Fill(truth_tlv.Eta())
            hists[ver]["fail_e"].Fill(truth_tlv.E())

        if b_trk:
            plotHelper.fillObjHists(hists[ver], "track", b_trk); stats[ver]["m_trk"] += 1
            stats[ver]["sum_precon"] += b_trk.Pt()
        if b_clu:
            plotHelper.fillObjHists(hists[ver], "cluster", b_clu); stats[ver]["m_clu"] += 1
            stats[ver]["sum_erecon"] += (b_clu.E() * calib_factor)

        hists[ver]["track_n"].Fill(n_t); hists[ver]["cluster_n"].Fill(n_c)
        stats[ver]["t_trk"] += n_t; stats[ver]["t_clu"] += n_c

    reader.close()

print(f"\n{'-'*10}\nFINAL AVERAGES (Calib: {calib_factor}, E_true > 10, E_reco > 2)\n{'-'*10}")
for v in config:
    n = stats[v]["evts"]
    if n > 0:
        print(f"{v}: Eff_T {stats[v]['m_trk']/n*100:.1f}%, Eff_C {stats[v]['m_clu']/n*100:.1f}%")
        print(f"   Residuals: pT_sum_miss: {stats[v]['sum_ptrue']-stats[v]['sum_precon']:.2f} GeV | E_sum_miss: {stats[v]['sum_etrue']-stats[v]['sum_erecon']:.2f} GeV")

res_map = {v: hists[v]["e_res"] for v in config}
plotHelper.plotHistograms(res_map, f"{out_dir}/energy_response_calibrated.png", xlabel="(E_reco_calib - E_true) / E_true", ylabel="Events", atltext=["Muon Collider", "Comparison", "Calibrated"])

for feat in ["fail_eta", "fail_e"]:
    f_map = {v: hists[v][feat] for v in config}
    plotHelper.plotHistograms(f_map, f"{out_dir}/{feat}.png", xlabel=feat, ylabel="Failure Counts", atltext=["Muon Collider", "Failures Only"])

print(f"\nDone. Results in {out_dir}/")
