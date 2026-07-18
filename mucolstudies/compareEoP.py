import os, glob, ROOT, pyLCIO
exec(open("./plotHelper.py").read()); ROOT.gROOT.SetBatch(True); os.makedirs("plots", exist_ok=True)

samples = {
    "Sample 1": "/scratch/jwatts/mucol/v2.11/reco/electronGun_pT_0_50/electronGun_pT_0_50_reco_Final.slcio",
    "Sample 2": "/scratch/jwatts/mucol/v11Container/reco/10kelectron0to50_reco.slcio",
}

h_1d_map = {}
h_2d = ROOT.TH2F("2d", ";Track p [GeV];Cluster E [GeV]", 100, 0, 200, 100, 0, 200)
r = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()

for i_sample, (sample_name, sample_path) in enumerate(samples.items()):
    h_1d = ROOT.TH1F(sample_name, ";E/p;Entries", 100, 0, 3)
    h_1d_map[sample_name] = h_1d
    ev_idx = 0
    for f in glob.glob(sample_path):
        if ev_idx >= 10000: break
        r.open(f)
        for ev in r:
            if ev_idx >= 10000: break
            try:
                tracks   = list(ev.getCollection("SelectedTracks"))
                clusters = list(ev.getCollection("PandoraClusters"))
                for mcp in ev.getCollection("MCParticle"):
                    if abs(mcp.getPDG()) != 11: continue
                    mcp_tlv = getTLV(mcp)
                    if mcp_tlv.E() < 10.0: continue          # PFO Energy cut (kept as-is)

                    # --- match track: highest-pT track within dR < 0.2 ---
                    best_trk, best_trk_pt = None, 0.0
                    for trk in tracks:
                        trk_tlv = getTrackTLV(trk)
                        if trk_tlv.DeltaR(mcp_tlv) < 0.2:
                            if trk_tlv.Perp() > best_trk_pt:
                                best_trk    = trk
                                best_trk_pt = trk_tlv.Perp()

                    # --- match cluster: highest-pT (E) cluster within dR < 0.2 ---
                    best_clu, best_clu_e = None, 0.0
                    for clu in clusters:
                        clu_tlv = getClusterTLV(clu)
                        if clu_tlv.DeltaR(mcp_tlv) < 0.2:
                            if clu_tlv.Perp() > best_clu_e:
                                best_clu  = clu
                                best_clu_e = clu_tlv.Perp()

                    if best_trk is None or best_clu is None: continue
                    p = getTrackTLV(best_trk).P()     # Actual Track Momentum magnitude
                    E = best_clu.getEnergy()           # Actual Cluster Energy
                    if p < 2.0: continue               # Track Momentum cut
                    if E < 2.0: continue               # Cluster Energy cut
                    h_1d.Fill(E / p)
                    if i_sample == 0: h_2d.Fill(p, E)
                    if ev_idx < 100: print(f"{sample_name} | Event {ev_idx:02d} | Cluster E: {E:6.2f} GeV | Track p: {p:6.2f} GeV | E/p: {E/p:.3f}")
            except Exception as e: pass
            ev_idx += 1
        r.close()

plotHistograms(h_1d_map, "plots2026/electron_ep_ratio_helper.png", "E/p", "Entries", atltext=["Muon Collider", "Simulation, no BIB", "|#eta| < 2.4", "MAIA Detector Concept"])
c = ROOT.TCanvas(); c.SetRightMargin(0.15); ROOT.gStyle.SetOptStat(0); h_2d.Draw("COLZ")
l = ROOT.TLine(0, 0, 200, 200); l.SetLineColor(ROOT.kRed); l.SetLineStyle(2); l.Draw("SAME"); c.SaveAs("plots2026/electron_E_vs_p_2D.png")
