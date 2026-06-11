import os, glob, ROOT, pyLCIO
exec(open("./plotHelper.py").read()); ROOT.gROOT.SetBatch(True); os.makedirs("plots", exist_ok=True)
h_1d = ROOT.TH1F("0_50", ";E/p;Entries", 100, 0, 3)
h_2d = ROOT.TH2F("2d", ";Track p [GeV];Cluster E [GeV]", 100, 0, 200, 100, 0, 200)
r = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()
ev_idx = 0
for f in glob.glob("/scratch/jwatts/mucol/v2.11/reco/electronGun_pT_0_50/10kelectron0to50_recoCalib2.slcio"):

#for f in glob.glob("/scratch/jwatts/mucol/v11Container/reco/10kelectron0to50_reco.slcio"):
#for f in glob.glob("/scratch/jwatts/mucol/v2.11/reco/electronGun_pT_0_50/electronGun_pT_0_50_reco_0.slcio"):
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

                # --- match track: closest dR track within dR < 0.2 ---
                best_trk, best_trk_dr = None, 0.2
                for trk in tracks:
                    trk_tlv = getTrackTLV(trk)
                    dr = trk_tlv.DeltaR(mcp_tlv)
                    if dr < best_trk_dr:
                        best_trk    = trk
                        best_trk_dr = dr

                # --- match cluster: closest dR cluster within dR < 0.2 ---
                best_clu, best_clu_dr = None, 0.2
                for clu in clusters:
                    clu_tlv = getClusterTLV(clu)
                    dr = clu_tlv.DeltaR(mcp_tlv)
                    if dr < best_clu_dr:
                        best_clu  = clu
                        best_clu_dr = dr

                if best_trk is None or best_clu is None: continue
                p = getTrackTLV(best_trk).P()     # Actual Track Momentum magnitude
                E = best_clu.getEnergy()           # Actual Cluster Energy
                if p < 2.0: continue               # Track Momentum cut
                if E < 2.0: continue               # Cluster Energy cut
                h_1d.Fill(E / p); h_2d.Fill(p, E)
                if ev_idx < 100: print(f"Event {ev_idx:02d} | Cluster E: {E:6.2f} GeV | Track p: {p:6.2f} GeV | E/p: {E/p:.3f}")
        except Exception as e: pass
        ev_idx += 1
    r.close()
plotHistograms({"0-50 GeV": h_1d}, "plots2026/electron_ep_ratio_helper.png", "E/p", "Entries", atltext=["Muon Collider", "Simulation, no BIB", "|#eta| < 2.4", "MAIA Detector Concept"])
c = ROOT.TCanvas(); c.SetRightMargin(0.15); ROOT.gStyle.SetOptStat(0); h_2d.Draw("COLZ")
l = ROOT.TLine(0, 0, 200, 200); l.SetLineColor(ROOT.kRed); l.SetLineStyle(2); l.Draw("SAME"); c.SaveAs("plots2026/electron_E_vs_p_2D.png")
