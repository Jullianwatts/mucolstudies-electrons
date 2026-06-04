import math
import glob
import ROOT
import pyLCIO
import os

exec(open("./plotHelper.py").read())
ROOT.gROOT.SetBatch()

PLOT_DIR = "/scratch/jwatts/mucol/mucolstudies/plots2026"
os.makedirs(PLOT_DIR, exist_ok=True)

samples = glob.glob("/scratch/jwatts/mucol/v11Container/reco/1kelectron25_reco.slcio")
#samples = glob.glob("/scratch/jwatts/mucol/v2.11/reco/electronGun_pT_0_50/electronGun_pT_0_50_reco_0.slcio")
#samples = glob.glob("/scratch/jwatts/mucol/data/reco/electronGun_pT_0_50/electronGun_pT_0_50_reco_9.slcio")
#files = {"electronGun_pT_0_50": samples}
files = {"1kelectron25": samples}

hists = {}
for s in files:
    hists[s] = {}
    for obj in ["mcp_el", "trk_el_match", "clusters_el_match", "pfo_el_match", "reco_el_match"]:
        # Updated range for eta plots
        hists[s][obj + "_eta"] = ROOT.TH1F(f"{s}_{obj}_eta", ";#eta;Counts", 50, -2.5, 2.5)

reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()

for s in files:
    for f in files[s]:
        reader.open(f)
        for event in reader:
            mcps = event.getCollection("MCParticle")
            trks = event.getCollection("SelectedTracks")
            clusters = event.getCollection("PandoraClusters")
            pfos = event.getCollection("PandoraPFOs")

            mcp_electrons = []
            for mcp in mcps:
                if mcp.getGeneratorStatus() != 1 or abs(mcp.getPDG()) != 11: continue
                tlv = getTLV(mcp)
                eta = tlv.Eta()

                # Apply the strict window cut: |eta| <= 2.4
                if abs(eta) > 2.4: continue

                hists[s]["mcp_el_eta"].Fill(eta)
                mcp_electrons.append(tlv)

            trk_electrons = [getTrackTLV(t) for t in trks]
            clu_electrons = [getClusterTLV(c) for c in clusters]
            pfo_electrons = [getTLV(p) for p in pfos if abs(p.getType()) == 11]

            for mcp_el in mcp_electrons:
                m_eta = mcp_el.Eta()

                best_trk_dr = min([mcp_el.DeltaR(t) for t in trk_electrons], default=999)
                best_clu_dr = min([mcp_el.DeltaR(c) for c in clu_electrons], default=999)
                best_pfo_dr = min([mcp_el.DeltaR(p) for p in pfo_electrons], default=999)

                has_track = best_trk_dr < 0.2
                has_cluster = best_clu_dr < 0.2
                has_pfo = best_pfo_dr < 0.2

                if has_track:
                    hists[s]["trk_el_match_eta"].Fill(m_eta)

                if has_cluster:
                    hists[s]["clusters_el_match_eta"].Fill(m_eta)

                if has_pfo:
                    hists[s]["pfo_el_match_eta"].Fill(m_eta)

                if has_track and has_cluster:
                    hists[s]["reco_el_match_eta"].Fill(m_eta)

        reader.close()

print("\n--- Efficiency Results ---")
for s in hists:
    # Keeps printing individual metrics to terminal, but limits the canvas to PFOs and Track + Cluster
    display_names = {
        "trk_el_match": "SelectedTracks",
        "clusters_el_match": "PandoraClusters",
        "pfo_el_match": "PandoraPFOs (PDGID=11)",
        "reco_el_match": "Track + Cluster"
    }

    eff_map = {}
    denom = hists[s]["mcp_el_eta"]

    for obj in display_names:
        num = hists[s][obj + "_eta"]

        n_num = num.GetEntries()
        n_den = denom.GetEntries()
        efficiency = n_num / n_den if n_den > 0 else 0

        print(f"Object: {display_names[obj]:<25} | Eff: {efficiency:.4f} ({int(n_num)}/{int(n_den)})")

        # Only add PFO and combined Track + Cluster to the final canvas dictionary
        if obj in ["pfo_el_match", "reco_el_match"]:
            teff = ROOT.TEfficiency(num, denom)
            graph = teff.CreateGraph()
            eff_map[display_names[obj]] = graph

    plotEfficiencies(eff_map, os.path.join(PLOT_DIR, "EFFICIENCY_COMBINED_ETA.png"),
                     xlabel="#eta", ylabel="Efficiency")

print(f"\nPlots saved to {PLOT_DIR}")
