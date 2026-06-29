import os
import math
import glob
import pyLCIO
from pyLCIO import EVENT, IOIMPL
import ROOT
import plotHelper

ROOT.gROOT.SetBatch()

PLOT_DIR = "/scratch/jwatts/mucol/mucolstudies/plots2026"
os.makedirs(PLOT_DIR, exist_ok=True)

samples = glob.glob("/scratch/jwatts/mucol/v11Container/reco/electron25Fixed_reco.slcio")
#samples = glob.glob("/scratch/jwatts/mucol/v2.11/reco/electronGun_pT_0_50/10kelectron0to50_recoCalib2.slcio")
files = {"10kelectron0to50": samples}

# Electron ID cuts applied to the best matched cluster
ECAL_ONLY     = False   # HCAL energy must be zero
MAX_WIDTH_MM  = 0.0   # transverse shower width < 25 mm

hists = {}
for s in files:
    hists[s] = {}
    for obj in ["mcp_el", "trk_el_match", "clusters_el_match",
                "el_id_match", "pfo_el_match", "reco_el_match", "reco_el_id_match"]:
        hists[s][obj + "_eta"] = ROOT.TH1F(f"{s}_{obj}_eta", ";#eta;Counts", 50, -2.5, 2.5)

reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()

for s in files:
    for f in files[s]:

        reader.open(f)

        for event in reader:

            mcps     = event.getCollection("MCParticle")
            trks     = event.getCollection("SelectedTracks")
            clusters = event.getCollection("PandoraClusters")
            pfos     = event.getCollection("PandoraPFOs")

            # Truth electrons
            mcp_electrons = []
            for mcp in mcps:
                if mcp.getGeneratorStatus() != 1 or abs(mcp.getPDG()) != 11:
                    continue
                tlv = plotHelper.getTLV(mcp)
                if abs(tlv.Eta()) > 2.4:
                    continue
                hists[s]["mcp_el_eta"].Fill(tlv.Eta())
                mcp_electrons.append(tlv)

            # Reco objects
            trk_tlvs = [plotHelper.getTrackTLV(t) for t in trks]

            # Keep cluster objects alongside their TLVs for cut evaluation
            cluster_pairs = [(c, plotHelper.getClusterTLV(c)) for c in clusters]

            pfo_electrons = [plotHelper.getTLV(p) for p in pfos if abs(p.getType()) == 11]

            for mcp_el in mcp_electrons:

                m_eta = mcp_el.Eta()

                # --- track match: highest pT in cone ---
                trks_in_cone = [t for t in trk_tlvs if mcp_el.DeltaR(t) < 0.2]
                best_trk = max(trks_in_cone, key=lambda t: t.Perp(), default=None)
                has_track = best_trk is not None

                # --- cluster match: highest E in cone ---
                pairs_in_cone = [(c_obj, c_tlv) for c_obj, c_tlv in cluster_pairs
                                 if mcp_el.DeltaR(c_tlv) < 0.2]
                best_pair = max(pairs_in_cone, key=lambda x: x[0].getEnergy(), default=None)
                has_cluster = best_pair is not None

                # --- electron ID cuts on best matched cluster ---
                passes_id = False
                if has_cluster:
                    best_clu_obj, best_clu_tlv = best_pair

                    sub_E  = best_clu_obj.getSubdetectorEnergies()
                    hcal_E = sub_E[1] if len(sub_E) > 1 else 0.0
                    ecal_only_pass = (hcal_E == 0.0)

                    hits  = best_clu_obj.getCalorimeterHits()
                    width = 999.0
                    if len(hits) > 0:
                        axis    = best_clu_tlv.Vect().Unit()
                        sum_E   = 0.0
                        sum_Ed2 = 0.0
                        for hit in hits:
                            pos   = hit.getPosition()
                            hit_E = hit.getEnergy()
                            hit_vec = ROOT.TVector3(pos[0], pos[1], pos[2])
                            d_perp  = (hit_vec - (hit_vec.Dot(axis)) * axis).Mag()
                            sum_Ed2 += hit_E * d_perp**2
                            sum_E   += hit_E
                        if sum_E > 0:
                            width = math.sqrt(sum_Ed2 / sum_E)

                    passes_id = ecal_only_pass and (width < MAX_WIDTH_MM)

                # --- PFO match: closest dR in cone ---
                pfo_drs = [(mcp_el.DeltaR(p), p) for p in pfo_electrons if mcp_el.DeltaR(p) < 0.2]
                best_pfo_dr = min(pfo_drs, key=lambda x: x[0], default=(999, None))
                has_pfo = best_pfo_dr[0] < 0.2

                # --- fill histograms ---
                if has_track:
                    hists[s]["trk_el_match_eta"].Fill(m_eta)
                if has_cluster:
                    hists[s]["clusters_el_match_eta"].Fill(m_eta)
                if passes_id:
                    hists[s]["el_id_match_eta"].Fill(m_eta)
                if has_pfo:
                    hists[s]["pfo_el_match_eta"].Fill(m_eta)
                if has_track and has_cluster:
                    hists[s]["reco_el_match_eta"].Fill(m_eta)
                if has_track and passes_id:
                    hists[s]["reco_el_id_match_eta"].Fill(m_eta)

        reader.close()

# --- Print summary ---
print("\n--- Efficiency Results ---")
for s in hists:

    display_names = {
        "trk_el_match":      "SelectedTracks",
        "clusters_el_match": "PandoraClusters",
        "el_id_match":       "Cluster + Electron ID",
        "pfo_el_match":      "PandoraPFOs (PDGID=11)",
        "reco_el_match":     "Track + Cluster",
        "reco_el_id_match":  "Track + Cluster + Electron ID"
    }

    denom = hists[s]["mcp_el_eta"]
    eff_map = {}

    for obj, label in display_names.items():
        num = hists[s][obj + "_eta"]
        n_num = num.GetEntries()
        n_den = denom.GetEntries()
        eff   = n_num / n_den if n_den > 0 else 0
        print(f"  {label:<35} | Eff: {eff:.4f} ({int(n_num)}/{int(n_den)})")

        if obj in ["pfo_el_match", "reco_el_match", "reco_el_id_match", "el_id_match"]:
            teff  = ROOT.TEfficiency(num, denom)
            graph = teff.CreateGraph()
            eff_map[label] = graph

    plotHelper.plotEfficiencies(
        eff_map,
        os.path.join(PLOT_DIR, "EFFICIENCY_COMBINED_ETA.png"),
        xlabel="#eta",
        ylabel="Efficiency"
    )

print(f"\nPlots saved to {PLOT_DIR}")
