import os
import math
# Samples: swap in fixed theta=pi/2, pt=25 files when ready
import pyLCIO
from pyLCIO import EVENT, IOIMPL
import ROOT
import plotHelper

PLOT_DIR = "/scratch/jwatts/mucol/mucolstudies/plots2026/pionshowers"
os.makedirs(PLOT_DIR, exist_ok=True)

# Samples: PDG ID is used to filter PFOs — compares what Pandora actually called each particle
samples = {
    "Electrons": {"path": "/scratch/jwatts/mucol/v11Container/reco/electron25Fixed_reco.slcio", "pdg": 11},
    "Pions":     {"path": "/scratch/jwatts/mucol/v11Container/reco/pion25Fixed_reco.slcio",     "pdg": 211}
}

max_events = 1000

shower_vars = {
    "clus_E":    {"nbins": 40, "xmin": 0,    "xmax": 40,   "label": "Leading Cluster Energy [GeV]"},
    "n_clus":    {"nbins": 10, "xmin": 0,    "xmax": 10,   "label": "Number of Clusters per Event"},
    "ecal_frac": {"nbins": 25, "xmin": 0,    "xmax": 1.05, "label": "ECAL Energy Fraction"},
    "hcal_E":    {"nbins": 30, "xmin": 0,    "xmax": 30,   "label": "HCAL Energy [GeV]"},
    "n_hits":    {"nbins": 30, "xmin": 0,    "xmax": 300,  "label": "Number of Hits in Leading Cluster"},
    "start_r":   {"nbins": 30, "xmin": 1400, "xmax": 4000, "label": "Shower Start Radius [mm]"},
    "depth_r":   {"nbins": 30, "xmin": 1400, "xmax": 4000, "label": "Energy-Weighted Shower Depth [mm]"},
    "width":     {"nbins": 30, "xmin": 0,    "xmax": 150,  "label": "Transverse Shower Width [mm]"},
    "E_over_p":  {"nbins": 30, "xmin": 0,    "xmax": 3.0,  "label": "E_{clus}/p_{trk}"}
}

# Book histograms
hists = {}
for s in samples:
    hists[s] = {}
    for var in shower_vars:
        v = shower_vars[var]
        hists[s][var] = ROOT.TH1F(s + "_" + var, s, v["nbins"], v["xmin"], v["xmax"])

for s in samples:
    print(f"\nProcessing sample: {s}")
    target_pdg = samples[s]["pdg"]
    event_count = 0
    n_found = 0

    reader = IOIMPL.LCFactory.getInstance().createLCReader()
    reader.open(samples[s]["path"])

    for event in reader:

        event_count += 1

        try:
            pfo_coll = event.getCollection("PandoraPFOs")
        except:
            continue
        for pfo in pfo_coll:
            print(f"  PFO type: {pfo.getType()}, energy: {pfo.getEnergy():.2f}")
        # Only PFOs Pandora labeled with the target PDG
        matched_pfos = [pfo for pfo in pfo_coll if abs(pfo.getType()) == target_pdg]
        hists[s]["n_clus"].Fill(len(matched_pfos))

        if len(matched_pfos) < 1:
            continue

        n_found += 1

        # Pick the highest-momentum matched PFO
        best_pfo = None
        best_p   = 0.0

        for pfo in matched_pfos:
            trks = pfo.getTracks()
            if len(trks) < 1:
                continue
            p = plotHelper.getP(trks[0])
            if p > best_p:
                best_p   = p
                best_pfo = pfo

        if best_pfo is None:
            continue

        clus_list = best_pfo.getClusters()
        if len(clus_list) < 1:
            continue

        # Leading cluster by energy from the matched PFO
        lead_clus = max(clus_list, key=lambda c: c.getEnergy())
        clus_E = lead_clus.getEnergy()
        hists[s]["clus_E"].Fill(clus_E)

        # ECAL/HCAL split
        sub_E = lead_clus.getSubdetectorEnergies()
        if len(sub_E) > 1 and (sub_E[0] + sub_E[1]) > 0:
            hists[s]["ecal_frac"].Fill(sub_E[0] / (sub_E[0] + sub_E[1]))
            hists[s]["hcal_E"].Fill(sub_E[1])

        # Hit-level shower shape
        hits = lead_clus.getCalorimeterHits()

        if len(hits) > 0:

            hists[s]["n_hits"].Fill(len(hits))

            clus_tlv = plotHelper.getClusterTLV(lead_clus)
            axis = clus_tlv.Vect().Unit()

            min_r   = 1e9
            sum_E   = 0.0
            sum_Er  = 0.0
            sum_Ed2 = 0.0

            for hit in hits:

                pos   = hit.getPosition()
                hit_E = hit.getEnergy()
                r = math.sqrt(pos[0]**2 + pos[1]**2)

                if r < min_r:
                    min_r = r

                sum_E  += hit_E
                sum_Er += hit_E * r

                hit_vec = ROOT.TVector3(pos[0], pos[1], pos[2])
                d_perp = (hit_vec - (hit_vec.Dot(axis)) * axis).Mag()
                sum_Ed2 += hit_E * d_perp**2

            if sum_E > 0:
                hists[s]["start_r"].Fill(min_r)
                hists[s]["depth_r"].Fill(sum_Er / sum_E)
                hists[s]["width"].Fill(math.sqrt(sum_Ed2 / sum_E))

        # E/p for the best matched PFO
        if best_p > 0:
            e_calo = sum(c.getEnergy() for c in clus_list)
            hists[s]["E_over_p"].Fill(e_calo / best_p)

        if event_count == max_events:
            break

    reader.close()
    print(f"  Processed {event_count} events")
    print(f"  Events with a PDG {target_pdg} PFO: {n_found} ({100*n_found/event_count:.1f}%)")

# Normalize to unit area for shape comparison
for s in samples:
    for var in shower_vars:
        integral = hists[s][var].Integral()
        if integral > 0:
            hists[s][var].Scale(1.0 / integral)

# Make overlay plots
for var in shower_vars:
    h_map = {}
    for s in samples:
        h_map[s] = hists[s][var]
    plotHelper.plotHistograms(
        h_map,
        os.path.join(PLOT_DIR, f"shower_compare_pdgfilter_{var}.png"),
        xlabel=shower_vars[var]["label"],
        ylabel="Normalized Entries",
        atltext=["Muon Collider", "Simulation, no BIB", "#theta = #pi/2, p_{T} = 25 GeV"]
    )

print(f"\nDone. Plots saved to {PLOT_DIR}")
