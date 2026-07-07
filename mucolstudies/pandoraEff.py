import math
import glob
import ROOT
import pyLCIO
import os

exec(open("./plotHelper.py").read())

ROOT.gROOT.SetBatch()

PLOT_DIR = "/scratch/jwatts/mucol/mucolstudies/plots2026"
os.makedirs(PLOT_DIR, exist_ok=True)

samples = glob.glob("/scratch/jwatts/mucol/v11Container/reco/10kelectron0to50_reco.slcio")

files = {"10kelectron0to50": samples}

hists = {}
for s in files:
    hists[s] = {}
    for obj in ["mcp_el", "pfo_el_match"]:
        hists[s][obj + "_eta"] = ROOT.TH1F(f"{s}_{obj}_eta", ";#eta;Counts", 50, -2.5, 2.5)

reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()

for s in files:
    for f in files[s]:
        reader.open(f)
        for event in reader:

            mcps = event.getCollection("MCParticle")
            pfos = event.getCollection("PandoraPFOs")

            mcp_electrons = []
            for mcp in mcps:
                if mcp.getGeneratorStatus() != 1 or abs(mcp.getPDG()) != 11: continue
                tlv = getTLV(mcp)
                eta = tlv.Eta()
                if abs(eta) > 2.4: continue
                hists[s]["mcp_el_eta"].Fill(eta)
                mcp_electrons.append(tlv)

            pfo_electrons = [getTLV(p) for p in pfos if abs(p.getType()) == 11]

            for mcp_el in mcp_electrons:
                m_eta = mcp_el.Eta()

                # pfos in cone -> closest dR
                pfo_drs_in_cone = [(mcp_el.DeltaR(p), p) for p in pfo_electrons if mcp_el.DeltaR(p) < 0.2]
                best_pfo = min(pfo_drs_in_cone, key=lambda x: x[0], default=(999, None))

                if best_pfo[0] < 0.2:
                    hists[s]["pfo_el_match_eta"].Fill(m_eta)

        reader.close()

print("\n--- Efficiency Results ---")
for s in hists:

    eff_map = {}
    denom = hists[s]["mcp_el_eta"]
    num = hists[s]["pfo_el_match_eta"]

    n_num = num.GetEntries()
    n_den = denom.GetEntries()
    efficiency = n_num / n_den if n_den > 0 else 0
    print(f"Object: {'PandoraPFOs (PDGID=11)':<25} | Eff: {efficiency:.4f} ({int(n_num)}/{int(n_den)})")

    teff = ROOT.TEfficiency(num, denom)
    graph = teff.CreateGraph()
    eff_map["PandoraPFOs (PDGID=11)"] = graph

    plotEfficiencies(eff_map, os.path.join(PLOT_DIR, "EFFICIENCY_PFO_ETA.svg"),
                     xlabel="#eta", ylabel="Efficiency")

print(f"\nPlots saved to {PLOT_DIR}")
