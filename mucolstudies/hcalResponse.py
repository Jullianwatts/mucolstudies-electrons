import math
import glob
import ROOT
import pyLCIO
from pyLCIO import EVENT, IOIMPL
import os

exec(open("./plotHelper.py").read())

ROOT.gROOT.SetBatch(True)

PLOT_DIR = "/scratch/jwatts/mucol/mucolstudies/plots2026"
os.makedirs(PLOT_DIR, exist_ok=True)

RECO_DIR = "/scratch/jwatts/mucol/v2.11/reco_v2"

samples = {}
#samples["electrons_0_50"]             = {"glob": f"{RECO_DIR}/electrons_0_50/*.slcio",             "pdg": 11}
#samples["photonGun_transitionRegion"] = {"glob": f"{RECO_DIR}/photonGun_transitionRegion/*.slcio", "pdg": 22}
samples["pions_0_50"]                 = {"glob": f"{RECO_DIR}/pions_0_50/*.slcio",                 "pdg": 211}
samples["pions_50_250"]               = {"glob": f"{RECO_DIR}/pions_50_250/*.slcio",               "pdg": 211}
samples["pions_250_1000"]             = {"glob": f"{RECO_DIR}/pions_250_1000/*.slcio",             "pdg": 211}
samples["pions_1000_5000"]            = {"glob": f"{RECO_DIR}/pions_1000_5000/*.slcio",            "pdg": 211}

max_events = 20000  # per sample
eta_max = 2.4

# hadronic scale factors applied to raw subdetector energies
# (getSubdetectorEnergies is raw regardless of the steering constants,
# which only affect cluster.getEnergy)
ecal_to_had = 1.38
hcal_to_had = 1.25

energy_bins = [0, 50, 100, 150, 200, 250, 500, 1000, 2500, 5000]

h2 = {}
prof = {}
h_res = {}
for s in samples:
    h2[s] = {}
    prof[s] = {}
    h_res[s] = {}
    for i in range(len(energy_bins) - 1):
        name = f"{energy_bins[i]}_{energy_bins[i+1]}"
        h2[s][i] = ROOT.TH2F(f"h2_{s}_{name}", f"{s} {name};HCAL Fraction;Response", 20, 0, 1, 40, -1, 3)
        prof[s][name] = ROOT.TProfile(f"prof_{s}_{name}", "", 20, 0, 1)
        h_res[s][name] = ROOT.TH1F(f"res_{s}_{name}", "", 80, -2, 2)

reader = IOIMPL.LCFactory.getInstance().createLCReader()

for s in samples:

    files = sorted(glob.glob(samples[s]["glob"]))
    pdg = samples[s]["pdg"]
    print(f"{s}: found {len(files)} files")

    n_events = 0

    for f in files:
        if n_events >= max_events:
            break
        reader.open(f)
        for event in reader:

            if n_events >= max_events:
                break
            n_events += 1

            try:
                mcps = event.getCollection("MCParticle")
                clusters = event.getCollection("PandoraClusters")
            except:
                continue

            # truth gun particle
            mcp_tlv = None
            for mcp in mcps:
                if mcp.getGeneratorStatus() != 1 or abs(mcp.getPDG()) != pdg: continue
                tlv = getTLV(mcp)
                if abs(tlv.Eta()) > eta_max: continue
                mcp_tlv = tlv
                break
            if mcp_tlv is None:
                continue

            truth_e = mcp_tlv.E()
            ebin = -1
            for i in range(len(energy_bins) - 1):
                if energy_bins[i] <= truth_e < energy_bins[i+1]:
                    ebin = i
                    break
            if ebin < 0:
                continue

            # raw subdetector energies summed over all clusters in the event
            ecal_raw = 0.0
            hcal_raw = 0.0
            for cl in clusters:
                sub = cl.getSubdetectorEnergies()
                if len(sub) > 0:
                    ecal_raw += sub[0]
                if len(sub) > 1:
                    hcal_raw += sub[1]
            tot_raw = ecal_raw + hcal_raw
            if tot_raw <= 0:
                continue

            hcal_frac = hcal_raw/tot_raw
            measured = ecal_raw*ecal_to_had + hcal_raw*hcal_to_had
            response = measured/truth_e

            h2[s][ebin].Fill(hcal_frac, response)
            name = f"{energy_bins[ebin]}_{energy_bins[ebin+1]}"
            prof[s][name].Fill(hcal_frac, response)
            h_res[s][name].Fill(response - 1)

        reader.close()

    print(f"{s}: used {n_events} events")

def plotOverlay(h_map, save_name, xlabel, ylabel):
    # plain line overlay, no markers or error bars
    can = ROOT.TCanvas("can", "can", 800, 600)
    can.SetLeftMargin(0.12)
    can.SetBottomMargin(0.12)
    can.SetRightMargin(0.05)
    keys = list(h_map.keys())
    maxy = 1.3*max(h_map[k].GetMaximum() for k in keys)
    for j, k in enumerate(keys):
        h_map[k].SetLineColor(colors[j % len(colors)])
        h_map[k].SetLineWidth(2)
        if j == 0:
            h_map[k].SetTitle("")
            h_map[k].GetXaxis().SetTitle(xlabel)
            h_map[k].GetYaxis().SetTitle(ylabel)
            h_map[k].SetMinimum(0)
            h_map[k].SetMaximum(maxy)
            h_map[k].Draw("hist")
        else:
            h_map[k].Draw("hist same")
    leg = ROOT.TLegend(0.68, 0.65, 0.93, 0.88)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.035)
    for k in keys:
        leg.AddEntry(h_map[k], k, "l")
    leg.Draw()
    can.SaveAs(save_name)
    can.Close()
    return

# 2D response vs HCAL fraction with red profile + linear fit, one canvas per sample per energy bin
ROOT.gStyle.SetOptStat(0)
print(f"\n{'Sample':<28} | {'E bin':<10} | {'fit intercept':<13} | {'fit slope':<10}")
for s in samples:
    for i in range(len(energy_bins) - 1):
        if h2[s][i].GetEntries() == 0:
            continue
        name = f"{energy_bins[i]}_{energy_bins[i+1]}"
        can = ROOT.TCanvas("can", "can", 800, 600)
        can.SetLeftMargin(0.12)
        can.SetBottomMargin(0.12)
        can.SetRightMargin(0.14)
        h2[s][i].Draw("COLZ")
        p = h2[s][i].ProfileX(f"pfx_{s}_{name}")
        p.SetLineColor(ROOT.kRed)
        p.SetLineWidth(2)
        fit = ROOT.TF1(f"fit_{s}_{name}", "pol1", 0, 1)
        p.Fit(fit, "QN")
        p.Draw("hist same")
        fit.SetLineColor(ROOT.kRed)
        fit.SetLineStyle(2)
        fit.Draw("same")
        print(f"{s:<28} | {name:<10} | {fit.GetParameter(0):<13.3f} | {fit.GetParameter(1):<10.3f}")
        can.SaveAs(os.path.join(PLOT_DIR, f"RESPONSE_VS_HCALFRAC_{s}_{name}.png"))
        can.Close()

# per sample: overlay of average response vs HCAL fraction for its energy bins
for s in samples:
    p_map = {}
    for name, p in prof[s].items():
        if p.GetEntries() == 0:
            continue
        h = ROOT.TH1F(f"line_{s}_{name}", "", p.GetNbinsX(), 0, 1)
        for b in range(1, p.GetNbinsX() + 1):
            h.SetBinContent(b, p.GetBinContent(b))
        p_map[name] = h
    if len(p_map) == 0:
        continue
    plotOverlay(p_map, os.path.join(PLOT_DIR, f"RESPONSE_VS_HCALFRAC_PROFILE_{s}.png"),
                "HCAL Fraction", "Response")

# per sample: 1D (E_meas - E_true)/E_true distributions, overlaid per energy bin
for s in samples:
    r_map = {name: h for name, h in h_res[s].items() if h.GetEntries() > 0}
    if len(r_map) == 0:
        continue
    plotOverlay(r_map, os.path.join(PLOT_DIR, f"ERES_{s}.png"),
                "(E_{meas} - E_{true})/E_{true}", "Counts")

print(f"\nPlots saved to {PLOT_DIR}")
