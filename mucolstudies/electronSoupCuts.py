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
samples = {}
samples["Electrons"] = {"files": glob.glob("/scratch/jwatts/mucol/v11Container/reco/10kelectron0to50_reco.slcio"), "pdg": 11}
samples["Pions"]     = {"files": glob.glob("/scratch/jwatts/mucol/v2.11/reco_v2/pions_0_50/*.slcio"),              "pdg": 211}
max_events = 1000
dr_match   = 0.2
eta_max    = 2.4
def getClusterRMS(cluster):
    pos = cluster.getPosition()
    r = math.sqrt(pos[0]**2 + pos[1]**2 + pos[2]**2)
    if r == 0:
        return -1
    ax = [pos[0]/r, pos[1]/r, pos[2]/r]
    sum_e = 0.0
    sum_ed2 = 0.0
    for hit in cluster.getCalorimeterHits():
        hp = hit.getPosition()
        e = hit.getEnergy()
        dot = hp[0]*ax[0] + hp[1]*ax[1] + hp[2]*ax[2]
        d2 = hp[0]**2 + hp[1]**2 + hp[2]**2 - dot**2
        sum_e += e
        sum_ed2 += e*max(d2, 0)
    if sum_e == 0:
        return -1
    return math.sqrt(sum_ed2/sum_e)
def getSoupVariables(pfo):
    ecal_raw = 0.0
    hcal_raw = 0.0
    for cl in pfo.getClusters():
        sub = cl.getSubdetectorEnergies()
        if len(sub) > 0:
            ecal_raw += sub[0]
        if len(sub) > 1:
            hcal_raw += sub[1]
    tot_raw = ecal_raw + hcal_raw
    hcal_frac = hcal_raw/tot_raw if tot_raw > 0 else -1
    lead_cl = max(pfo.getClusters(), key=lambda c: c.getEnergy())
    rms = getClusterRMS(lead_cl)
    p = getP(pfo.getTracks()[0])
    eop = lead_cl.getEnergy()/p if p > 0 else -1
    return hcal_frac, rms, eop
profs = {}
for s in samples:
    profs[s] = {}
    profs[s]["hcal_frac"] = ROOT.TProfile(f"{s}_hcal_frac", "", 25, 0, 100)
    profs[s]["rms"]       = ROOT.TProfile(f"{s}_rms",       "", 25, 0, 100)
    profs[s]["eop"]       = ROOT.TProfile(f"{s}_eop",       "", 25, 0, 100)
reader = IOIMPL.LCFactory.getInstance().createLCReader()
for s in samples:
    files = sorted(samples[s]["files"])
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
                pfos = event.getCollection("PandoraPFOs")
            except:
                continue
            mcp_tlv = None
            for mcp in mcps:
                if mcp.getGeneratorStatus() != 1 or abs(mcp.getPDG()) != pdg: continue
                tlv = getTLV(mcp)
                if abs(tlv.Eta()) > eta_max: continue
                mcp_tlv = tlv
                break
            if mcp_tlv is None:
                continue
            # matched PFO: track + cluster, within dr_match, highest energy
            in_cone = []
            for pfo in pfos:
                if len(pfo.getTracks()) < 1 or len(pfo.getClusters()) < 1: continue
                if getTLV(pfo).DeltaR(mcp_tlv) < dr_match:
                    in_cone.append(pfo)
            if len(in_cone) == 0:
                continue
            best_pfo = max(in_cone, key=lambda x: x.getEnergy())
            hcal_frac, rms, eop = getSoupVariables(best_pfo)
            if hcal_frac >= 0: profs[s]["hcal_frac"].Fill(mcp_tlv.E(), hcal_frac)
            if rms >= 0:       profs[s]["rms"].Fill(mcp_tlv.E(), rms)
            if eop >= 0:       profs[s]["eop"].Fill(mcp_tlv.E(), eop)
        reader.close()
    print(f"\n{s}: used {n_events} events")

ROOT.gStyle.SetOptStat(0)

# overlay the TProfiles directly, no markers or error bars
labels = {"hcal_frac": "average HCAL fraction", "rms": "average shower RMS [mm]", "eop": "average E/p"}
for var in labels:
    c = ROOT.TCanvas("can", "can")
    keys = list(samples.keys())
    maxy = 1.3*max(profs[s][var].GetMaximum() for s in keys)
    for j, s in enumerate(keys):
        profs[s][var].SetLineColor(colors[j % len(colors)])
        profs[s][var].SetLineWidth(2)
        if j == 0:
            profs[s][var].SetTitle("")
            profs[s][var].GetXaxis().SetTitle("Truth E [GeV]")
            profs[s][var].GetYaxis().SetTitle(labels[var])
            profs[s][var].SetMinimum(0)
            profs[s][var].SetMaximum(maxy)
            profs[s][var].Draw("hist")
        else:
            profs[s][var].Draw("hist same")
    leg = ROOT.TLegend(0.68, 0.72, 0.93, 0.88)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.035)
    for s in keys:
        leg.AddEntry(profs[s][var], s, "l")
    leg.Draw()
    c.SaveAs(os.path.join(PLOT_DIR, f"AVG_{var}.png"))
    c.Close()

print(f"\nPlots saved to {PLOT_DIR}")
