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

samples = glob.glob("/scratch/jwatts/mucol/v11Container/reco/10kelectron0to50_reco.slcio")
#samples       = glob.glob("/scratch/jwatts/mucol/v2.11/reco/electronGun_pT_0_50/electronGun_pT_0_50_reco_0.slcio")
dr_match      = 0.2
eta_max       = 2.4
hcal_frac_max = 0.1     # raw HCAL / (raw ECAL + raw HCAL)
rms_max       = 30.0    # mm, energy-weighted transverse cluster RMS
eop_window    = 0.6     # |E/p - 1| < this

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

def passSoup(pfo):
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
    if hcal_frac < 0 or hcal_frac > hcal_frac_max:
        return False
    lead_cl = max(pfo.getClusters(), key=lambda c: c.getEnergy())
    rms = getClusterRMS(lead_cl)
    if rms < 0 or rms > rms_max:
        return False
    p = getP(pfo.getTracks()[0])
    eop = lead_cl.getEnergy()/p if p > 0 else -1
    if abs(eop - 1) > eop_window:
        return False
    return True

h_denom = ROOT.TH1F("mcp_el_eta",  ";#eta;Counts", 50, -2.5, 2.5)
h_soup  = ROOT.TH1F("soup_el_eta", ";#eta;Counts", 50, -2.5, 2.5)
h_pand  = ROOT.TH1F("pand_el_eta", ";#eta;Counts", 50, -2.5, 2.5)   # delete if soup only

reader = IOIMPL.LCFactory.getInstance().createLCReader()

files = sorted(samples)
print(f"found {len(files)} files")

n_events = 0

for f in files:
    reader.open(f)
    for event in reader:

        n_events += 1

        try:
            mcps = event.getCollection("MCParticle")
            pfos = event.getCollection("PandoraPFOs")
        except:
            continue

        mcp_electrons = []
        for mcp in mcps:
            if mcp.getGeneratorStatus() != 1 or abs(mcp.getPDG()) != 11: continue
            tlv = getTLV(mcp)
            if abs(tlv.Eta()) > eta_max: continue
            h_denom.Fill(tlv.Eta())
            mcp_electrons.append(tlv)

        good_pfos = []
        for pfo in pfos:
            if len(pfo.getTracks()) < 1 or len(pfo.getClusters()) < 1: continue
            good_pfos.append((getTLV(pfo), pfo))

        for mcp_el in mcp_electrons:

            in_cone = [(tlv, pfo) for tlv, pfo in good_pfos if mcp_el.DeltaR(tlv) < dr_match]
            if len(in_cone) == 0:
                continue
            best_tlv, best_pfo = max(in_cone, key=lambda x: x[1].getEnergy())

            if passSoup(best_pfo):
                h_soup.Fill(mcp_el.Eta())
            if abs(best_pfo.getType()) == 11:                            # delete if soup only
                h_pand.Fill(mcp_el.Eta())                                # delete if soup only

    reader.close()

print(f"\nProcessed {n_events} events")

print("\n--- Efficiency Results ---")
n_den = h_denom.GetEntries()
for name, h in [("Electron soup", h_soup), ("Pandora ID=11", h_pand)]:
    n_num = h.GetEntries()
    eff = n_num/n_den if n_den > 0 else 0
    print(f"{name:<15} | Eff: {eff:.4f} ({int(n_num)}/{int(n_den)})")

eff_map = {}
eff_map["Electron soup"] = ROOT.TEfficiency(h_soup, h_denom).CreateGraph()
eff_map["Pandora ID=11"] = ROOT.TEfficiency(h_pand, h_denom).CreateGraph()   # delete if soup only

plotEfficiencies(eff_map, os.path.join(PLOT_DIR, "EFFICIENCY_SOUP_ETA.svg"),
                 xlabel="#eta", ylabel="Efficiency")

print(f"\nPlots saved to {PLOT_DIR}")
