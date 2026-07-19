import math
import glob
import ROOT
import pyLCIO
from pyLCIO import EVENT, IOIMPL, UTIL
import os
exec(open("./plotHelper.py").read())
ROOT.gROOT.SetBatch(True)
PLOT_DIR = "/scratch/jwatts/mucol/mucolstudies/plots2026"
os.makedirs(PLOT_DIR, exist_ok=True)
samples = {}
samples["Electrons"] = {"files": glob.glob("/scratch/jwatts/mucol/v11Container/reco/10kelectron0to50_reco.slcio"), "pdg": 11}
samples["Pions"]     = {"files": glob.glob("/scratch/jwatts/mucol/v2.11/reco_v2/pions_0_50/*.slcio"),              "pdg": 211}
max_events = 10000
dr_match   = 0.2
eta_max    = 2.4

# ---- Shower profile constants (from LCShowerProfilePlugin.cc) ----
LONG_PROFILE_BIN_WIDTH = 0.5        # X0
LONG_PROFILE_N_BINS = 100
LONG_PROFILE_MIN_COS_ANGLE = 0.3
LONG_PROFILE_CRITICAL_ENERGY = 0.08 # GeV
LONG_PROFILE_PARAMETER_0 = 1.25
LONG_PROFILE_PARAMETER_1 = 0.5
LONG_PROFILE_MAX_DIFFERENCE = 0.1

# ---- From steer_reco.py (DDMarlinPandora block) ----
ECAL_TO_EM_GEV = 1.02373335516      # ECalToEMGeVCalibration
ECAL_TO_MIP = 181.818               # ECalToMipCalibration
ECAL_MIP_THRESHOLD = 0.5            # ECalMipThreshold

# ---- From MAIA_v0 ECalBarrel_o2_v01_02.xml: 50 layers, 2.20mm W ----
X0_PER_LAYER = 2.20/3.504 + 0.012   # ~0.64 X0 per layer
INNER_LAYER_OFFSET = 1

ecal_collections = ["EcalBarrelCollectionSel", "EcalEndcapCollectionSel"]

def get_ecal_hit_info(event):
    """Map (cellID0, cellID1) -> (layer, is_barrel) for ECal hits Pandora actually used."""
    info = {}
    for cname in ecal_collections:
        try:
            coll = event.getCollection(cname)
        except:
            continue
        encoding = coll.getParameters().getStringVal(EVENT.LCIO.CellIDEncoding)
        decoder = UTIL.BitField64(encoding)
        is_barrel = "Barrel" in cname
        for hit in coll:
            cellid = ((hit.getCellID1() & 0xffffffff) << 32) | (hit.getCellID0() & 0xffffffff)
            decoder.setValue(cellid)
            layer = decoder["layer"].value()
            info[(hit.getCellID0(), hit.getCellID1())] = (layer, is_barrel)
    return info

def compute_shower_profile(cluster, hit_info):
    """Mirror of LCShowerProfilePlugin::CalculateLongitudinalProfile.
    Returns (e_ecal, s_start, s_disc) or None."""
    kept = []
    sume, sumx, sumy, sumz = 0.0, 0.0, 0.0, 0.0
    for hit in cluster.getCalorimeterHits():
        key = (hit.getCellID0(), hit.getCellID1())
        if key not in hit_info:
            continue
        e = hit.getEnergy()
        if e * ECAL_TO_MIP < ECAL_MIP_THRESHOLD:
            continue
        layer, is_barrel = hit_info[key]
        hp = hit.getPosition()
        kept.append((layer, is_barrel, e, (hp[0], hp[1], hp[2])))
        sume += e
        sumx += e*hp[0]
        sumy += e*hp[1]
        sumz += e*hp[2]
    if sume <= 0 or not kept:
        return None
    norm = math.sqrt(sumx**2 + sumy**2 + sumz**2)
    if norm == 0:
        return None
    cdir = (sumx/norm, sumy/norm, sumz/norm)
    hits_by_layer = {}
    for layer, is_barrel, e, hp in kept:
        if is_barrel:
            r = math.sqrt(hp[0]**2 + hp[1]**2)
            if r == 0:
                continue
            normal = (hp[0]/r, hp[1]/r, 0.0)
        else:
            normal = (0.0, 0.0, 1.0 if hp[2] > 0 else -1.0)
        cos_angle = abs(normal[0]*cdir[0] + normal[1]*cdir[1] + normal[2]*cdir[2])
        cos_angle = max(cos_angle, LONG_PROFILE_MIN_COS_ANGLE)
        hits_by_layer.setdefault(layer, []).append((e * ECAL_TO_EM_GEV, X0_PER_LAYER / cos_angle))
    if not hits_by_layer:
        return None
    profile = [0.0] * LONG_PROFILE_N_BINS
    e_ecal = 0.0
    n_rad = 0.0
    n_rad_last = 0.0
    inner_layer = min(hits_by_layer)
    for layer in range(inner_layer, max(hits_by_layer) + 1):
        if layer not in hits_by_layer:
            n_rad += n_rad_last
            continue
        layer_hits = hits_by_layer[layer]
        e_layer = sum(h[0] for h in layer_hits)
        x0_layer = sum(h[1] for h in layer_hits) / len(layer_hits)
        e_ecal += e_layer
        n_rad_last = x0_layer
        n_rad += x0_layer
        if layer == inner_layer:
            n_rad *= float(inner_layer + INNER_LAYER_OFFSET)
        end_pos = n_rad / LONG_PROFILE_BIN_WIDTH
        end_bin = min(int(end_pos), LONG_PROFILE_N_BINS - 1)
        delta_pos = x0_layer / LONG_PROFILE_BIN_WIDTH
        start_pos = end_pos - delta_pos
        start_bin = int(start_pos)
        for ibin in range(start_bin, end_bin + 1):
            if ibin >= LONG_PROFILE_N_BINS:
                break
            delta = 1.0
            if ibin == start_bin:
                delta -= start_pos - start_bin
            elif ibin == end_bin:
                delta -= 1.0 - end_pos + end_bin
            profile[ibin] += e_layer * (delta / delta_pos)
    profile_end_bin = min(int(n_rad / LONG_PROFILE_BIN_WIDTH), LONG_PROFILE_N_BINS)
    if profile_end_bin == 0 or e_ecal <= 0:
        return None
    a = LONG_PROFILE_PARAMETER_0 + LONG_PROFILE_PARAMETER_1 * math.log(e_ecal / LONG_PROFILE_CRITICAL_ENERGY)
    gamma_a = math.exp(math.lgamma(a))
    expected = []
    t = 0.0
    for ibin in range(LONG_PROFILE_N_BINS):
        t += LONG_PROFILE_BIN_WIDTH
        expected.append(e_ecal / 2.0 * (t/2.0)**(a - 1.0) * math.exp(-t/2.0) * LONG_PROFILE_BIN_WIDTH / gamma_a)
    min_diff = float("inf")
    best_offset = 0
    for offset in range(profile_end_bin):
        diff = 0.0
        for ibin in range(profile_end_bin):
            if ibin < offset:
                diff += profile[ibin]
            else:
                diff += abs(expected[ibin - offset] - profile[ibin])
        if diff < min_diff:
            min_diff = diff
            best_offset = offset
        if diff - min_diff > LONG_PROFILE_MAX_DIFFERENCE:
            break
    return e_ecal, best_offset * LONG_PROFILE_BIN_WIDTH, min_diff / e_ecal

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

def getSoupVariables(pfo, hit_info):
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
    result = compute_shower_profile(lead_cl, hit_info)
    if result is None:
        sstart, sdisc = -1, -1
    else:
        e_ecal, sstart, sdisc = result
    return hcal_frac, rms, eop, sstart, sdisc

profs = {}
for s in samples:
    profs[s] = {}
    profs[s]["hcal_frac"] = ROOT.TProfile(f"{s}_hcal_frac", "", 25, 0, 100)
    profs[s]["rms"]       = ROOT.TProfile(f"{s}_rms",       "", 25, 0, 100)
    profs[s]["eop"]       = ROOT.TProfile(f"{s}_eop",       "", 25, 0, 100)
    profs[s]["sstart"]    = ROOT.TProfile(f"{s}_sstart",    "", 25, 0, 100)
    profs[s]["sdisc"]     = ROOT.TProfile(f"{s}_sdisc",     "", 25, 0, 100)
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
            hit_info = get_ecal_hit_info(event)
            # matched PFO: track + cluster, within dr_match, highest energy
            in_cone = []
            for pfo in pfos:
                if len(pfo.getTracks()) < 1 or len(pfo.getClusters()) < 1: continue
                if getTLV(pfo).DeltaR(mcp_tlv) < dr_match:
                    in_cone.append(pfo)
            if len(in_cone) == 0:
                continue
            best_pfo = max(in_cone, key=lambda x: x.getEnergy())
            hcal_frac, rms, eop, sstart, sdisc = getSoupVariables(best_pfo, hit_info)
            if hcal_frac >= 0: profs[s]["hcal_frac"].Fill(mcp_tlv.E(), hcal_frac)
            if rms >= 0:       profs[s]["rms"].Fill(mcp_tlv.E(), rms)
            if eop >= 0:       profs[s]["eop"].Fill(mcp_tlv.E(), eop)
            if sstart >= 0:    profs[s]["sstart"].Fill(mcp_tlv.E(), sstart)
            if sdisc >= 0:     profs[s]["sdisc"].Fill(mcp_tlv.E(), sdisc)
        reader.close()
    print(f"\n{s}: used {n_events} events")
ROOT.gStyle.SetOptStat(0)
# overlay the TProfiles directly, no markers or error bars
labels = {"hcal_frac": "average HCAL fraction", "rms": "average shower RMS [mm]", "eop": "average E/p",
          "sstart": "average shower profile start [X0]", "sdisc": "average shower profile discrepancy"}
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
