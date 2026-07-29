import math
import glob
import ROOT
import pyLCIO
from pyLCIO import EVENT, IOIMPL, UTIL
import os
from array import array
exec(open("./plotHelper.py").read())
ROOT.gROOT.SetBatch(True)
PLOT_DIR = "/scratch/jwatts/mucol/mucolstudies/plots2026"
os.makedirs(PLOT_DIR, exist_ok=True)
samples = glob.glob("/scratch/jwatts/mucol/v11Container/reco/10kelectron0to50_reco.slcio")
#samples = glob.glob("/scratch/jwatts/mucol/v2.11/reco/electronGun_pT_0_50/electronGun_pT_0_50_reco_Final.slcio")
pion_samples = glob.glob("/scratch/jwatts/mucol/v2.11/reco_v2/pions_0_50/*.slcio")
max_events    = 10000
eta_max       = 2.4
# ---- Electron ID cuts (applied cumulatively, in this order) ----
hcal_frac_max = 0.001   # raw HCAL / (raw ECAL + raw HCAL)
rms_max       = 26.0    # mm, energy-weighted transverse cluster RMS
sstart_max    = 6.0     # X0, shower profile start
# ---- Threshold scan values (printed, to verify cut choices) ----
scan_values = {
    "hcal_frac": [0.0005, 0.001, 0.005, 0.01, 0.05, 0.1],
    "rms":       [24.0, 26.0, 28.0, 30.0],
    "sstart":    [4.5, 6.0, 8.0],
}
# ---- Truth energy ranges for distribution plots: label -> (min, max) ----
e_ranges = {
    "E_0_15":  (0,  15),
    "E_15_35": (15, 35),
    "E_35_50": (35, 50),
}
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
def getCutVariables(pfo, hit_info):
    """Returns (hcal_frac, rms, sstart). -1 values where computation not possible."""
    if len(pfo.getClusters()) < 1:
        return -1, -1, -1
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
    result = compute_shower_profile(lead_cl, hit_info)
    if result is None:
        sstart = -1
    else:
        e_ecal, sstart, s_disc = result
    return hcal_frac, rms, sstart
def getBestPFO(pfos):
    """Highest-pT charged PFO (>= 1 track, |eta| < eta_max) in the event, or None."""
    best_pfo, best_pt = None, -1
    for pfo in pfos:
        if len(pfo.getTracks()) < 1:
            continue
        tlv = getTLV(pfo)
        if abs(tlv.Eta()) > eta_max:
            continue
        pt = tlv.Perp()
        if pt > best_pt:
            best_pfo = pfo
            best_pt = pt
    return best_pfo
def passStage(stage, hcal_frac, rms, sstart):
    """Cumulative stages: 1 = HCAL, 2 = HCAL+RMS, 3 = HCAL+RMS+sstart."""
    if not (0 <= hcal_frac <= hcal_frac_max):
        return False
    if stage >= 2 and not (0 <= rms <= rms_max):
        return False
    if stage >= 3 and not (0 <= sstart <= sstart_max):
        return False
    return True
def passSingle(cut, hcal_frac, rms, sstart):
    """Each cut applied entirely on its own."""
    if cut == "hcal_frac": return (0 <= hcal_frac <= hcal_frac_max)
    if cut == "rms":       return (0 <= rms <= rms_max)
    if cut == "sstart":    return (0 <= sstart <= sstart_max)
    return False
# ---- Histograms: denominator + cumulative stages + individual cuts ----
stage_names  = {1: "HCAL frac", 2: "+ shower RMS", 3: "+ shower start"}
single_names = {"hcal_frac": "HCAL frac only", "rms": "shower RMS only", "sstart": "shower start only"}
h_denom = ROOT.TH1F("mcp_el_eta", ";#eta;Counts", 50, -2.5, 2.5)
h_stage = {}
for stage in stage_names:
    h_stage[stage] = ROOT.TH1F(f"stage{stage}_el_eta", ";#eta;Counts", 50, -2.5, 2.5)
h_single = {}
for cut in single_names:
    h_single[cut] = ROOT.TH1F(f"single_{cut}_el_eta", ";#eta;Counts", 50, -2.5, 2.5)
# ---- Normalized 1D distributions: variable, binning, cut value, x label, log axes ----
# hcal_frac uses log-spaced x bins from 1e-4 to 1; values below 1e-4 (incl. exactly 0) clamped into first bin
dist_vars = {
    "hcal_frac": {"nbins": 40, "xmin": 1e-4, "xmax": 1,   "cut": hcal_frac_max, "label": "HCAL fraction",             "logy": True,  "logx": True},
    "rms":       {"nbins": 50, "xmin": 0,    "xmax": 100, "cut": rms_max,       "label": "shower RMS [mm]",           "logy": False, "logx": False},
    "sstart":    {"nbins": 40, "xmin": 0,    "xmax": 20,  "cut": sstart_max,    "label": "shower profile start [X0]", "logy": True,  "logx": False},
}
def makeDistHist(name, v):
    if v["logx"]:
        edges = array('d', [10**(math.log10(v["xmin"]) + i*(math.log10(v["xmax"])-math.log10(v["xmin"]))/v["nbins"]) for i in range(v["nbins"]+1)])
        return ROOT.TH1F(name, "", v["nbins"], edges)
    return ROOT.TH1F(name, "", v["nbins"], v["xmin"], v["xmax"])
h_dist = {}
for erange in e_ranges:
    h_dist[erange] = {}
    for ptype in ["Electrons", "Pions"]:
        h_dist[erange][ptype] = {}
        for var in dist_vars:
            h_dist[erange][ptype][var] = makeDistHist(f"{erange}_{ptype}_{var}_dist", dist_vars[var])
def getERangeLabel(energy):
    for label, (emin, emax) in e_ranges.items():
        if emin <= energy < emax:
            return label
    return None
def fillDistHists(erange, ptype, hcal_frac, rms, sstart):
    if erange is None: return
    if hcal_frac >= 0:
        h_dist[erange][ptype]["hcal_frac"].Fill(max(hcal_frac, 1.2e-4))   # clamp 0 into first log bin
    if rms >= 0:       h_dist[erange][ptype]["rms"].Fill(rms)
    if sstart >= 0:    h_dist[erange][ptype]["sstart"].Fill(sstart)
# ---- Raw variable values, kept for threshold scans and fake rates ----
el_vals = []    # (hcal_frac, rms, sstart) per truth electron with a charged PFO
pi_vals = []    # (hcal_frac, rms, sstart) per truth pion with a charged PFO
reader = IOIMPL.LCFactory.getInstance().createLCReader()
files = sorted(samples)
print(f"found {len(files)} electron files")
n_events = 0
n_el_truth = 0
n_el_nopfo = 0
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
        hit_info = get_ecal_hit_info(event)
        mcp_electrons = []
        for mcp in mcps:
            if mcp.getGeneratorStatus() != 1 or abs(mcp.getPDG()) != 11: continue
            tlv = getTLV(mcp)
            if abs(tlv.Eta()) > eta_max: continue
            h_denom.Fill(tlv.Eta())
            mcp_electrons.append(tlv)
        best_pfo = getBestPFO(pfos)
        for mcp_el in mcp_electrons:
            n_el_truth += 1
            if best_pfo is None:
                n_el_nopfo += 1
                continue
            hcal_frac, rms, sstart = getCutVariables(best_pfo, hit_info)
            el_vals.append((hcal_frac, rms, sstart))
            fillDistHists(getERangeLabel(mcp_el.E()), "Electrons", hcal_frac, rms, sstart)
            for stage in stage_names:
                if passStage(stage, hcal_frac, rms, sstart):
                    h_stage[stage].Fill(mcp_el.Eta())
            for cut in single_names:
                if passSingle(cut, hcal_frac, rms, sstart):
                    h_single[cut].Fill(mcp_el.Eta())
    reader.close()
# ---- Pion loop: same truth-blind selection, for scans, fake rates, distributions ----
pion_files = sorted(pion_samples)
print(f"found {len(pion_files)} pion files")
n_pion_events = 0
for f in pion_files:
    if n_pion_events >= max_events:
        break
    reader.open(f)
    for event in reader:
        if n_pion_events >= max_events:
            break
        n_pion_events += 1
        try:
            mcps = event.getCollection("MCParticle")
            pfos = event.getCollection("PandoraPFOs")
        except:
            continue
        hit_info = get_ecal_hit_info(event)
        mcp_pions = []
        for mcp in mcps:
            if mcp.getGeneratorStatus() != 1 or abs(mcp.getPDG()) != 211: continue
            tlv = getTLV(mcp)
            if abs(tlv.Eta()) > eta_max: continue
            mcp_pions.append(tlv)
        if len(mcp_pions) == 0:
            continue
        best_pfo = getBestPFO(pfos)
        if best_pfo is None:
            continue
        hcal_frac, rms, sstart = getCutVariables(best_pfo, hit_info)
        for mcp_pi in mcp_pions:
            pi_vals.append((hcal_frac, rms, sstart))
            fillDistHists(getERangeLabel(mcp_pi.E()), "Pions", hcal_frac, rms, sstart)
    reader.close()
print(f"\nProcessed {n_events} electron events, {n_pion_events} pion events")
print(f"Truth electrons: {n_el_truth} | no charged PFO in event: {n_el_nopfo}")
print(f"Pions with charged PFO: {len(pi_vals)}")
# ---- Threshold scans: fraction of electrons / pions passing each candidate cut ----
def scanFrac(vals, idx, cut):
    if len(vals) == 0:
        return 0
    return sum(1 for v in vals if 0 <= v[idx] <= cut) / len(vals)
var_idx = {"hcal_frac": 0, "rms": 1, "sstart": 2}
print("\n--- Threshold Scans (fraction passing, electrons | pions) ---")
for var in scan_values:
    print(f"\n-- {var} --")
    for cut in scan_values[var]:
        f_el = scanFrac(el_vals, var_idx[var], cut)
        f_pi = scanFrac(pi_vals, var_idx[var], cut)
        print(f"  cut <= {cut:<8} | el: {f_el:.4f} | pi: {f_pi:.4f}")
# ---- Cumulative stage results: electron efficiency and pion fake rate ----
print("\n--- Cumulative Cut Results ---")
n_den = h_denom.GetEntries()
print(f"{'Stage':<18} | {'El eff':<10} | {'Pi fake rate':<12}")
for stage in stage_names:
    n_num = h_stage[stage].GetEntries()
    eff = n_num/n_den if n_den > 0 else 0
    n_pi_pass = sum(1 for v in pi_vals if passStage(stage, v[0], v[1], v[2]))
    fake = n_pi_pass/len(pi_vals) if len(pi_vals) > 0 else 0
    print(f"{stage_names[stage]:<18} | {eff:<10.4f} | {fake:<12.4f}")
print("\n--- Individual Cut Results (each cut alone) ---")
print(f"{'Cut':<18} | {'El eff':<10} | {'Pi fake rate':<12}")
for cut in single_names:
    n_num = h_single[cut].GetEntries()
    eff = n_num/n_den if n_den > 0 else 0
    n_pi_pass = sum(1 for v in pi_vals if passSingle(cut, v[0], v[1], v[2]))
    fake = n_pi_pass/len(pi_vals) if len(pi_vals) > 0 else 0
    print(f"{single_names[cut]:<18} | {eff:<10.4f} | {fake:<12.4f}")
# ---- Efficiency plots: cumulative overlay + one plot per individual cut ----
eff_map = {}
for stage in stage_names:
    eff_map[stage_names[stage]] = ROOT.TEfficiency(h_stage[stage], h_denom).CreateGraph()
plotEfficiencies(eff_map, os.path.join(PLOT_DIR, "EFFICIENCY_CUMULATIVE_ETA.png"),
                 xlabel="#eta", ylabel="Efficiency")
for cut in single_names:
    single_map = {}
    single_map[single_names[cut]] = ROOT.TEfficiency(h_single[cut], h_denom).CreateGraph()
    plotEfficiencies(single_map, os.path.join(PLOT_DIR, f"EFFICIENCY_{cut.upper()}_ONLY_ETA.png"),
                     xlabel="#eta", ylabel="Efficiency")
# ---- Normalized 1D distribution plots with cut lines, one per energy range ----
ROOT.gStyle.SetOptStat(0)
for var in dist_vars:
    v = dist_vars[var]
    for erange, (emin, emax) in e_ranges.items():
        c = ROOT.TCanvas("can", "can", 800, 600)
        c.SetLeftMargin(0.12)
        c.SetBottomMargin(0.12)
        c.SetRightMargin(0.05)
        if v["logy"]:
            c.SetLogy(1)
        if v["logx"]:
            c.SetLogx(1)
        for j, ptype in enumerate(["Electrons", "Pions"]):
            h = h_dist[erange][ptype][var]
            if h.Integral() > 0:
                h.Scale(1.0/h.Integral())
            h.SetLineColor(colors[j % len(colors)])
            h.SetLineWidth(2)
        if v["logy"]:
            miny = 1e-4
            maxy = 10*max(h_dist[erange][p][var].GetMaximum() for p in ["Electrons", "Pions"])
        else:
            miny = 0
            maxy = 1.3*max(h_dist[erange][p][var].GetMaximum() for p in ["Electrons", "Pions"])
        if maxy <= miny:
            maxy = 1
        for j, ptype in enumerate(["Electrons", "Pions"]):
            h = h_dist[erange][ptype][var]
            if j == 0:
                h.SetTitle("")
                h.GetXaxis().SetTitle(v["label"])
                h.GetYaxis().SetTitle("Normalized Entries")
                h.GetXaxis().SetTitleSize(0.045)
                h.GetYaxis().SetTitleSize(0.045)
                h.GetXaxis().SetLabelSize(0.04)
                h.GetYaxis().SetLabelSize(0.04)
                h.SetMinimum(miny)
                h.SetMaximum(maxy)
                h.Draw("hist")
            else:
                h.Draw("hist same")
        cut_line = ROOT.TLine(v["cut"], miny, v["cut"], maxy)
        cut_line.SetLineStyle(2)
        cut_line.SetLineWidth(2)
        cut_line.SetLineColor(ROOT.kGray+2)
        cut_line.Draw("same")
        leg = ROOT.TLegend(0.68, 0.72, 0.93, 0.88)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetTextSize(0.035)
        for ptype in ["Electrons", "Pions"]:
            leg.AddEntry(h_dist[erange][ptype][var], ptype, "l")
        leg.AddEntry(cut_line, "cut value", "l")
        leg.Draw()
        text = ROOT.TLatex()
        text.SetNDC()
        text.SetTextSize(0.04)
        text.SetTextFont(62)
        text.DrawLatex(0.15, 0.88, "Muon Collider")
        text.SetTextFont(42)
        text.DrawLatex(0.15, 0.83, "Simulation, no BIB")
        text.DrawLatex(0.15, 0.78, f"Truth E: {emin}-{emax} GeV")
        c.SaveAs(os.path.join(PLOT_DIR, f"DIST_{var}_{erange}.png"))
        c.Close()
print(f"\nPlots saved to {PLOT_DIR}")
