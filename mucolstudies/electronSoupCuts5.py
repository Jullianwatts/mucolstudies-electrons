#makes cut plots and plots efficiencies that build on eachother
#efficiency denominator = leading charged PFO per event (track + cluster), |eta| < eta_max
#cut studies are iterative: each variable is studied only on PFOs passing the previous cuts
#good script
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
#cuts are applied in this order
hcal_frac_max = 0.001   # raw HCAL / (raw ECAL + raw HCAL)
rms_max       = 28.0    # mm, energy-weighted transverse cluster RMS
sstart_max    = 7     # X0, shower profile start
scan_values = {
    "hcal_frac": [0.0005, 0.001, 0.005, 0.01, 0.05, 0.1],
    "rms":       [24.0, 26.0, 28.0, 30.0],
    "sstart":    [4.5, 6.0, 6.5, 7.0, 8.0],
}
#testing diff shower start layer values
sstart_variants = [4.5, 5.0, 5.5, 6.0, 6.5, 7.0]
#Truth energy ranges for distribution plots
e_ranges = {
    "E_0_15":  (0,  15),
    "E_15_35": (15, 35),
    "E_35_50": (35, 50),
}
ptype_labels = {"Electrons": "Electron gun", "Pions": "Pion gun"}
#  Shower profile constants (from LCShowerProfilePlugin.cc)
LONG_PROFILE_BIN_WIDTH = 0.5        # X0
LONG_PROFILE_N_BINS = 100
LONG_PROFILE_MIN_COS_ANGLE = 0.3
LONG_PROFILE_CRITICAL_ENERGY = 0.08 # GeV
LONG_PROFILE_PARAMETER_0 = 1.25
LONG_PROFILE_PARAMETER_1 = 0.5
LONG_PROFILE_MAX_DIFFERENCE = 0.1
# from steer_reco
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
    """Highest-pT charged PFO (>= 1 track, >= 1 cluster, |eta| < eta_max) in the event, or None."""
    best_pfo, best_tlv, best_pt = None, None, -1
    for pfo in pfos:
        if len(pfo.getTracks()) < 1:
            continue
        if len(pfo.getClusters()) < 1:
            continue
        tlv = getTLV(pfo)
        if abs(tlv.Eta()) > eta_max:
            continue
        pt = tlv.Perp()
        if pt > best_pt:
            best_pfo = pfo
            best_tlv = tlv
            best_pt = pt
    return best_pfo, best_tlv
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
def passFullChain(hcal_frac, rms, sstart, s_cut):
    """All three cuts, with an alternative sstart cut value."""
    if not (0 <= hcal_frac <= hcal_frac_max):
        return False
    if not (0 <= rms <= rms_max):
        return False
    if not (0 <= sstart <= s_cut):
        return False
    return True
# ---- Histograms: denominator (leading charged PFO) + cumulative stages + individual cuts ----
stage_names  = {1: "HCAL frac", 2: "+ shower RMS", 3: "+ shower start"}
single_names = {"hcal_frac": "HCAL frac only", "rms": "shower RMS only", "sstart": "shower start only"}
h_denom = ROOT.TH1F("charged_pfo_eta", ";#eta;Counts", 50, -2.5, 2.5)
h_stage = {}
for stage in stage_names:
    h_stage[stage] = ROOT.TH1F(f"stage{stage}_pfo_eta", ";#eta;Counts", 50, -2.5, 2.5)
h_single = {}
for cut in single_names:
    h_single[cut] = ROOT.TH1F(f"single_{cut}_pfo_eta", ";#eta;Counts", 50, -2.5, 2.5)
# ---- Normalized 1D distributions: variable, binning, cut value, x label, log axes ----
# Distributions are now ITERATIVE:
#   hcal_frac -> filled for all leading charged PFOs
#   rms       -> filled only for PFOs passing the HCAL fraction cut
#   sstart    -> filled only for PFOs passing the HCAL fraction + shower RMS cuts
# hcal_frac: log-spaced x bins from 1e-4 to 1, plus a dedicated "= 0" bin below 1e-4.
# Exact zeros go in the zero bin; nonzero values below 1e-4 clamp into the first log bin.
# rms and sstart: overflow clamped into the last bin so normalization counts everything.
dist_vars = {
    "hcal_frac": {"nbins": 20, "xmin": 1e-4, "xmax": 1,   "cut": hcal_frac_max, "label": "HCAL fraction",             "logy": True,  "logx": True,  "zero_bin": True,  "note": "first bin: HCAL fraction = 0", "sel": "all charged PFOs"},
    "rms":       {"nbins": 50, "xmin": 0,    "xmax": 100, "cut": rms_max,       "label": "shower RMS [mm]",           "logy": False, "logx": False, "zero_bin": False, "note": "",                             "sel": "after HCAL frac cut"},
    "sstart":    {"nbins": 40, "xmin": 0,    "xmax": 20,  "cut": sstart_max,    "label": "shower profile start [X0]", "logy": True,  "logx": False, "zero_bin": False, "note": "",                             "sel": "after HCAL frac + RMS cuts"},
}
def makeDistHist(name, v):
    if v["logx"]:
        edges = [10**(math.log10(v["xmin"]) + i*(math.log10(v["xmax"])-math.log10(v["xmin"]))/v["nbins"]) for i in range(v["nbins"]+1)]
        if v["zero_bin"]:
            edges = [v["xmin"]*0.3] + edges   # dedicated bin for exact zeros
        return ROOT.TH1F(name, "", len(edges)-1, array('d', edges))
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
    """Iterative filling: each variable only sees PFOs that survived the earlier cuts."""
    if erange is None: return
    if hcal_frac >= 0:
        if hcal_frac == 0:
            h_dist[erange][ptype]["hcal_frac"].Fill(5e-5)                          # zero bin
        else:
            h_dist[erange][ptype]["hcal_frac"].Fill(min(max(hcal_frac, 1.2e-4), 0.999))  # clamp under/overflow
    pass_hcal = (0 <= hcal_frac <= hcal_frac_max)
    if not pass_hcal: return
    if rms >= 0:
        h_dist[erange][ptype]["rms"].Fill(min(rms, 99.9))                          # clamp overflow
    pass_rms = (0 <= rms <= rms_max)
    if not pass_rms: return
    if sstart >= 0:
        h_dist[erange][ptype]["sstart"].Fill(min(sstart, 19.9))                    # clamp overflow
el_vals = []    # (hcal_frac, rms, sstart) per leading charged PFO, electron sample
pi_vals = []    # (hcal_frac, rms, sstart) per leading charged PFO, pion sample
reader = IOIMPL.LCFactory.getInstance().createLCReader()
files = sorted(samples)
print(f"found {len(files)} electron files")
n_events = 0
n_events_nopfo = 0
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
        # leading truth electron energy only sets the energy bin for the distribution plots
        truth_e = -1
        for mcp in mcps:
            if mcp.getGeneratorStatus() != 1 or abs(mcp.getPDG()) != 11: continue
            tlv = getTLV(mcp)
            if abs(tlv.Eta()) > eta_max: continue
            if tlv.E() > truth_e:
                truth_e = tlv.E()
        erange = getERangeLabel(truth_e) if truth_e >= 0 else None
        best_pfo, best_tlv = getBestPFO(pfos)
        if best_pfo is None:
            n_events_nopfo += 1
            continue
        hit_info = get_ecal_hit_info(event)
        h_denom.Fill(best_tlv.Eta())
        hcal_frac, rms, sstart = getCutVariables(best_pfo, hit_info)
        el_vals.append((hcal_frac, rms, sstart))
        fillDistHists(erange, "Electrons", hcal_frac, rms, sstart)
        for stage in stage_names:
            if passStage(stage, hcal_frac, rms, sstart):
                h_stage[stage].Fill(best_tlv.Eta())
        for cut in single_names:
            if passSingle(cut, hcal_frac, rms, sstart):
                h_single[cut].Fill(best_tlv.Eta())
    reader.close()
# ---- Pion loop: same leading charged PFO selection, for scans, fake rates, distributions ----
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
        truth_e = -1
        for mcp in mcps:
            if mcp.getGeneratorStatus() != 1 or abs(mcp.getPDG()) != 211: continue
            tlv = getTLV(mcp)
            if abs(tlv.Eta()) > eta_max: continue
            if tlv.E() > truth_e:
                truth_e = tlv.E()
        erange = getERangeLabel(truth_e) if truth_e >= 0 else None
        best_pfo, best_tlv = getBestPFO(pfos)
        if best_pfo is None:
            continue
        hit_info = get_ecal_hit_info(event)
        hcal_frac, rms, sstart = getCutVariables(best_pfo, hit_info)
        pi_vals.append((hcal_frac, rms, sstart))
        fillDistHists(erange, "Pions", hcal_frac, rms, sstart)
    reader.close()
print(f"\nProcessed {n_events} electron events, {n_pion_events} pion events")
print(f"Charged PFOs (electron sample): {len(el_vals)} | events with no charged PFO: {n_events_nopfo}")
print(f"Charged PFOs (pion sample):     {len(pi_vals)}")
var_idx = {"hcal_frac": 0, "rms": 1, "sstart": 2}
n_el_all = len(el_vals)
n_pi_all = len(pi_vals)
def subsetPassing(vals, n_prev_stage):
    """PFOs passing the first n_prev_stage cumulative cuts. n_prev_stage = 0 -> everything."""
    if n_prev_stage == 0:
        return list(vals)
    return [v for v in vals if passStage(n_prev_stage, v[0], v[1], v[2])]
# ---- Iterative threshold scans: each variable scanned only on PFOs passing the earlier cuts ----
scan_order = [
    ("hcal_frac", 0, "all charged PFOs"),
    ("rms",       1, f"charged PFOs with HCAL frac <= {hcal_frac_max}"),
    ("sstart",    2, f"charged PFOs with HCAL frac <= {hcal_frac_max} and RMS <= {rms_max}"),
]
print("\n=== Iterative Threshold Scans ===")
print("rel = fraction of the surviving sample that passes this cut")
print("abs = fraction of ALL charged PFOs that passes this cut plus the earlier ones")
for var, n_prev, desc in scan_order:
    el_sub = subsetPassing(el_vals, n_prev)
    pi_sub = subsetPassing(pi_vals, n_prev)
    print(f"\n-- {var} scanned over {desc} --")
    print(f"   electrons in sample: {len(el_sub)}/{n_el_all}   pions in sample: {len(pi_sub)}/{n_pi_all}")
    print(f"   {'cut <=':<10} | {'el eff rel':<11} | {'pi fake rel':<12} | {'el eff abs':<11} | {'pi fake abs':<12}")
    for cut in scan_values[var]:
        n_el_pass = sum(1 for v in el_sub if 0 <= v[var_idx[var]] <= cut)
        n_pi_pass = sum(1 for v in pi_sub if 0 <= v[var_idx[var]] <= cut)
        el_rel = n_el_pass/len(el_sub) if len(el_sub) > 0 else 0
        pi_rel = n_pi_pass/len(pi_sub) if len(pi_sub) > 0 else 0
        el_abs = n_el_pass/n_el_all if n_el_all > 0 else 0
        pi_abs = n_pi_pass/n_pi_all if n_pi_all > 0 else 0
        print(f"   {cut:<10} | {el_rel:<11.4f} | {pi_rel:<12.4f} | {el_abs:<11.4f} | {pi_abs:<12.4f}")
# electron efficiency and pion fake rate, denominator = all charged PFOs
print("\n--- Cumulative Cut Results (denominator = all charged PFOs) ---")
print(f"{'Stage':<18} | {'El eff':<10} | {'Pi fake rate':<12}")
for stage in stage_names:
    n_num = h_stage[stage].GetEntries()
    eff = n_num/n_el_all if n_el_all > 0 else 0
    n_pi_pass = sum(1 for v in pi_vals if passStage(stage, v[0], v[1], v[2]))
    fake = n_pi_pass/n_pi_all if n_pi_all > 0 else 0
    print(f"{stage_names[stage]:<18} | {eff:<10.4f} | {fake:<12.4f}")
#  Full chain with alternative sstart cuts
print("\n--- Full Chain vs sstart_max Variant (HCAL + RMS + sstart) ---")
print(f"{'sstart_max':<12} | {'El eff':<10} | {'Pi fake rate':<12}")
for s_cut in sstart_variants:
    n_el_pass = sum(1 for v in el_vals if passFullChain(v[0], v[1], v[2], s_cut))
    n_pi_pass = sum(1 for v in pi_vals if passFullChain(v[0], v[1], v[2], s_cut))
    eff = n_el_pass/n_el_all if n_el_all > 0 else 0
    fake = n_pi_pass/n_pi_all if n_pi_all > 0 else 0
    print(f"{s_cut:<12} | {eff:<10.4f} | {fake:<12.4f}")
print("\n--- Individual Cut Results (each cut alone, denominator = all charged PFOs) ---")
print(f"{'Cut':<18} | {'El eff':<10} | {'Pi fake rate':<12}")
for cut in single_names:
    n_num = h_single[cut].GetEntries()
    eff = n_num/n_el_all if n_el_all > 0 else 0
    n_pi_pass = sum(1 for v in pi_vals if passSingle(cut, v[0], v[1], v[2]))
    fake = n_pi_pass/n_pi_all if n_pi_all > 0 else 0
    print(f"{single_names[cut]:<18} | {eff:<10.4f} | {fake:<12.4f}")
#eff plots
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
#1D distribtion plots
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
                h.GetYaxis().SetTitle("Fraction of PFOs")
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
            leg.AddEntry(h_dist[erange][ptype][var], ptype_labels[ptype], "l")
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
        text.SetTextSize(0.03)
        text.SetTextColor(ROOT.kGray+2)
        text.DrawLatex(0.15, 0.73, v["sel"])
        if v["note"] != "":
            text.DrawLatex(0.15, 0.69, v["note"])
        text.SetTextColor(ROOT.kBlack)
        c.SaveAs(os.path.join(PLOT_DIR, f"DIST_{var}_{erange}.png"))
        c.Close()
print(f"\nPlots saved to {PLOT_DIR}")
