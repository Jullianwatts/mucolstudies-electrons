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
samples = glob.glob("/scratch/jwatts/mucol/v2.11/reco/electronGun_pT_0_50/electronGun_pT_0_50_reco_Final.slcio")
#samples = glob.glob("/scratch/jwatts/mucol/v11Container/reco/10kelectron0to50_reco.slcio")
#samples       = glob.glob("/scratch/jwatts/mucol/v2.11/reco/electronGun_pT_0_50/electronGun_pT_0_50_reco_0.slcio")
dr_match      = 0.2
eta_max       = 2.4
hcal_frac_max = 0.1     # raw HCAL / (raw ECAL + raw HCAL)
rms_max       = 30.0    # mm, energy-weighted transverse cluster RMS
sstart_max    = 4.5     # X0, shower profile start
sdisc_max     = 1.0     # shower profile discrepancy

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

def getCutResults(pfo, hit_info):
    """Returns dict of cut name -> pass/fail for this PFO."""
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
    cuts = {}
    cuts["hcal_frac"] = (0 <= hcal_frac <= hcal_frac_max)
    cuts["rms"]       = (0 <= rms <= rms_max)
    if result is None:
        cuts["profile"] = False
        cuts["sstart"]  = False
        cuts["sdisc"]   = False
    else:
        e_ecal, s_start, s_disc = result
        cuts["profile"] = True
        cuts["sstart"]  = (s_start <= sstart_max)
        cuts["sdisc"]   = (s_disc <= sdisc_max)
    return cuts

h_denom = ROOT.TH1F("mcp_el_eta",  ";#eta;Counts", 50, -2.5, 2.5)
h_soup  = ROOT.TH1F("soup_el_eta", ";#eta;Counts", 50, -2.5, 2.5)
h_pand  = ROOT.TH1F("pand_el_eta", ";#eta;Counts", 50, -2.5, 2.5)   # delete if soup only

cut_names = ["hcal_frac", "rms", "profile", "sstart", "sdisc"]
cutflow = {"total": 0, "no_match": 0, "pass_all": 0}
for c in cut_names:
    cutflow["fail_" + c] = 0    # failed this cut (any others may also fail)
    cutflow["only_" + c] = 0    # N-1: failed ONLY this cut

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
        hit_info = get_ecal_hit_info(event)
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
            cutflow["total"] += 1
            in_cone = [(tlv, pfo) for tlv, pfo in good_pfos if mcp_el.DeltaR(tlv) < dr_match]
            if len(in_cone) == 0:
                cutflow["no_match"] += 1
                continue
            best_tlv, best_pfo = max(in_cone, key=lambda x: x[1].getEnergy())
            cuts = getCutResults(best_pfo, hit_info)
            failed = [c for c in cut_names if not cuts[c]]
            for c in failed:
                cutflow["fail_" + c] += 1
            if len(failed) == 1:
                cutflow["only_" + failed[0]] += 1
            if len(failed) == 0:
                cutflow["pass_all"] += 1
                h_soup.Fill(mcp_el.Eta())
            if abs(best_pfo.getType()) == 11:                            # delete if soup only
                h_pand.Fill(mcp_el.Eta())                                # delete if soup only
    reader.close()
print(f"\nProcessed {n_events} events")
print("\n--- Cut Flow (truth electrons) ---")
n_tot = cutflow["total"]
print(f"{'Category':<20} | {'Count':<8} | {'Frac of total':<12}")
print(f"{'total':<20} | {n_tot:<8} | {1.0:<12.4f}")
for key in ["no_match", "pass_all"]:
    print(f"{key:<20} | {cutflow[key]:<8} | {cutflow[key]/n_tot if n_tot else 0:<12.4f}")
print("\n-- failed cut (not exclusive) --")
for c in cut_names:
    print(f"{'fail_'+c:<20} | {cutflow['fail_'+c]:<8} | {cutflow['fail_'+c]/n_tot if n_tot else 0:<12.4f}")
print("\n-- N-1: failed ONLY this cut (loosen it to gain these directly) --")
for c in cut_names:
    print(f"{'only_'+c:<20} | {cutflow['only_'+c]:<8} | {cutflow['only_'+c]/n_tot if n_tot else 0:<12.4f}")
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
