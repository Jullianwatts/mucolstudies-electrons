import math
import csv
import pyLCIO
from pyLCIO import EVENT, IOIMPL, UTIL

file_path = " /scratch/jwatts/mucol/v11Container/reco/pion25Fixed_reco.slcio"
#file_path = "/scratch/jwatts/mucol/v11Container/reco/1kelectron25_reco.slcio"
#file_path = "/scratch/jwatts/mucol/v11Container/reco/10kelectron0to50_reco.slcio"
#file_path = "/scratch/jwatts/mucol/v2.11/reco/electronGun_pT_0_50/electronGun_pT_0_50_reco_0.slcio"
out_csv = "shower_profile_pfos_fixed_pion.csv"
max_events = -1   # -1 = all

# ---- From LCShowerProfilePlugin.cc (your paste, confirmed defaults) ----
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
ECAL_MIP_THRESHOLD = 0.5            # ECalMipThreshold (hits below dropped by Pandora)

# ---- From MAIA_v0 ECalBarrel_o2_v01_02.xml: 50 layers, 2.20mm W + thin stuff ----
X0_PER_LAYER = 2.20/3.504 + 0.012   # ~0.64 X0 per layer

# Pandora scales the inner-layer depth by (innerPseudoLayer + 1 - pseudoLayerAtIp).
# With LCIO layer L (0-based), factor = L + INNER_LAYER_OFFSET.
# Start with 1; if validation vs dump log shows sStart off by a constant ~0.64, use 2.
INNER_LAYER_OFFSET = 1

# Pandora input collections (Sel, NOT Rec -- see ECalCaloHitCollections in steer_reco.py)
ecal_collections = ["EcalBarrelCollectionSel", "EcalEndcapCollectionSel"]

# LCElectronId cuts, for summary printout
MAX_PROFILE_START = 4.5
MAX_PROFILE_DISCREPANCY = 0.6

pdg_map = {
    11: "Electron",
    13: "Muon",
    22: "Photon",
    211: "Pion",
    2112: "Neutron"
}

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

    # Collect ECal hits passing Pandora's MIP threshold
    kept = []
    sume, sumx, sumy, sumz = 0.0, 0.0, 0.0, 0.0
    for hit in cluster.getCalorimeterHits():
        key = (hit.getCellID0(), hit.getCellID1())
        if key not in hit_info:
            continue  # not an ECal hit (HCal etc.) -- Pandora breaks at coarse granularity
        e = hit.getEnergy()
        if e * ECAL_TO_MIP < ECAL_MIP_THRESHOLD:
            continue  # below Pandora's 0.5 MIP hit threshold
        layer, is_barrel = hit_info[key]
        hp = hit.getPosition()
        kept.append((layer, is_barrel, e, (hp[0], hp[1], hp[2])))
        sume += e
        sumx += e*hp[0]
        sumy += e*hp[1]
        sumz += e*hp[2]

    if sume <= 0 or not kept:
        return None

    # Cluster direction: energy-weighted centroid (approximates Pandora fit-to-all-hits)
    norm = math.sqrt(sumx**2 + sumy**2 + sumz**2)
    if norm == 0:
        return None
    cdir = (sumx/norm, sumy/norm, sumz/norm)

    # Group by layer with angle-corrected X0 per hit
    hits_by_layer = {}
    for layer, is_barrel, e, hp in kept:
        if is_barrel:
            r = math.sqrt(hp[0]**2 + hp[1]**2)
            if r == 0:
                continue
            normal = (hp[0]/r, hp[1]/r, 0.0)  # stave normal ~ radial
        else:
            normal = (0.0, 0.0, 1.0 if hp[2] > 0 else -1.0)
        cos_angle = abs(normal[0]*cdir[0] + normal[1]*cdir[1] + normal[2]*cdir[2])
        cos_angle = max(cos_angle, LONG_PROFILE_MIN_COS_ANGLE)
        hits_by_layer.setdefault(layer, []).append((e * ECAL_TO_EM_GEV, X0_PER_LAYER / cos_angle))

    if not hits_by_layer:
        return None

    # Observed profile: walk layers, fractional spreading across bins (as in Pandora)
    profile = [0.0] * LONG_PROFILE_N_BINS
    e_ecal = 0.0
    n_rad = 0.0
    n_rad_last = 0.0
    inner_layer = min(hits_by_layer)

    for layer in range(inner_layer, max(hits_by_layer) + 1):
        if layer not in hits_by_layer:
            n_rad += n_rad_last   # empty layer: assume same material as last
            continue
        layer_hits = hits_by_layer[layer]
        e_layer = sum(h[0] for h in layer_hits)
        x0_layer = sum(h[1] for h in layer_hits) / len(layer_hits)
        e_ecal += e_layer
        n_rad_last = x0_layer
        n_rad += x0_layer

        # Account for material before start of cluster
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

    # Expected EM profile: evaluated at bin UPPER edge, as in Pandora
    a = LONG_PROFILE_PARAMETER_0 + LONG_PROFILE_PARAMETER_1 * math.log(e_ecal / LONG_PROFILE_CRITICAL_ENERGY)
    gamma_a = math.exp(math.lgamma(a))
    expected = []
    t = 0.0
    for ibin in range(LONG_PROFILE_N_BINS):
        t += LONG_PROFILE_BIN_WIDTH
        expected.append(e_ecal / 2.0 * (t/2.0)**(a - 1.0) * math.exp(-t/2.0) * LONG_PROFILE_BIN_WIDTH / gamma_a)

    # Slide and minimize, with Pandora's early exit
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

stats = {}
rows = []
event_count = 0

reader = IOIMPL.LCFactory.getInstance().createLCReader()
reader.open(file_path)

for event in reader:

    event_count += 1
    if max_events > 0 and event_count > max_events:
        break

    hit_info = get_ecal_hit_info(event)

    try:
        pfo_coll = event.getCollection("PandoraPFOs")
    except:
        continue

    for pfo in pfo_coll:

        p_type = abs(pfo.getType())
        name = pdg_map.get(p_type, f"ID({p_type})")

        for cluster in pfo.getClusters():

            result = compute_shower_profile(cluster, hit_info)
            if result is None:
                continue
            e_ecal, s_start, s_disc = result

            p = pfo.getMomentum()
            p_mag = math.sqrt(p[0]**2 + p[1]**2 + p[2]**2)
            theta = math.acos(p[2]/p_mag) if p_mag > 0 else -1

            rows.append([e_ecal, s_start, s_disc, theta, name])

            if name not in stats:
                stats[name] = [0, 0, 0]
            stats[name][0] += 1
            if s_start <= MAX_PROFILE_START:
                stats[name][1] += 1
            if s_disc <= MAX_PROFILE_DISCREPANCY:
                stats[name][2] += 1

reader.close()

with open(out_csv, "w", newline="") as fout:
    writer = csv.writer(fout)
    writer.writerow(["energy", "sStart", "sDisc", "theta", "label"])
    writer.writerows(rows)

print(f"\n({event_count} EVENTS, {len(rows)} CLUSTERS)\n")
print(f"{'Particle Type':<15} | {'Count':<8} | {'pass sStart':<12} | {'pass sDisc':<12}")

for name in sorted(stats.keys()):

    count = stats[name][0]
    f_start = stats[name][1] / count if count > 0 else 0
    f_disc = stats[name][2] / count if count > 0 else 0

    print(
        f"{name:<15} | "
        f"{count:<8} | "
        f"{f_start:<12.3f} | "
        f"{f_disc:<12.3f}"
    )

print(f"\nWrote {len(rows)} rows to {out_csv}")
