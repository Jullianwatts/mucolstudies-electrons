import math
import pyLCIO
from pyLCIO import EVENT, IOIMPL, UTIL
import ROOT
exec(open("./plotHelper.py").read())
ROOT.gROOT.SetBatch(True)
#change script
file_path = "/scratch/jwatts/mucol/v11Container/reco/electron25Fixed_reco.slcio"

X0_PER_LAYER  = 2.2 / 3.504

CUT_MAX_ENERGY        = 5.0
CUT_MAX_INNER_LAYER   = 4
CUT_MAX_PROFILE_START = 4.5
CUT_MAX_EOP_RESIDUAL  = 0.2

reader = IOIMPL.LCFactory.getInstance().createLCReader()
reader.open(file_path)

encoding = None
for event in reader:
    for coll_name in ["EcalBarrelCollectionRec", "ECALBarrel", "ECalBarrel",
                      "EcalBarrelCollection", "HCALBarrel", "HCalBarrelCollection"]:
        try:
            coll     = event.getCollection(coll_name)
            encoding = coll.getParameters().getStringVal("CellIDEncoding")
            print(f"Got CellIDEncoding from '{coll_name}': {encoding}")
            break
        except:
            continue
    if encoding: break
reader.close()

if encoding is None:
    raise RuntimeError("Could not find CellIDEncoding.")

decoder = UTIL.BitField64(encoding)

def get_layer(hit):
    decoder.setValue(hit.getCellID0())
    return int(decoder["layer"])

def get_system(hit):
    decoder.setValue(hit.getCellID0())
    return int(decoder["system"])

# scan all system IDs across 1000 events so we can set ECAL_SYSTEM_IDS correctly
print("Scanning system IDs across 1000 events...")
all_system_ids = set()
scan_reader = IOIMPL.LCFactory.getInstance().createLCReader()
scan_reader.open(file_path)
for i, event in enumerate(scan_reader):
    if i >= 1000: break
    try:
        pfos = event.getCollection("PandoraPFOs")
    except:
        continue
    for pfo in pfos:
        for cluster in pfo.getClusters():
            for hit in cluster.getCalorimeterHits():
                all_system_ids.add(get_system(hit))
scan_reader.close()
print(f"All system IDs in PFO cluster hits: {sorted(all_system_ids)}")
print(f"-> update ECAL_SYSTEM_IDS below if needed\n")

# ---- update this once you see the system IDs above ----
ECAL_SYSTEM_IDS = {29}
# -------------------------------------------------------

reader = IOIMPL.LCFactory.getInstance().createLCReader()
reader.open(file_path)

event_count = 0
n_matched   = 0
records     = []

for event in reader:
    if event_count >= 1000: break
    event_count += 1

    try:
        pfos = event.getCollection("PandoraPFOs")
        mcps = event.getCollection("MCParticle")
    except:
        continue

    truth_e = [getTLV(m) for m in mcps if abs(m.getPDG()) == 11 and m.getGeneratorStatus() == 1]
    if not truth_e: continue

    for pfo in pfos:

        tracks   = pfo.getTracks()
        clusters = pfo.getClusters()
        if len(tracks) == 0 or len(clusters) == 0: continue

        cluster = clusters[0]
        track   = tracks[0]
        if track.getOmega() == 0: continue

        trk_tlv = getTrackTLV(track)
        p = trk_tlv.P()
        if p < 0.5: continue

        pfo_p   = pfo.getMomentum()
        pfo_tlv = ROOT.TLorentzVector(pfo_p[0], pfo_p[1], pfo_p[2], pfo.getEnergy())
        if not any(pfo_tlv.DeltaR(e) < 0.1 for e in truth_e): continue

        n_matched += 1

        em_energy = cluster.getEnergy()
        ep        = em_energy / p

        sub_e       = cluster.getSubdetectorEnergies()
        ecal_energy = sub_e[0] if len(sub_e) > 0 else 0.0
        hcal_energy = sub_e[1] if len(sub_e) > 1 else 0.0
        total_sub_e = ecal_energy + hcal_energy
        ecal_frac   = ecal_energy / total_sub_e if total_sub_e > 0 else 0.0

        layer_energy = {}
        for hit in cluster.getCalorimeterHits():
            if get_system(hit) in ECAL_SYSTEM_IDS:
                layer = get_layer(hit)
                layer_energy[layer] = layer_energy.get(layer, 0.0) + hit.getEnergy()

        if layer_energy:
            min_ecal_layer     = min(layer_energy.keys())
            inner_pseudo_layer = min_ecal_layer + 1
            profile_start_est  = min_ecal_layer * X0_PER_LAYER
        else:
            inner_pseudo_layer = 9999
            profile_start_est  = 9999.0

        records.append({
            "em_energy":          em_energy,
            "ep":                 ep,
            "abs_ep_minus_1":     abs(ep - 1),
            "inner_pseudo_layer": inner_pseudo_layer,
            "profile_start_est":  profile_start_est,
            "ecal_energy":        ecal_energy,
            "hcal_energy":        hcal_energy,
            "ecal_frac":          ecal_frac,
            "pfo_type":           abs(pfo.getType()),
        })

reader.close()

def stats(vals):
    if not vals: return "  no data"
    mean = sum(vals) / len(vals)
    return f"  mean={mean:8.3f}  min={min(vals):8.3f}  max={max(vals):8.3f}"

def extract(recs, key):
    return [r[key] for r in recs]

def n_fail(vals, cut):
    return sum(1 for v in vals if v > cut)

def n_fail_abs(vals, cut):
    return sum(1 for v in vals if abs(v - 1) > cut)

called_e = [r for r in records if r["pfo_type"] == 11]
missed_e = [r for r in records if r["pfo_type"] != 11]

print(f"\n{event_count} events | {n_matched} truth-matched electron PFOs")
print(f"  Pandora called electron (type=11): {len(called_e)} / {n_matched}")
print(f"  Pandora did NOT call electron:     {len(missed_e)} / {n_matched}")

for label, recs in [("CALLED ELECTRON (type=11)", called_e),
                    ("NOT CALLED ELECTRON",        missed_e)]:

    n = len(recs)
    if n == 0:
        print(f"\n--- {label} : none ---")
        continue

    em  = extract(recs, "em_energy")
    ep  = extract(recs, "ep")
    aep = extract(recs, "abs_ep_minus_1")
    il  = extract(recs, "inner_pseudo_layer")
    ps  = extract(recs, "profile_start_est")
    ef  = extract(recs, "ecal_frac")
    ee  = extract(recs, "ecal_energy")
    he  = extract(recs, "hcal_energy")

    print(f"\n{'='*70}")
    print(f"  {label}  ({n} PFOs)")
    print(f"{'='*70}")
    print(f"  {'Variable':<44} Stats")
    print(f"  {'-'*67}")

    print(f"  {'em_energy [GeV]  (cut: >5 if !IsEmShower)':<44}{stats(em)}")
    print(f"      -> {n_fail(em, CUT_MAX_ENERGY)}/{n} exceed 5 GeV")

    print(f"  {'inner_pseudo_layer  (cut: >4 if !IsEmShower)':<44}{stats(il)}")
    print(f"      pseudo-layer 4 = ECal layer 3 = {3*X0_PER_LAYER:.2f} X0 from ECal face")
    print(f"      -> {n_fail(il, CUT_MAX_INNER_LAYER)}/{n} exceed pseudo-layer 4")
    print(f"      NOTE: 9999 means no ECal hits found with ECAL_SYSTEM_IDS={ECAL_SYSTEM_IDS}")

    print(f"  {'profile_start_est [X0]  (cut: >4.5)':<44}{stats(ps)}")
    print(f"      = min_ecal_layer * 0.628  (rough proxy; Pandora uses template fit)")
    print(f"      -> {n_fail(ps, CUT_MAX_PROFILE_START)}/{n} estimated exceed 4.5 X0")

    print(f"  {'E/p':<44}{stats(ep)}")
    print(f"  {'|E/p - 1|  (cut: >0.2 kills)':<44}{stats(aep)}")
    print(f"      -> {n_fail_abs(ep, CUT_MAX_EOP_RESIDUAL)}/{n} fail E/p cut")

    print(f"  {'ECal energy [GeV]  (subdetectorEnergies[0])':<44}{stats(ee)}")
    print(f"  {'HCal energy [GeV]  (subdetectorEnergies[1])':<44}{stats(he)}")
    print(f"  {'ECal energy fraction':<44}{stats(ef)}")
