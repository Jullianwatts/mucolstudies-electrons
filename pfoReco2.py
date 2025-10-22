import os
import math
import ROOT
import glob
import pyLCIO
from pyLCIO import IOIMPL

exec(open("./plotHelper.py").read())
ROOT.gROOT.SetBatch()

# Make sure plots folder exists
os.makedirs("plots", exist_ok=True)

# Solenoid field (tune this with your printed combined mean)
B_FIELD = 2.65  # Tesla

max_events = -1

samples = glob.glob("/data/fmeloni/DataMuC_MAIA_v0/v5/reco/electronGun*")

files = {}
slices = ["0_50", "50_250", "250_1000", "1000_5000"]
for s in slices:
    files[f"electronGun_pT_{s}"] = []

for s in samples:
    sname = s.split("/")[-1]
    if sname not in files:
        continue
    files[sname] = glob.glob(f"{s}/*.slcio")

match_counts = {s: 0 for s in files}

# --- matching & E/p helpers ---

def match_track_to_cluster(track_tlv, cluster_tlv, dR_cut=0.2):
    """Return True if track and cluster are matched within ΔR."""
    if track_tlv.E() <= 0 or cluster_tlv.E() <= 0:
        return False
    return track_tlv.DeltaR(cluster_tlv) < dR_cut

def calculate_e_over_p_metrics(cluster_energy, track_momentum):
    """Return (E/p, (E-p)/p, |E-p|/p). For invalid p, return sentinels."""
    if track_momentum <= 0:
        return -1.0, -999.0, -1.0
    e_over_p = cluster_energy / track_momentum
    e_minus_p_over_p = (cluster_energy - track_momentum) / track_momentum
    abs_e_minus_p_over_p = abs(e_minus_p_over_p)
    return e_over_p, e_minus_p_over_p, abs_e_minus_p_over_p

def makeTrackTLV_from_params(trk, b_field, m=0.000511):
    """
    Build a TLorentzVector for the track using LCIO parameters.
    pT = 0.0003 * B / |omega|  (GeV),  pz = pT * tanLambda
    """
    omega = trk.getOmega()
    if abs(omega) < 1e-12:
        tlv = ROOT.TLorentzVector(); tlv.SetPxPyPzE(0, 0, 0, m)
        return tlv
    phi  = trk.getPhi()        # LCIO phi0
    tanL = trk.getTanLambda()  # LCIO tan(lambda)

    pT = 0.0003 * b_field / abs(omega)  # GeV
    pz = pT * tanL
    p  = math.hypot(pT, pz)
    px = pT * math.cos(phi)
    py = pT * math.sin(phi)
    E  = math.sqrt(p * p + m * m)

    tlv = ROOT.TLorentzVector()
    tlv.SetPxPyPzE(px, py, pz, E)
    return tlv

def find_closest_track(cluster_tlv, track_list, dR_cut=0.2):
    if cluster_tlv.E() <= 0:
        return None
    closest_track = None
    min_dR = float('inf')
    for track_tlv in track_list:
        if track_tlv.E() <= 0:
            continue
        dR = cluster_tlv.DeltaR(track_tlv)
        if dR < dR_cut and dR < min_dR:
            min_dR = dR
            closest_track = track_tlv
    return closest_track

# --- hist containers ---
hists = {}
hists2d = {}

for s in files:
    hists[s] = {}
    hists2d[s] = {}
    for obj in ["pfo","pfo_e","pfo_e_match", "mcp", "mcp_e", "trk", "trk_e", "trk_e_match", "clusters", "clusters_e_match"]:
        for vtype in ["obj", "evt"]:
            for var in variables[vtype]:
                hists[s][obj+"_"+var] = ROOT.TH1F(
                    s+"_"+obj+"_"+var, s,
                    variables[vtype][var]["nbins"],
                    variables[vtype][var]["xmin"],
                    variables[vtype][var]["xmax"]
                )

    # E/p ratio hists (existing style)
    hists[s]["pfo_e_match_obj_ep_ratio"]      = ROOT.TH1F(s+"_pfo_e_match_obj_ep_ratio", s, 100, 0, 3)
    hists[s]["clusters_e_match_obj_ep_ratio"] = ROOT.TH1F(s+"_clusters_e_match_obj_ep_ratio", s, 100, 0, 3)

    # New: explicit E/p analysis hists (like your older block)
    hists[s]["electron_E_over_p"]             = ROOT.TH1F(f"{s}_electron_E_over_p",             f"E/p Distribution {s};E/p;Entries", 100, 0, 3)
    hists[s]["electron_E_minus_p_over_p"]     = ROOT.TH1F(f"{s}_electron_E_minus_p_over_p",     f"(E-p)/p Distribution {s};(E-p)/p;Entries", 100, -1, 1)
    hists[s]["electron_abs_E_minus_p_over_p"] = ROOT.TH1F(f"{s}_electron_abs_E_minus_p_over_p", f"|E-p|/p Distribution {s};|E-p|/p;Entries", 100, 0, 2)

reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()
reader.setReadCollectionNames(["MCParticle","SiTracks", "AllTracks", "PandoraPFOs", "SeedTracks","PandoraClusters"])

for s in files:
    print(f"Working on sample {s}")
    i = 0

    for f in files[s]:
        if max_events > 0 and i >= max_events:
            break

        reader.open(f)

        for event in reader:
            if max_events > 0 and i >= max_events:
                break
            if i % 100 == 0:
                print(f"\tProcessing event: {i}")

            try:
                mcps = event.getCollection("MCParticle")
            except:
                mcps = []
                print("No MCP")

            try:
                clusters = event.getCollection("PandoraClusters")
                tracks   = event.getCollection("SiTracks")
            except:
                i += 1
                continue

            # Build TLVs
            cluster_tlvs = []
            track_tlvs   = []

            for cluster in clusters:
                cluster_tlv = getClusterTLV(cluster)
                if cluster_tlv.E() < 10:
                    continue
                cluster_tlvs.append(cluster_tlv)

            for track in tracks:
                # Track TLV from LCIO params (omega, phi, tanLambda) with chosen B_FIELD
                track_tlv = makeTrackTLV_from_params(track, b_field=B_FIELD, m=0.000511)
                if track_tlv.E() < 10:
                    continue
                track_tlvs.append(track_tlv)

            # Debug (first 100 events)
            if i < 100:
                print(f"\n--- {s} Event {i} ---")
                print(f"Clusters ({len(cluster_tlvs)} total):")
                for j, cluster_tlv in enumerate(cluster_tlvs):
                    print(f"  Cluster {j}: E = {cluster_tlv.E():.2f} GeV")
                print(f"Tracks ({len(track_tlvs)} total):")
                for j, track_tlv in enumerate(track_tlvs):
                    print(f"  Track {j}: p = {track_tlv.P():.2f} GeV")

            # E/p filling: match cluster ↔ closest track and fill both styles of E/p hists
            for cluster_tlv in cluster_tlvs:
                closest_track = find_closest_track(cluster_tlv, track_tlvs, dR_cut=0.2)
                if closest_track is None:
                    continue
                if not match_track_to_cluster(closest_track, cluster_tlv, dR_cut=0.2):
                    continue

                match_counts[s] += 1

                E = cluster_tlv.E()
                p = closest_track.P()
                if p <= 0:
                    continue

                # explicit metrics
                ep_ratio, e_minus_p_over_p, abs_e_minus_p_over_p = calculate_e_over_p_metrics(E, p)

                # fill both old + new hists
                fillObjHists(hists[s], "pfo_e_match",      closest_track)
                fillObjHists(hists[s], "clusters_e_match", cluster_tlv)

                hists[s]["pfo_e_match_obj_ep_ratio"].Fill(ep_ratio)
                hists[s]["clusters_e_match_obj_ep_ratio"].Fill(ep_ratio)

                hists[s]["electron_E_over_p"].Fill(ep_ratio)
                hists[s]["electron_E_minus_p_over_p"].Fill(e_minus_p_over_p)
                hists[s]["electron_abs_E_minus_p_over_p"].Fill(abs_e_minus_p_over_p)

                if i < 10:
                    print(f"{s}: E={E:.1f} GeV, p={p:.1f} GeV, E/p={ep_ratio:.3f}")

            i += 1

        reader.close()

print(f"Total matches: {sum(match_counts.values())}")

# Labels keyed by electron samples
label_map = {
    "electronGun_pT_0_50": "0-50 GeV",
    "electronGun_pT_50_250": "50-250 GeV",
    "electronGun_pT_250_1000": "250-1000 GeV",
    "electronGun_pT_1000_5000": "1000-5000 GeV"
}

# --- Plot E/p ratios comparison (existing style) ---
ep_hists = {}
for sample_name in files:
    if "pfo_e_match_obj_ep_ratio" in hists[sample_name] and hists[sample_name]["pfo_e_match_obj_ep_ratio"].GetEntries() > 0:
        ep_hists[sample_name] = hists[sample_name]["pfo_e_match_obj_ep_ratio"]

if ep_hists:
    plotHistograms(
        ep_hists,
        "plots/electron_ep_ratio_comparison.png",
        xlabel="E/p", ylabel="Entries",
        atltext=["Muon Collider", "Simulation, no BIB", "|#eta| < 2.4", "MAIA Detector Concept"]
    )

# --- Plot variables for matched electrons ---
for var in ["pt", "eta", "phi", "E"]:
    var_hists = {}
    for sample_name in files:
        hist_name = f"pfo_e_match_obj_{var}"
        if hist_name in hists[sample_name] and hists[sample_name][hist_name].GetEntries() > 0:
            var_hists[sample_name] = hists[sample_name][hist_name]

    if var_hists:
        plotHistograms(
            var_hists,
            f"plots/pfo_e_match_{var}.png",
            xlabel=variables["obj"][var]["label"], ylabel="Entries",
            atltext=["Muon Collider", "Simulation, no BIB", "|#eta| < 2.4", "MAIA Detector Concept"]
        )

# --- Plot cluster E/p ratios (existing style) ---
cluster_ep_hists = {}
for sample_name in files:
    hist_name = "clusters_e_match_obj_ep_ratio"
    if hist_name in hists[sample_name] and hists[sample_name][hist_name].GetEntries() > 0:
        cluster_ep_hists[sample_name] = hists[sample_name][hist_name]

if cluster_ep_hists:
    plotHistograms(
        cluster_ep_hists,
        "plots/electron_cluster_ep_ratio_comparison.png",
        xlabel="E/p", ylabel="Entries",
        atltext=["Muon Collider", "Simulation, no BIB", "|#eta| < 2.4", "MAIA Detector Concept"]
    )

# --- Optional: save the new E/p histograms too ---
new_ep_hists = {}
new_epres_hists = {}
new_epabs_hists = {}
for sample_name in files:
    if hists[sample_name]["electron_E_over_p"].GetEntries() > 0:
        new_ep_hists[sample_name] = hists[sample_name]["electron_E_over_p"]
    if hists[sample_name]["electron_E_minus_p_over_p"].GetEntries() > 0:
        new_epres_hists[sample_name] = hists[sample_name]["electron_E_minus_p_over_p"]
    if hists[sample_name]["electron_abs_E_minus_p_over_p"].GetEntries() > 0:
        new_epabs_hists[sample_name] = hists[sample_name]["electron_abs_E_minus_p_over_p"]

if new_ep_hists:
    plotHistograms(
        new_ep_hists,
        "plots/electron_E_over_p.png",
        xlabel="E/p", ylabel="Entries",
        atltext=["Muon Collider", "Simulation, no BIB", "|#eta| < 2.4", "MAIA Detector Concept"]
    )
if new_epres_hists:
    plotHistograms(
        new_epres_hists,
        "plots/electron_E_minus_p_over_p.png",
        xlabel="(E-p)/p", ylabel="Entries",
        atltext=["Muon Collider", "Simulation, no BIB", "|#eta| < 2.4", "MAIA Detector Concept"]
    )
if new_epabs_hists:
    plotHistograms(
        new_epabs_hists,
        "plots/electron_abs_E_minus_p_over_p.png",
        xlabel="|E-p|/p", ylabel="Entries",
        atltext=["Muon Collider", "Simulation, no BIB", "|#eta| < 2.4", "MAIA Detector Concept"]
    )

# --- Print stats (existing) ---
for sample_name in files:
    if sample_name not in hists or match_counts[sample_name] == 0:
        continue
    print(f"\n{sample_name}:")
    print(f"  Total matches: {match_counts[sample_name]}")

    ep_hist = hists[sample_name]["pfo_e_match_obj_ep_ratio"]
    if ep_hist.GetEntries() > 0:
        total_entries = ep_hist.GetEntries()
        near_one = 0
        for i in range(1, ep_hist.GetNbinsX() + 1):
            bin_center = ep_hist.GetBinCenter(i)
            if 0.8 <= bin_center <= 1.2:
                near_one += ep_hist.GetBinContent(i)
        percent_near_one = (near_one / total_entries) * 100 if total_entries > 0 else 0
        print(f"  E/p entries with 0.8 < E/p < 1.2: {percent_near_one:.1f}%")

# --- Save all 2D histograms if present ---
for s in hists2d:
    for h in hists2d[s]:
        c = ROOT.TCanvas("c", "c", 800, 600)
        hists2d[s][h].Draw("COLZ")
        c.SaveAs(f"plots/{hists2d[s][h].GetName()}.png")
        c.Close()

# --- Combined E/p suggestion (weighted) to refine B_FIELD ---
target = 1.5
total_entries = 0
weighted_mean_sum = 0.0
for sample_name in files:
    h = hists.get(sample_name, {}).get("pfo_e_match_obj_ep_ratio")
    if not h:
        continue
    n = int(h.GetEntries())
    if n == 0:
        continue
    weighted_mean_sum += h.GetMean() * n
    total_entries += n

if total_entries > 0:
    combined_mean = weighted_mean_sum / total_entries
    suggested_B = B_FIELD * (combined_mean / target)
    print(f"[ALL] combined entries={total_entries}  mean(E/p)={combined_mean:.3f}")
    print(f"[ALL] Suggested single B_FIELD ≈ {suggested_B:.2f} T  (current {B_FIELD:.2f} T, target {target})")
else:
    print("[ALL] No E/p entries to compute combined mean. Increase max_events or loosen matching.")

