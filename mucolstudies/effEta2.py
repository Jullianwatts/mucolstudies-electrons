import os
import math
import glob
import pyLCIO
from pyLCIO import EVENT, IOIMPL
import ROOT
import plotHelper
ROOT.gROOT.SetBatch()

PLOT_DIR = "/scratch/jwatts/mucol/mucolstudies/plots2026"
os.makedirs(PLOT_DIR, exist_ok=True)

samples = glob.glob("/scratch/jwatts/mucol/v2.11/reco_v2/electrons_0_50/electron_0_50_reco.slcio")
files = {"electron_0_50": samples}

# --- Custom electron ID cuts ---
MIN_TRK_PT_GEV = 1.0   # minimum track pT to reject soft secondaries
MIN_CLU_E_GEV  = 1.0   # minimum cluster energy
MAX_WIDTH_MM   = 20.0   # transverse shower width cut; set <= 0 to disable
# ECAL-only (hcal_E == 0) always required

# --- Longitudinal profile parameters (match Pandora defaults) ---
# X0_PER_MM: MAIA ECAL has 2.2mm W absorber (X0_W = 3.504mm) in ~6mm layer pitch
# → 0.63 X0/layer ÷ 6mm/layer = 0.105 X0/mm along shower axis.
# NOTE: This is an approximation using hit positions. Pandora uses
#       getNCellRadiationLengths() from DD4hep geometry per hit.
X0_PER_MM       = 0.105   # X0 per mm depth along shower axis in MAIA ECAL barrel
PROF_BIN_WIDTH  = 0.5     # X0 per bin (Pandora default: 0.5)
PROF_N_BINS     = 100     # number of bins (covers 0–50 X0)
PROF_PARAM0     = 1.25    # Gamma function a = a0 + a1*ln(E/Ec)
PROF_PARAM1     = 0.5
PROF_CRITICAL_E = 0.08    # critical energy in silicon [GeV]
PROF_MAX_DIFF   = 0.1     # early termination threshold (Pandora default)

# Pandora cut values shown as reference lines on plots
CUT_PROFILE_START = 4.5   # m_maxProfileStart [X0]
CUT_PROFILE_DISC  = 0.6   # m_maxProfileDiscrepancy

PT_NBINS = 20
PT_MIN   = 0.0
PT_MAX   = 80.0

EFF_OBJECTS = [
    "mcp_el",
    "trk_match",
    "clu_match",
    "clu_id_match",
    "custom_el",
    "pfo_el",
]

hists = {}
for s in files:
    hists[s] = {}
    for obj in EFF_OBJECTS:
        hists[s][obj + "_eta"] = ROOT.TH1F(f"{s}_{obj}_eta", ";#eta;Counts",           50,      -2.5,   2.5)
        hists[s][obj + "_pt"]  = ROOT.TH1F(f"{s}_{obj}_pt",  ";p_{T} [GeV];Counts", PT_NBINS, PT_MIN, PT_MAX)
    hists[s]["clu_width"]  = ROOT.TH1F(f"{s}_clu_width",  ";Shower Width [mm];Clusters",      100, 0,  100)
    hists[s]["prof_start"] = ROOT.TH1F(f"{s}_prof_start", ";Profile Start [X_{0}];Clusters",  100, 0,   20)
    hists[s]["prof_disc"]  = ROOT.TH1F(f"{s}_prof_disc",  ";Profile Discrepancy;Clusters",    100, 0,    2)


def compute_width(hits):
    """Energy-weighted 2D transverse RMS (Pandora-style)."""
    sum_E = cx = cy = cz = 0.0
    for hit in hits:
        e = hit.getEnergy(); p = hit.getPosition()
        sum_E += e; cx += e*p[0]; cy += e*p[1]; cz += e*p[2]
    if sum_E <= 0: return -1.0
    cx /= sum_E; cy /= sum_E; cz /= sum_E
    mag_xy = math.sqrt(cx*cx + cy*cy)
    u_axis = ROOT.TVector3(-cy/mag_xy, cx/mag_xy, 0.0) if mag_xy > 1e-6 else ROOT.TVector3(1,0,0)
    c_unit = ROOT.TVector3(cx, cy, cz).Unit()
    v_axis = u_axis.Cross(c_unit)
    su = sv = suu = svv = 0.0
    for hit in hits:
        e = hit.getEnergy(); p = hit.getPosition()
        rel = ROOT.TVector3(p[0]-cx, p[1]-cy, p[2]-cz)
        u = rel.Dot(u_axis); v = rel.Dot(v_axis)
        su += e*u; sv += e*v; suu += e*u*u; svv += e*v*v
    ubar = su/sum_E; vbar = sv/sum_E
    return math.sqrt(max(0.0, suu/sum_E + svv/sum_E - ubar**2 - vbar**2))


def compute_longitudinal_profile(hits, cluster_energy):
    """
    Re-implementation of Pandora's longitudinal shower profile algorithm.

    Depth is approximated by projecting hit positions along the cluster direction
    (unit vector from IP to energy-weighted centroid), then converting mm → X0
    using X0_PER_MM for the MAIA ECAL.

    Returns (profileStart [X0], profileDiscrepancy).
    profileStart  = offset (in X0) at which the Gamma function best matches
                    the observed profile via a sliding L1 scan.
    profileDiscrepancy = minimum L1 difference / cluster_energy (normalized).
    """
    if not hits or cluster_energy <= 0:
        return -1.0, -1.0

    # Energy-weighted centroid
    sum_E = cx = cy = cz = 0.0
    for hit in hits:
        e = hit.getEnergy(); p = hit.getPosition()
        sum_E += e; cx += e*p[0]; cy += e*p[1]; cz += e*p[2]
    if sum_E <= 0:
        return -1.0, -1.0
    cx /= sum_E; cy /= sum_E; cz /= sum_E

    # Cluster direction: unit vector from IP toward centroid
    c_dir = ROOT.TVector3(cx, cy, cz).Unit()

    # Project each hit along cluster direction; shift so minimum depth = 0
    depths_mm = [ROOT.TVector3(h.getPosition()[0]-cx,
                               h.getPosition()[1]-cy,
                               h.getPosition()[2]-cz).Dot(c_dir)
                 for h in hits]
    min_depth = min(depths_mm)

    # Build observed profile: energy in 0.5 X0 bins
    profile = [0.0] * PROF_N_BINS
    for hit, d_mm in zip(hits, depths_mm):
        ibin = int((d_mm - min_depth) * X0_PER_MM / PROF_BIN_WIDTH)
        if 0 <= ibin < PROF_N_BINS:
            profile[ibin] += hit.getEnergy()

    # Last non-empty bin
    profile_end = next((i+1 for i in range(PROF_N_BINS-1, -1, -1) if profile[i] > 0), 0)
    if profile_end == 0:
        return -1.0, -1.0

    # Expected Gamma distribution profile: f(t) ∝ (t/2)^(a-1) * exp(-t/2)
    a = PROF_PARAM0 + PROF_PARAM1 * math.log(cluster_energy / PROF_CRITICAL_E)
    if a <= 0:
        return -1.0, -1.0
    gamma_a = math.exp(math.lgamma(a))
    expected = []
    t = 0.0
    for _ in range(PROF_N_BINS):
        t += PROF_BIN_WIDTH
        expected.append(cluster_energy/2.0 * (t/2.0)**(a-1) * math.exp(-t/2.0) * PROF_BIN_WIDTH / gamma_a)

    # Sliding offset scan: find best-fit offset minimizing L1 difference
    # Matches Pandora's CalculateLongitudinalProfile exactly.
    min_diff    = float('inf')
    best_offset = 0
    for i_off in range(profile_end):
        diff = 0.0
        for i_bin in range(profile_end):
            diff += profile[i_bin] if i_bin < i_off else abs(expected[i_bin - i_off] - profile[i_bin])
        if diff < min_diff:
            min_diff    = diff
            best_offset = i_off
        elif diff - PROF_MAX_DIFF * cluster_energy > min_diff:
            break  # Pandora early termination

    profile_start       = best_offset * PROF_BIN_WIDTH
    profile_discrepancy = min_diff / cluster_energy
    return profile_start, profile_discrepancy


def plot_with_vline(hist, save_name, xlabel, ylabel, vline_x, vline_label):
    """Plot histogram with a vertical red dashed line marking the Pandora cut."""
    can = ROOT.TCanvas("can_vl", "can_vl", 800, 600)
    can.SetLeftMargin(0.12); can.SetBottomMargin(0.12)
    can.SetTopMargin(0.12);  can.SetRightMargin(0.05)
    hist.SetLineColor(plotHelper.colors[0])
    hist.SetMarkerColor(plotHelper.colors[0])
    hist.SetMarkerStyle(20)
    hist.SetTitle("")
    hist.GetXaxis().SetTitle(xlabel); hist.GetYaxis().SetTitle(ylabel)
    hist.GetXaxis().SetTitleSize(0.045); hist.GetYaxis().SetTitleSize(0.045)
    hist.GetXaxis().SetLabelSize(0.04);  hist.GetYaxis().SetLabelSize(0.04)
    hist.SetMaximum(hist.GetMaximum() * 1.45)
    hist.Draw("E hist")
    ymax = hist.GetMaximum() * 0.88
    line = ROOT.TLine(vline_x, 0, vline_x, ymax)
    line.SetLineColor(ROOT.kRed); line.SetLineWidth(2); line.SetLineStyle(2)
    line.Draw()
    # Label in NDC
    xmin_ax = hist.GetXaxis().GetXmin(); xmax_ax = hist.GetXaxis().GetXmax()
    x_ndc = can.GetLeftMargin() + (vline_x - xmin_ax) / (xmax_ax - xmin_ax) \
            * (1 - can.GetLeftMargin() - can.GetRightMargin())
    tex = ROOT.TLatex()
    tex.SetNDC(True); tex.SetTextSize(0.033); tex.SetTextColor(ROOT.kRed)
    tex.DrawLatex(min(x_ndc + 0.01, 0.82), 0.78, vline_label)
    ROOT.gStyle.SetOptStat(0)
    can.SaveAs(save_name)
    can.Close()


reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()

prof_disc_sum   = {s: 0.0 for s in files}
prof_disc_count = {s: 0   for s in files}

for s in files:
    for f in files[s]:
        reader.open(f)
        for event in reader:

            mcps     = event.getCollection("MCParticle")
            trks     = event.getCollection("SelectedTracks")
            clusters = event.getCollection("PandoraClusters")
            pfos     = event.getCollection("PandoraPFOs")

            # Truth electrons
            mcp_electrons = []
            for mcp in mcps:
                if mcp.getGeneratorStatus() != 1 or abs(mcp.getPDG()) != 11: continue
                tlv = plotHelper.getTLV(mcp)
                if abs(tlv.Eta()) > 2.4: continue
                hists[s]["mcp_el_eta"].Fill(tlv.Eta())
                hists[s]["mcp_el_pt"].Fill(tlv.Perp())
                mcp_electrons.append(tlv)

            trk_tlvs      = [plotHelper.getTrackTLV(t) for t in trks]
            cluster_pairs = [(c, plotHelper.getClusterTLV(c)) for c in clusters]
            pfo_el_tlvs   = [plotHelper.getTLV(p) for p in pfos if abs(p.getType()) == 11]

            for mcp_el in mcp_electrons:
                m_eta = mcp_el.Eta(); m_pt = mcp_el.Perp()

                def fill(key):
                    hists[s][key+"_eta"].Fill(m_eta)
                    hists[s][key+"_pt"].Fill(m_pt)

                # Track
                trks_in_cone = [t for t in trk_tlvs
                                if mcp_el.DeltaR(t) < 0.2 and t.Perp() > MIN_TRK_PT_GEV]
                has_track = len(trks_in_cone) > 0

                # Raw cluster match: highest E in cone
                pairs_in_cone = [(c_obj, c_tlv) for c_obj, c_tlv in cluster_pairs
                                 if mcp_el.DeltaR(c_tlv) < 0.2]
                best_pair   = max(pairs_in_cone, key=lambda x: x[0].getEnergy(), default=None)
                has_cluster = best_pair is not None

                # Longitudinal profile on best matched cluster (no ID cuts applied)
                if has_cluster:
                    c_obj, c_tlv = best_pair
                    hits     = list(c_obj.getCalorimeterHits())
                    p_start, p_disc = compute_longitudinal_profile(hits, c_obj.getEnergy())
                    if p_start >= 0:
                        hists[s]["prof_start"].Fill(p_start)
                    if p_disc >= 0:
                        hists[s]["prof_disc"].Fill(p_disc)
                        prof_disc_sum[s]   += p_disc
                        prof_disc_count[s] += 1

                # Custom cluster ID: filter first, then pick highest E
                custom_candidates = []
                for c_obj, c_tlv in pairs_in_cone:
                    if c_obj.getEnergy() < MIN_CLU_E_GEV: continue
                    sub_E  = c_obj.getSubdetectorEnergies()
                    hcal_E = sub_E[1] if len(sub_E) > 1 else 0.0
                    if hcal_E != 0.0: continue
                    hits  = list(c_obj.getCalorimeterHits())
                    width = compute_width(hits)
                    hists[s]["clu_width"].Fill(width)
                    if MAX_WIDTH_MM > 0 and (width < 0 or width >= MAX_WIDTH_MM): continue
                    custom_candidates.append((c_obj, c_tlv))
                has_custom_cluster = len(custom_candidates) > 0

                # PFO type=11 (comparison only)
                has_pfo_el = any(mcp_el.DeltaR(p) < 0.2 for p in pfo_el_tlvs)

                if has_track:                        fill("trk_match")
                if has_cluster:                      fill("clu_match")
                if has_custom_cluster:               fill("clu_id_match")
                if has_track and has_custom_cluster: fill("custom_el")
                if has_pfo_el:                       fill("pfo_el")

        reader.close()

# --- Print summary ---
print("\n--- Efficiency Results ---")
for s in hists:
    display_names = {
        "trk_match":    "Track in cone",
        "clu_match":    "Cluster in cone (no cuts)",
        "clu_id_match": "Cluster passing custom ID cuts",
        "custom_el":    "Track + Cluster (custom ID)",
        "pfo_el":       "PandoraPFOs type=11 [Pandora]",
    }
    denom_eta = hists[s]["mcp_el_eta"]
    denom_pt  = hists[s]["mcp_el_pt"]
    n_den     = denom_eta.GetEntries()
    print(f"\n[{s}]  {int(n_den)} truth electrons")
    print(f"  Cuts: trk pT > {MIN_TRK_PT_GEV} GeV | clu E > {MIN_CLU_E_GEV} GeV | width < {MAX_WIDTH_MM} mm | ECAL-only")

    # Profile diagnostics
    if prof_disc_count[s] > 0:
        avg_disc = prof_disc_sum[s] / prof_disc_count[s]
        h_disc   = hists[s]["prof_disc"]
        h_start  = hists[s]["prof_start"]
        frac_disc_above  = h_disc.Integral(h_disc.FindBin(CUT_PROFILE_DISC), h_disc.GetNbinsX()) \
                           / max(prof_disc_count[s], 1)
        frac_start_above = h_start.Integral(h_start.FindBin(CUT_PROFILE_START), h_start.GetNbinsX()) \
                           / max(h_start.GetEntries(), 1)
        print(f"\n  --- Longitudinal Profile Diagnostics (approx, X0_PER_MM={X0_PER_MM}) ---")
        print(f"  Average profile discrepancy : {avg_disc:.4f}  (Pandora cut: <= {CUT_PROFILE_DISC})")
        print(f"  Fraction failing discrepancy cut : {frac_disc_above:.4f}")
        print(f"  Fraction failing profileStart cut (> {CUT_PROFILE_START} X0): {frac_start_above:.4f}")

    eff_map_eta = {}; eff_map_pt = {}
    for obj, label in display_names.items():
        num_eta = hists[s][obj+"_eta"]
        n_num   = num_eta.GetEntries()
        eff     = n_num / n_den if n_den > 0 else 0
        print(f"  {label:<40} | Eff: {eff:.4f} ({int(n_num)}/{int(n_den)})")
        teff_eta = ROOT.TEfficiency(num_eta, denom_eta)
        eff_map_eta[label] = teff_eta.CreateGraph()
        teff_pt  = ROOT.TEfficiency(hists[s][obj+"_pt"], denom_pt)
        eff_map_pt[label]  = teff_pt.CreateGraph()

    plotHelper.plotEfficiencies(eff_map_eta,
        os.path.join(PLOT_DIR, "EFFICIENCY_CUSTOM_ETA.png"),
        xlabel="#eta", ylabel="Efficiency")
    plotHelper.plotEfficiencies(eff_map_pt,
        os.path.join(PLOT_DIR, "EFFICIENCY_CUSTOM_PT.png"),
        xlabel="p_{T} [GeV]", ylabel="Efficiency")
    plotHelper.plotHistograms({s: hists[s]["clu_width"]},
        os.path.join(PLOT_DIR, "SHOWER_WIDTH.png"),
        xlabel="Shower Width [mm]", ylabel="Clusters / bin")
    plot_with_vline(hists[s]["prof_start"],
        os.path.join(PLOT_DIR, "PROFILE_START.png"),
        "Profile Start [X_{0}]", "Clusters / bin",
        CUT_PROFILE_START, f"Pandora cut: {CUT_PROFILE_START} X_{{0}}")
    plot_with_vline(hists[s]["prof_disc"],
        os.path.join(PLOT_DIR, "PROFILE_DISCREPANCY.png"),
        "Profile Discrepancy", "Clusters / bin",
        CUT_PROFILE_DISC, f"Pandora cut: {CUT_PROFILE_DISC}")

print(f"\nPlots saved to {PLOT_DIR}")
