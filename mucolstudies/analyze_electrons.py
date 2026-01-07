import os, math, glob
import ROOT, pyLCIO
from pyLCIO import IOIMPL
from array import array

ROOT.gROOT.SetBatch(True)
exec(open("plotHelper.py").read())

# ------------------
# Config
# ------------------
B_FIELD    = 5.0
MAX_EVENTS = -1
MIN_E_GEV  = 10.0
DR_CUT     = 0.2

files = {
    "v5_electronGun": ["/data/fmeloni/DataMuC_MAIA_v0/v5/reco/electronGun_pT_0_50/*.slcio"
        ]
}


OUTDIR = "mucolstudies/recoPlots/v5"
os.makedirs(OUTDIR, exist_ok=True)

# ------------------
# Utilities
# ------------------
def expand_file_list(patterns):
    out = []
    for p in patterns:
        out += sorted(glob.glob(p))
    return out

def dphi(phi1, phi2):
    d = phi1 - phi2
    while d >  math.pi: d -= 2*math.pi
    while d < -math.pi: d += 2*math.pi
    return d

def deltaR(t1, t2):
    return math.sqrt((t1.Eta()-t2.Eta())**2 + dphi(t1.Phi(), t2.Phi())**2)

def makeTrackTLV(trk, m=0.000511):
    if trk.getOmega() == 0:
        return None
    pt = 3e-4 * abs(B_FIELD / trk.getOmega())
    px = pt * math.cos(trk.getPhi())
    py = pt * math.sin(trk.getPhi())
    pz = pt * trk.getTanLambda()
    p  = math.sqrt(px*px + py*py + pz*pz)
    E  = math.sqrt(p*p + m*m)
    tlv = ROOT.TLorentzVector()
    tlv.SetPxPyPzE(px, py, pz, E)
    return tlv

def find_closest(tlv, cands):
    best, best_dr = None, 1e9
    for c in cands:
        dr = deltaR(tlv, c)
        if dr < best_dr:
            best, best_dr = c, dr
    return best if best_dr < DR_CUT else None

def makeMCTLV(mcp, m=0.000511):
    p3 = mcp.getMomentum()
    px = float(p3[0])
    py = float(p3[1])
    pz = float(p3[2])
    p  = math.sqrt(px*px + py*py + pz*pz)
    E  = math.sqrt(p*p + m*m)
    tlv = ROOT.TLorentzVector()
    tlv.SetPxPyPzE(px, py, pz, E)
    return tlv

def pick_truth_electron(mcps):
    best, best_pt = None, -1
    for mcp in mcps:
        if abs(int(mcp.getPDG())) != 11:
            continue
        tlv = makeMCTLV(mcp)
        if tlv.Pt() > best_pt:
            best, best_pt = tlv, tlv.Pt()
    return best

# ------------------
# Histograms
# ------------------
hists = {}

# pT bins capped at 50 GeV
pt_bins = array("d", [0, 10, 20, 30, 40, 50])
nb_pt = len(pt_bins) - 1

# eta bins
eta_bins = array("d", [-2.4, -1.6, -0.8, 0.0, 0.8, 1.6, 2.4])
nb_eta = len(eta_bins) - 1

for s in files:
    hists[s] = {}

    # Resolution / validation
    hists[s]["E_over_p"] = ROOT.TH1F(f"{s}_E_over_p", "", 80, 0, 2.0)
    hists[s]["pt_reco_over_true"] = ROOT.TH1F(f"{s}_pt_reco_over_true", "", 120, 0, 2.0)
    hists[s]["pt_reco_minus_true_over_true"] = ROOT.TH1F(
        f"{s}_pt_reco_minus_true_over_true", "", 160, -1, 1
    )

    # Efficiency bookkeeping (ETA)
    hists[s]["truth_eta"]    = ROOT.TH1F(f"{s}_truth_eta", "", nb_eta, eta_bins)
    hists[s]["track_eta"]    = ROOT.TH1F(f"{s}_track_eta", "", nb_eta, eta_bins)
    hists[s]["cluster_eta"]  = ROOT.TH1F(f"{s}_cluster_eta", "", nb_eta, eta_bins)
    hists[s]["electron_eta"] = ROOT.TH1F(f"{s}_electron_eta","", nb_eta, eta_bins)

# ------------------
# Event loop
# ------------------
reader = IOIMPL.LCFactory.getInstance().createLCReader()
reader.setReadCollectionNames(["MCParticle", "SiTracks", "AllTracks", "PandoraClusters"])

for s in files:
    for fpath in expand_file_list(files[s]):
        reader.open(fpath)

        for i, event in enumerate(reader):
            if MAX_EVENTS > 0 and i >= MAX_EVENTS:
                break

            try:
                clusters = event.getCollection("PandoraClusters")
                tracks   = event.getCollection("SiTracks")
                mcps     = event.getCollection("MCParticle")
            except:
                continue

            cluster_tlvs = []
            for c in clusters:
                ct = getClusterTLV(c)
                if ct.E() > MIN_E_GEV:
                    cluster_tlvs.append(ct)

            track_tlvs = []
            for t in tracks:
                tt = makeTrackTLV(t)
                if tt and tt.E() > MIN_E_GEV:
                    track_tlvs.append(tt)

            truth = pick_truth_electron(mcps)
            if not truth:
                continue

            eta_true = truth.Eta()
            hists[s]["truth_eta"].Fill(eta_true)

            reco_trk = find_closest(truth, track_tlvs)
            reco_cl  = find_closest(truth, cluster_tlvs)

            if reco_trk:
                hists[s]["track_eta"].Fill(eta_true)
            if reco_cl:
                hists[s]["cluster_eta"].Fill(eta_true)
            if reco_trk and reco_cl:
                hists[s]["electron_eta"].Fill(eta_true)
                hists[s]["E_over_p"].Fill(reco_cl.E() / reco_trk.P())

            if reco_trk and truth.Pt() > 0:
                hists[s]["pt_reco_over_true"].Fill(reco_trk.Pt() / truth.Pt())
                hists[s]["pt_reco_minus_true_over_true"].Fill(
                    (reco_trk.Pt() - truth.Pt()) / truth.Pt()
                )

        reader.close()

# ------------------
# Build efficiencies
# ------------------
for s in files:
    hists[s]["eff_track"]    = hists[s]["track_eta"].Clone(f"{s}_eff_track")
    hists[s]["eff_cluster"]  = hists[s]["cluster_eta"].Clone(f"{s}_eff_cluster")
    hists[s]["eff_electron"] = hists[s]["electron_eta"].Clone(f"{s}_eff_electron")

    hists[s]["eff_track"].Divide(hists[s]["truth_eta"])
    hists[s]["eff_cluster"].Divide(hists[s]["truth_eta"])
    hists[s]["eff_electron"].Divide(hists[s]["truth_eta"])

# ------------------
# Plotting
# ------------------
def save_overlay(key, name, xlabel, ylabel):
    plotHistograms(
        {s: hists[s][key] for s in hists},
        os.path.join(OUTDIR, name),
        xlabel=xlabel,
        ylabel=ylabel,
        atltext=["Muon Collider", "Simulation, no BIB", "|#eta| < 2.4", "MAIA Detector Concept"]
    )

save_overlay("E_over_p", "E_over_p.png", "E/p", "Entries")
save_overlay("pt_reco_over_true", "pt_reco_over_true.png", "p_{T}^{reco}/p_{T}^{true}", "Entries")
save_overlay(
    "pt_reco_minus_true_over_true",
    "pt_reco_minus_true_over_true.png",
    "(p_{T}^{reco}-p_{T}^{true})/p_{T}^{true}",
    "Entries"
)

save_overlay("eff_track", "track_reco_eff_vs_eta.png", "#eta^{true}", "Efficiency")
save_overlay("eff_cluster", "cluster_reco_eff_vs_eta.png", "#eta^{true}", "Efficiency")
save_overlay("eff_electron", "electron_reco_eff_vs_eta.png", "#eta^{true}", "Efficiency")

# ------------------
# Print efficiency stats
# ------------------
for s in files:
    n_truth    = hists[s]["truth_eta"].Integral()
    n_track    = hists[s]["track_eta"].Integral()
    n_cluster  = hists[s]["cluster_eta"].Integral()
    n_electron = hists[s]["electron_eta"].Integral()

    print("\n===================================")
    print(f"Total truth electrons: {int(n_truth)}")

    if n_truth > 0:
        print(f"Track reconstruction efficiency:    {100.0 * n_track / n_truth:.2f}%")
        print(f"Cluster reconstruction efficiency:  {100.0 * n_cluster / n_truth:.2f}%")
        print(f"Electron reconstruction efficiency: {100.0 * n_electron / n_truth:.2f}%")
    else:
        print("No truth electrons found.")
    print("===================================")

print(f"[done] plots in {OUTDIR}")

