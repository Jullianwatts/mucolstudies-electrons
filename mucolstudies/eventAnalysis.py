import pyLCIO, ROOT, os, plotHelper, glob
from plotHelper import plotEfficiencies

ROOT.gROOT.SetBatch()

#samples = glob.glob("/scratch/jwatts/mucol/electronGun_pT_0_50_reco_0_calibtrpid.slcio")
samples = glob.glob("/scratch/jwatts/mucol/data/reco/electronGun_pT_0_50/electronGun_pT_0_50_reco_7.slcio")

PLOT_DIR = "plots2026"
if not os.path.exists(PLOT_DIR):
    os.makedirs(PLOT_DIR)

DR_CUT  = 0.2
TRK_CUT = 2.0
CLU_CUT = 2.0
MCP_CUT = 10.0

hists = {}

hists["all"] = {
    "mcp_el_eta"  : ROOT.TH1F("mcp_el_eta", "", 25, -2.5, 2.5),
    "trk_el_eta"  : ROOT.TH1F("trk_el_eta", "", 25, -2.5, 2.5),
    "clu_el_eta"  : ROOT.TH1F("clu_el_eta", "", 25, -2.5, 2.5),
    "reco_el_eta" : ROOT.TH1F("reco_el_eta", "", 25, -2.5, 2.5),
}

reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()

for f in samples:

    reader.open(f)

    for event in reader:

        electrons = [
            m for m in event.getCollection("MCParticle")
            if abs(m.getPDG()) == 11
            and m.getGeneratorStatus() == 1
            and m.getEnergy() >= MCP_CUT
        ]

        if not electrons:
            continue

        trk_tlvs = [
            plotHelper.getTrackTLV(t)
            for t in event.getCollection("SelectedTracks")
            if plotHelper.getTrackTLV(t).E() >= TRK_CUT
        ]

        clu_tlvs = [
            plotHelper.getClusterTLV(c)
            for c in event.getCollection("PandoraClusters")
            if plotHelper.getClusterTLV(c).E() >= CLU_CUT
        ]

        for mcp in electrons:

            truth = plotHelper.getTLV(mcp)
            eta = truth.Eta()

            hists["all"]["mcp_el_eta"].Fill(eta)

            best_trk_dr = min(
                [truth.DeltaR(t) for t in trk_tlvs],
                default=999
            )

            best_clu_dr = min(
                [truth.DeltaR(c) for c in clu_tlvs],
                default=999
            )

            has_track = best_trk_dr < DR_CUT
            has_clu   = best_clu_dr < DR_CUT

            if has_track:
                hists["all"]["trk_el_eta"].Fill(eta)

            if has_clu:
                hists["all"]["clu_el_eta"].Fill(eta)

            if has_track and has_clu:
                hists["all"]["reco_el_eta"].Fill(eta)

    reader.close()

display_names = {
    "trk_el"  : "SelectedTracks",
    "clu_el"  : "PandoraClusters",
    "reco_el" : "Track + Cluster"
}

eff_map = {}

for obj in display_names:

    num = hists["all"][obj + "_eta"]
    den = hists["all"]["mcp_el_eta"]

    teff = ROOT.TEfficiency(num, den)
    graph = teff.CreateGraph()

    eff_map[display_names[obj]] = graph

plotEfficiencies(
    eff_map,
    os.path.join(PLOT_DIR, "EFFICIENCY_COMBINED.png"),
    xlabel="#eta",
    ylabel="Efficiency"
)

for obj, label in [
    ("trk_el",  "TRACK"),
    ("clu_el",  "CLUSTER"),
    ("reco_el", "TRACK+CLUSTER")
]:

    m = hists["all"][obj + "_eta"].GetEntries()
    t = hists["all"]["mcp_el_eta"].GetEntries()

    print(f"{label}: matched = {int(m)}, total = {int(t)}, eff = {m/t if t>0 else 0:.3f}")
