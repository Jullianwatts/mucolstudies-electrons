import math, glob, ROOT, pyLCIO, os
exec(open("./plotHelper.py").read())

ROOT.gROOT.SetBatch()

PLOT_DIR = "/scratch/jwatts/mucol/mucolstudies/plots2026"
os.makedirs(PLOT_DIR, exist_ok=True)

files = {
    "electronGun_pT_0_50": glob.glob(
        "/scratch/jwatts/mucol/data/reco/electronGun_pT_0_50/electronGun_pT_0_50_reco_7.slcio"
    )
}

DR_MATCH = 0.2

ENDCAP1_MIN = 0.139
ENDCAP1_MAX = 0.698
ENDCAP2_MIN = 2.443
ENDCAP2_MAX = 3.002

h_ncl   = ROOT.TH1F("ncl",   ";N clusters;Counts", 20, 0, 50)
h_sumE  = ROOT.TH1F("sumE",  ";Summed E [GeV];Counts", 50, 0, 100)
h_trueE = ROOT.TH1F("trueE", ";True E [GeV];Counts", 50, 0, 100)
h_maxE  = ROOT.TH1F("maxE",  ";Max cluster E [GeV];Counts", 50, 0, 100)
h_dr    = ROOT.TH1F("dr",    ";#DeltaR(max cluster, track);Counts", 50, 0, 1)

reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()

for sample in files:
    for f in files[sample]:
        print("Opening", f)
        reader.open(f)

        for event in reader:

            mcps = event.getCollection("MCParticle")
            trks = event.getCollection("SelectedTracks")
            clus = event.getCollection("PandoraClusters")

            electrons = [m for m in mcps if abs(m.getPDG()) == 11 and m.getGeneratorStatus() == 1]
            if not electrons:
                continue

            mcp = electrons[0]
            mcp_tlv = getTLV(mcp)
            trueE = mcp_tlv.E()
            theta = mcp_tlv.Theta()

            # --- ONLY ENDCAP ---
            if not ((ENDCAP1_MIN < theta < ENDCAP1_MAX) or (ENDCAP2_MIN < theta < ENDCAP2_MAX)):
                continue
            trk_tlvs = [getTrackTLV(t) for t in trks]
            #trk_tlvs = [getTrackTLV(t) for t in trks if getTrackTLV(t).Pt() > 2]

            cluster_tlvs = []
            cluster_energies = []

            for c in clus:
                if c.getEnergy() < 2:
                    continue

                tlv = getClusterTLV(c)

                cluster_tlvs.append(tlv)
                cluster_energies.append(tlv.E())

            # --- ONLY EVENTS WITH >1 CLUSTER ---
            if len(cluster_tlvs) <= 1:
                continue

            has_track = False
            has_cluster = False

            for t in trk_tlvs:
                if mcp_tlv.DeltaR(t) < DR_MATCH:
                    has_track = True
                    break

            for c in cluster_tlvs:
                if mcp_tlv.DeltaR(c) < DR_MATCH:
                    has_cluster = True
                    break
            if has_track and has_cluster:
                continue

            sumE = sum(cluster_energies)
            ncl  = len(cluster_energies)

            h_ncl.Fill(ncl)
            h_sumE.Fill(sumE)
            h_trueE.Fill(trueE)

            if ncl > 0:

                max_idx = cluster_energies.index(max(cluster_energies))
                max_cl  = cluster_tlvs[max_idx]

                h_maxE.Fill(max_cl.E())

                min_dr = 999.0

                if trk_tlvs:

                    for t in trk_tlvs:

                        dr = max_cl.DeltaR(t)

                        if dr < min_dr:
                            min_dr = dr

                    h_dr.Fill(min_dr)

                print(f"\nFAILED BASELINE EVENT")
                print(f"trueE = {trueE:.2f} GeV")
                print(f"theta = {theta:.3f}")
                print(f"nClusters = {ncl}")
                print(f"sumClusterE = {sumE:.2f} GeV")
                print(f"maxClusterE = {max_cl.E():.2f} GeV")
                print(f"nTracks = {len(trk_tlvs)}")
                print(f"min dR(max cluster, track) = {min_dr:.3f}")

                print("Cluster energies:", [round(e,2) for e in cluster_energies])

        reader.close()

out = os.path.join(PLOT_DIR, "FAILED_BASELINE_ENDCAP_MULTICLUSTER.png")

plotHistograms(
    {
        "N clusters": h_ncl,
        "Sum E": h_sumE,
        "True E": h_trueE,
        "Max cluster E": h_maxE
    },
    out,
    xlabel="Value",
    ylabel="Counts"
)

print("Saved:", out)
