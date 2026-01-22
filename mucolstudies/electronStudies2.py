import math
import glob
import ROOT
import pyLCIO
import os

exec(open("./plotHelper.py").read())
ROOT.gROOT.SetBatch()

# PLOT OUTPUT DIRECTORY
PLOT_DIR = "/scratch/jwatts/mucol/mucolstudies/plots2026"
os.makedirs(PLOT_DIR, exist_ok=True)

max_events = -1
samples = glob.glob("/scratch/jwatts/mucol/data/reco/electronGun_pT_0_50/electronGun_pT_0_50_reco_4.slcio")

files = {}
slices = ["0_50", "50_250", "250_1000", "1000_5000"]

for sl in slices: files[f"electronGun_pT_{sl}"] = []

for s in samples:
    if s.endswith(".slcio"):
        for sl in slices:
            if f"pT_{sl}" in s:
                files[f"electronGun_pT_{sl}"].append(s)
    else:
        for sl in slices:
            if f"pT_{sl}" in s: files[f"electronGun_pT_{sl}"] = glob.glob(f"{s}/*.slcio")

mcp_el_counts = {s: 0 for s in files}

label_map = {"electronGun_pT_0_50": "0-50 GeV", "electronGun_pT_50_250": "50-250 GeV", "electronGun_pT_250_1000": "250-1000 GeV", "electronGun_pT_1000_5000": "1000- 5000 GeV"}

# HISTOGRAM SETUP
hists = {}
for s in files:
    hists[s] = {}
    for obj in ["pfo","pfo_el","pfo_el_match", "mcp", "mcp_el", "trk", "trk_el", "trk_el_match", "clusters", "clusters_el_match"]:
        for vtype in ["obj", "evt"]:
            for var in variables[vtype]:
                hists[s][obj+"_"+var] = ROOT.TH1F(s+"_"+obj+"_"+var,s,
                    variables[vtype][var]["nbins"],
                    variables[vtype][var]["xmin"],
                    variables[vtype][var]["xmax"])
        hists[s]["cluster_nhits"] = ROOT.TH1F(s+"_cluster_nhits", s, 100, 0, 100)
        hists[s]["cluster_r"] = ROOT.TH1F(s+"_cluster_r", s, 100, 0, 2000)

# MATCHING FUNCTION
def isMatched(tlv1, tlv2): return tlv1.DeltaR(tlv2) < 0.1

reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()
reader.setReadCollectionNames(["MCParticle","SiTracks","AllTracks","PandoraPFOs","SeedTracks","PandoraClusters"])

# EVENT LOOP
for s in files: print("Working on sample", s)
    i = 0

    for f in files[s]:
        if max_events > 0 and i >= max_events: break

        reader.open(f)

        for event in reader:
            if max_events > 0 and i >= max_events: break
            if i % 100 == 0: print("\tProcessing event:", i)

            n_trk = 0
            n_matched_trk = 0
            n_pfo = 0
            n_mcp = 0
            n_pfo_el = 0
            n_mcp_el = 0
            n_matched_clusters = 0
            n_matched_el = 0
            n_clusters = 0

            try: mcps = event.getCollection("MCParticle")
            except: mcps = []

            pfos     = event.getCollection("PandoraPFOs")
            trks     = event.getCollection("SiTracks")
            clusters = event.getCollection("PandoraClusters")

            mcp_electrons     = []
            pfo_electrons     = []
            trk_electrons     = []
            cluster_electrons = []

            # MCP LOOP
            for mcp in mcps:
                if mcp.getGeneratorStatus() != 1: continue
                mcp_tlv = getTLV(mcp)
                if mcp_tlv.E() < 20: continue
                if abs(mcp_tlv.Eta()) > 2.4: continue

                fillObjHists(hists[s], "mcp", mcp_tlv)
                n_mcp += 1

                if abs(mcp.getPDG()) == 11:
                    fillObjHists(hists[s], "mcp_el", mcp_tlv)
                    mcp_electrons.append(mcp_tlv)
                    n_mcp_el += 1

            if n_mcp < 1: continue

            hists[s]["mcp_n"].Fill(n_mcp)
            hists[s]["mcp_el_n"].Fill(n_mcp_el)

            # PFO LOOP
            for pfo in pfos:
                pfo_tlv = getTLV(pfo)
                if pfo_tlv.E() < 10: continue
                fillObjHists(hists[s], "pfo", pfo_tlv)
                n_pfo += 1

                if abs(pfo.getType()) == 11:
                    fillObjHists(hists[s], "pfo_el", pfo_tlv)
                    pfo_electrons.append(pfo_tlv)
                    n_pfo_el += 1

            hists[s]["pfo_n"].Fill(n_pfo)
            hists[s]["pfo_el_n"].Fill(n_pfo_el)

            # CLUSTER LOOP
            for cluster in clusters:
                cluster_E = cluster.getEnergy()
                if cluster_E < 10: continue

                cluster_position = cluster.getPosition()
                cluster_x = cluster_position[0]
                cluster_y = cluster_position[1]
                cluster_z = cluster_position[2]

                cluster_r = math.sqrt(cluster_x**2 + cluster_y**2)

                cluster_vec = ROOT.TVector3(cluster_x, cluster_y, cluster_z)
                direction = cluster_vec.Unit()

                cluster_tlv = ROOT.TLorentzVector()
                cluster_tlv.SetVect(direction * cluster_E)
                cluster_tlv.SetE(cluster_E)

                fillObjHists(hists[s], "clusters", cluster_tlv)
                cluster_electrons.append(cluster_tlv)
                n_clusters += 1
                if "cluster_r" in hists[s]:
                   hists[s]["cluster_r"].Fill(cluster_r)

                hists[s]["clusters_n"].Fill(n_clusters)
                hists[s]["cluster_nhits"].Fill(
                len(cluster.getCalorimeterHits()) if cluster.getCalorimeterHits() else 0)

            # TRACK LOOP
            for trk in trks:
                trk_tlv = getTrackTLV(trk, m=0.0005)
                if trk_tlv.E() < 10: continue

                fillObjHists(hists[s], "trk", trk_tlv)
                trk_electrons.append(trk_tlv)
                n_trk += 1

            hists[s]["trk_n"].Fill(n_trk)

            # MATCHING
            for mcp_el in mcp_electrons:
                best_track = None
                best_track_dr = 20

                for track in trk_electrons:
                    dr = mcp_el.DeltaR(track)
                    if dr < 0.1 and dr < best_track_dr:
                        best_track = track
                        best_track_dr = dr

                if best_track:
                    fillObjHists(hists[s], "trk_el_match", mcp_el)
                    n_matched_trk += 1

                best_pfo = None
                best_pfo_dr = 20
                for pfo in pfo_electrons:
                    dr = mcp_el.DeltaR(pfo)
                    if dr < 0.1 and dr < best_pfo_dr:
                        best_pfo = pfo
                        best_pfo_dr = dr

                if best_pfo:
                    fillObjHists(hists[s], "pfo_el_match", mcp_el)
                    n_matched_el += 1

                best_cluster = None
                best_cluster_dr = 20
                for cluster in cluster_electrons:
                    dr = mcp_el.DeltaR(cluster)
                    if dr < 0.1 and dr < best_cluster_dr:
                        best_cluster = cluster
                        best_cluster_dr = dr

                if best_cluster:
                    fillObjHists(hists[s], "clusters_el_match", mcp_el)
                    n_matched_clusters += 1

            hists[s]["trk_el_match_n"].Fill(n_matched_trk)
            hists[s]["clusters_el_match_n"].Fill(n_matched_clusters)

            i += 1

        reader.close()

print(f"\nSaving plots to {PLOT_DIR}...")

for s in hists:
    denom = hists[s]["mcp_el_eta"]
    eff_map = {}
    
    for obj in ["trk_el_match", "pfo_el_match", "clusters_el_match"]:
        num = hists[s][obj + "_eta"]
        
        teff = ROOT.TEfficiency(num, denom)
        
        eff_graph = teff.CreateGraph()
        eff_map[obj] = eff_graph
        
        single_map = {obj: eff_graph}
        save_path = os.path.join(PLOT_DIR, f"EFFICIENCY_{obj}_{s}.png")
        plotEfficiencies(single_map, save_path, xlabel="#eta", ylabel="Efficiency")

    combined_path = os.path.join(PLOT_DIR, f"EFFICIENCY_COMBINED_{s}.png")
    plotEfficiencies(eff_map, combined_path, xlabel="#eta", ylabel="Efficiency")

print("TRACK MATCHING EFFICIENCY SUMMARY")
for s in hists:
    matched = hists[s]["trk_el_match_eta"].GetEntries() if "trk_el_match_eta" in hists[s] else 0
    total   = hists[s]["mcp_el_eta"].GetEntries() if "mcp_el_eta" in hists[s] else 0
    eff     = matched / total if total > 0 else 0
    print(f"{label_map.get(s, s)}: matched = {int(matched)}, total = {int(total)}, eff = {eff:.3f}")

print(" PFO MATCHING EFFICIENCY SUMMARY")
for s in hists:
    matched = hists[s]["pfo_el_match_eta"].GetEntries() if "pfo_el_match_eta" in hists[s] else 0
    total   = hists[s]["mcp_el_eta"].GetEntries() if "mcp_el_eta" in hists[s] else 0
    eff     = matched / total if total > 0 else 0
    print(f"{label_map.get(s, s)}: matched = {int(matched)}, total = {int(total)}, eff = {eff:.3f}")

print("CLUSTER MATCHING EFFICIENCY")
for s in hists:
    matched = hists[s]["clusters_el_match_eta"].GetEntries() if "clusters_el_match_eta" in hists[s] else 0
    total   = hists[s]["mcp_el_eta"].GetEntries() if "mcp_el_eta" in hists[s] else 0
    eff     = matched / total if total > 0 else 0
    print(f"{label_map.get(s, s)}: matched = {int(matched)}, total = {int(total)}, eff = {eff:.3f}")

