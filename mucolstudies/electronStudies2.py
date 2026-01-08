import math
import glob
import ROOT
import pyLCIO

exec(open("./plotHelper.py").read())
ROOT.gROOT.SetBatch()

max_events = -1
samples = glob.glob("/scratch/jwatts/mucol/data/reco/electronGun_pT_0_50/electronGun_pT_0_50_reco_0.slcio")
files = {}
slices = ["0_50", "50_250", "250_1000", "1000_5000"]
for s in slices: files[f"electronGun_pT_{s}"] = []
for s in samples:
    sname = s.split("/")[-1]
    if sname not in files: continue
    files[sname] = glob.glob(f"{s}/*.slcio")
mcp_el_counts = {s: 0 for s in files}

# Set up histograms
hists = {}
for s in files:
    hists[s] = {}
    for obj in ["pfo","pfo_el","pfo_el_match", "mcp", "mcp_el", "trk", "trk_el", "trk_el_match",  "clusters", "clusters_el_match"]:
        for vtype in ["obj", "evt"]:
            for var in variables[vtype]:
                hists[s][obj+"_"+var] = ROOT.TH1F(
                    s+"_"+obj+"_"+var, s,
                    variables[vtype][var]["nbins"],
                    variables[vtype][var]["xmin"],
                    variables[vtype][var]["xmax"]
                )
    hists[s]["cluster_nhits"] = ROOT.TH1F(s+"_cluster_nhits", s, 50, 0, 100)
    hists[s]["cluster_r"] = ROOT.TH1F(s+"_cluster_r", s, 50, 0, 3000)

# Perform matching between two TLVs
def isMatched(tlv1, tlv2):
    if tlv1.DeltaR(tlv2) < 0.1:
        return True
    return False

reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()
reader.setReadCollectionNames(["MCParticle","SiTracks", "AllTracks", "PandoraPFOs", "SeedTracks","PandoraClusters"])

for s in files:
    print("Working on sample", s)
    i = 0

    for f in files[s]:
        if max_events > 0 and i >= max_events: break
        reader.open(f)

        for event in reader:
            if max_events > 0 and i >= max_events: break
            if i%100 == 0: print("\tProcessing event:", i)

            n_trk = 0
            n_matched_trk = 0
            n_pfo = 0
            n_mcp = 0
            n_pfo_el = 0
            n_mcp_el = 0
            n_matched_clusters = 0
            n_matched_el = 0
            n_clusters = 0
            no_mcp_count = 0

            try:
                mcps = event.getCollection("MCParticle")
            except:
                mcps = []
                no_mcp_count += 1
                print("No MCP")

            pfos = event.getCollection("PandoraPFOs")
            trks = event.getCollection("SiTracks")
            clusters = event.getCollection("PandoraClusters")

            mcp_electrons = []
            pfo_electrons = []
            trk_electrons = []
            cluster_electrons = []

            for mcp in mcps:
                if not mcp.getGeneratorStatus() == 1: continue
                if mcp == 0: continue 
                mcp_tlv = getTLV(mcp)
                if mcp_tlv.E() < 20: continue
                if abs(mcp_tlv.Eta())>2.4: continue
                fillObjHists(hists[s], "mcp", mcp_tlv)
                n_mcp += 1

                if abs(mcp.getPDG()) == 11:
                    fillObjHists(hists[s], "mcp_el", mcp_tlv)
                    mcp_electrons.append(mcp_tlv)
                    n_mcp_el += 1

            if n_mcp < 1: continue

            hists[s]["mcp_n"].Fill(n_mcp)
            hists[s]["mcp_el_n"].Fill(n_mcp_el)

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

            for cluster in clusters:
                cluster_position = cluster.getPosition()
                cluster_E = cluster.getEnergy()
                if cluster_E < 10: continue

                cluster_nhits = len(cluster.getCalorimeterHits()) if cluster.getCalorimeterHits() else 0
                cluster_x, cluster_y, cluster_z = cluster_position
                cluster_r = math.sqrt(cluster_x**2 + cluster_y**2)

                cluster_vec = ROOT.TVector3(cluster_x, cluster_y, cluster_z)
                direction = cluster_vec.Unit()
                cluster_tlv = ROOT.TLorentzVector()
                cluster_tlv.SetVect(direction * cluster_E)
                cluster_tlv.SetE(cluster_E)

                fillObjHists(hists[s], "clusters", cluster_tlv)
                cluster_electrons.append(cluster_tlv)
                n_clusters += 1
                hists[s]["clusters_n"].Fill(n_clusters)
                hists[s]["cluster_nhits"].Fill(cluster_nhits)

            hists[s]["cluster_r"].Fill(cluster_r)

            for trk in trks:
                trk_tlv = getTrackTLV(trk, m=0.0005)
                if trk_tlv.E() < 10: continue
                fillObjHists(hists[s], "trk", trk_tlv)
                trk_electrons.append(trk_tlv)
                n_trk += 1

            hists[s]["trk_n"].Fill(n_trk)

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

