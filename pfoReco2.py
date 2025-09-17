import os
import ROOT
import glob
import pyLCIO
from pyLCIO import IOIMPL

exec(open("./plotHelper.py").read())
ROOT.gROOT.SetBatch()


max_events = 10

samples = glob.glob("/data/fmeloni/DataMuC_MAIA_v0/v5/reco/pionGun*")

files = {}
slices = ["0_50", "50_250", "250_1000", "1000_5000"]
for s in slices: 
    files[f"pionGun_pT_{s}"] = []

for s in samples:
    sname = s.split("/")[-1]
    if sname not in files: 
        continue
    files[sname] = glob.glob(f"{s}/*.slcio")

match_counts = {s: 0 for s in files}

hists = {}
for s in files:
    hists[s] = {}
    for obj in ["pfo","pfo_pi","pfo_pi_match", "mcp", "mcp_pi", "trk", "trk_pi", "trk_pi_match", "clusters", "clusters_pi_match"]:
        for vtype in ["obj", "evt"]:
            for var in variables[vtype]:
                hists[s][obj+"_"+var] = ROOT.TH1F(s+"_"+obj+"_"+var, s, 
                                                variables[vtype][var]["nbins"], 
                                                variables[vtype][var]["xmin"], 
                                                variables[vtype][var]["xmax"])
    
    # E/p ratio hists for pions
    hists[s]["pfo_pi_match_obj_ep_ratio"] = ROOT.TH1F(s+"_pfo_pi_match_obj_ep_ratio", s, 100, 0, 3)
    hists[s]["clusters_pi_match_obj_ep_ratio"] = ROOT.TH1F(s+"_clusters_pi_match_obj_ep_ratio", s, 100, 0, 3)

def find_closest_track(cluster_tlv, track_list, dR_cut=0.1):
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
            if i%100 == 0: 
                print(f"\tProcessing event: {i}")
            
            try:
                mcps = event.getCollection("MCParticle")
            except:
                mcps = []
                print("No MCP")
            
            try:
                clusters = event.getCollection("PandoraClusters")
                tracks = event.getCollection("SiTracks")
            except:
                i += 1
                continue
            
            # Convert to TLorentzVectors
            cluster_tlvs = []
            track_tlvs = []
            
            for cluster in clusters:
                cluster_tlv = getClusterTLV(cluster)
                if cluster_tlv.E() < 10: 
                    continue
                cluster_tlvs.append(cluster_tlv)
            
            for track in tracks:
                # pion mass for TLV
                track_tlv = getTrackTLV(track, m=0.13957)
                if track_tlv.E() < 10: 
                    continue
                track_tlvs.append(track_tlv)
            
            # Print first 100 events
            if i < 100:
                print(f"\n--- {s} Event {i} ---")
                print(f"Clusters ({len(cluster_tlvs)} total):")
                for j, cluster_tlv in enumerate(cluster_tlvs):
                    print(f"  Cluster {j}: E = {cluster_tlv.E():.2f} GeV")
                print(f"Tracks ({len(track_tlvs)} total):")
                for j, track_tlv in enumerate(track_tlvs):
                    print(f"  Track {j}: p = {track_tlv.P():.2f} GeV")
            
            # Process each cluster for E/p
            for cluster_tlv in cluster_tlvs:
                closest_track = find_closest_track(cluster_tlv, track_tlvs)
                if closest_track is None:
                    continue
                
                match_counts[s] += 1
                
                E = cluster_tlv.E()
                p = closest_track.P()
                if p <= 0: 
                    continue
                ep_ratio = E / p
                
                # fill pion hists
                fillObjHists(hists[s], "pfo_pi_match", closest_track)
                fillObjHists(hists[s], "clusters_pi_match", cluster_tlv)
                hists[s]["pfo_pi_match_obj_ep_ratio"].Fill(ep_ratio)
                hists[s]["clusters_pi_match_obj_ep_ratio"].Fill(ep_ratio)
                
                if i < 10:
                    print(f"{s}: E={E:.1f} GeV, p={p:.1f} GeV, E/p={ep_ratio:.3f}")
            
            i += 1
        
        reader.close()

print(f"Total matches: {sum(match_counts.values())}")

# Labels keyed by pion samples
label_map = {
    "pionGun_pT_0_50": "0-50 GeV",
    "pionGun_pT_50_250": "50-250 GeV", 
    "pionGun_pT_250_1000": "250-1000 GeV",
    "pionGun_pT_1000_5000": "1000-5000 GeV"
}

# Plot E/p ratios comparison
ep_hists = {}
for sample_name in files:
    if "pfo_pi_match_obj_ep_ratio" in hists[sample_name] and hists[sample_name]["pfo_pi_match_obj_ep_ratio"].GetEntries() > 0:
        ep_hists[sample_name] = hists[sample_name]["pfo_pi_match_obj_ep_ratio"]

if ep_hists:
    plotHistograms(ep_hists, "plots/pi_ep_ratio_comparison.png", 
                  xlabel="E/p", ylabel="Entries",
                  atltext=["Muon Collider", "Simulation, no BIB", "|#eta| < 2.4", "MAIA Detector Concept"])

# Plot variables
for var in ["pt", "eta", "phi", "E"]:
    var_hists = {}
    for sample_name in files:
        hist_name = f"pfo_pi_match_obj_{var}"
        if hist_name in hists[sample_name] and hists[sample_name][hist_name].GetEntries() > 0:
            var_hists[sample_name] = hists[sample_name][hist_name]
    
    if var_hists:
        plotHistograms(var_hists, f"plots/pfo_pi_match_{var}.png",
                      xlabel=variables["obj"][var]["label"], ylabel="Entries",
                      atltext=["Muon Collider", "Simulation, no BIB", "|#eta| < 2.4", "MAIA Detector Concept"])

# Plot cluster E/p ratios
cluster_ep_hists = {}
for sample_name in files:
    hist_name = "clusters_pi_match_obj_ep_ratio"
    if hist_name in hists[sample_name] and hists[sample_name][hist_name].GetEntries() > 0:
        cluster_ep_hists[sample_name] = hists[sample_name][hist_name]

if cluster_ep_hists:
    plotHistograms(cluster_ep_hists, "plots/pi_cluster_ep_ratio_comparison.png",
                  xlabel="E/p", ylabel="Entries", 
                  atltext=["Muon Collider", "Simulation, no BIB", "|#eta| < 2.4", "MAIA Detector Concept"])

# Print stats
for sample_name in files:
    if sample_name not in hists or match_counts[sample_name] == 0:
        continue
        
    print(f"\n{sample_name}:")
    print(f"  Total matches: {match_counts[sample_name]}")
    
    ep_hist = hists[sample_name]["pfo_pi_match_obj_ep_ratio"]
    if ep_hist.GetEntries() > 0:
        total_entries = ep_hist.GetEntries()
        near_one = 0
        for i in range(1, ep_hist.GetNbinsX() + 1):
            bin_center = ep_hist.GetBinCenter(i)
            if 0.8 <= bin_center <= 1.2:
                near_one += ep_hist.GetBinContent(i)
        
        percent_near_one = (near_one / total_entries) * 100 if total_entries > 0 else 0
        print(f"  E/p entries with 0.8 < E/p < 1.2: {percent_near_one:.1f}%")

