# Note: This script was to test out opening k4hep files with ROOT. An actual example can be found here:
# https://gitlab.desy.de/ftx-sft-key4hep/tutorials/-/blob/main/edm4hep_analysis/edm4hep_python.ipynb

# Before running this script be sure to run the commands needed to access the software:
# apptainer build k4toroid.sif docker://madbaron/k4test-ubuntu:latest
# apptainer run --no-home -B /collab/project/snowmass21/data/muonc:/data -B /home/$USER k4toroid.sif
# source /setup.sh
import math
import glob
import ROOT
#from ROOT import edm4hep
#from podio.root_io import Reader
import pyLCIO

exec(open("./plotHelper.py").read())
ROOT.gROOT.SetBatch()

# Set up some options
max_events = -1
# Open the edm4hep files with ROOT
#samples = glob.glob("/data/fmeloni/DataMuC_MuColl10_v0A/k4reco/electronGun*")
#samples = glob.glob("/data/fmeloni/DataMuC_MuColl10_v0A/reco/electronGun*")
#samples = glob.glob("/data/fmeloni/DataMuC_MuColl10_v0A/recoBIB/electronGun*")
#samples = glob.glob("/data/fmeloni/DataMuC_MuColl10_v0A/reco_highrange/electronGun*")
#samples = glob.glob("/data/fmeloni/DataMuC_MAIA_v0/v3/reco/electronGun*")
#samples = glob.glob("/data/fmeloni/DataMuC_MAIA_v0/v3/electronGun*")
#samples = glob.glob("/data/fmeloni/DataMuC_MAIA_v0/v4rotated/electronGun*")
#samples = glob.glob("/data/fmeloni/DataMuC_MAIA_v0/v4/electronGun*")
#samples = glob.glob("/data/fmeloni/DataMuC_MuColl10_v0A/v0/reco/electronGun*")
samples = glob.glob("/data/fmeloni/DataMuC_MAIA_v0/v5/reco/electronGun*")
#samples = glob.glob("/data/fmeloni/DataMuC_MAIA_v0/v2/reco/electronGun*")
files = {}
print("=== FILE DEBUG ===")
for s in files:
    print(f"Sample {s}: {len(files[s])} files found")
    if len(files[s]) > 0:
        print(f"  First file: {files[s][0]}")
print("==================")
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
                hists[s][obj+"_"+var] = ROOT.TH1F(s+"_"+obj+"_"+var, s, variables[vtype][var]["nbins"], variables[vtype][var]["xmin"], variables[vtype][var]["xmax"])
    hists[s]["cluster_nhits"] = ROOT.TH1F(s+"_cluster_nhits", s, 50, 0, 100)
    hists[s]["cluster_r"] = ROOT.TH1F(s+"_cluster_r", s, 50, 0, 3000)

hists2d = {}
for s in files:
    hists2d[s] = {}

    hists2d[s]["trk_eta_v_trk_pt"] = ROOT.TH2F(f"trk_eta_v_trk_pt_{s}", f"trk_eta_v_trk_pt_{s}", 30,-3,3,30,0,3000)
    hists2d[s]["trk_eta_v_trk_phi"] = ROOT.TH2F(f"trk_eta_v_trk_phi_{s}", f"trk_eta_v_trk_phi_{s}", 30,-3,3,30,-3,3)
    hists2d[s]["trk_eta_v_trk_n"] = ROOT.TH2F(f"trk_eta_v_trk_n_{s}", f"trk_eta_v_trk_n_{s}", 30,-3,3,20,0,20)
    hists2d[s]["trk_eta_v_mcp_eta"] = ROOT.TH2F(f"trk_eta_v_mcp_eta_{s}", f"trk_eta_v_mcp_eta_{s}", 30,-3,3,30,-3,3)
    hists2d[s]["trk_eta_v_mcp_pt"] = ROOT.TH2F(f"trk_eta_v_mcp_pt_{s}", f"trk_eta_v_mcp_pt_{s}", 30,-3,3,30,0,3000)
    hists2d[s]["trk_eta_v_mcp_phi"] = ROOT.TH2F(f"trk_eta_v_mcp_phi_{s}", f"trk_eta_v_mcp_phi_{s}", 30,-3,3,30,-3,3)
    hists2d[s]["trk_pt_v_mcp_pt"] = ROOT.TH2F(f"trk_pt_v_mcp_pt_{s}", f"trk_pt_v_mcp_pt_{s}", 30,0,3000,30,0,3000)
    hists2d[s]["pfo_pt_v_mcp_pt"] = ROOT.TH2F(f"pfo_pt_v_mcp_pt_{s}", f"pfo_pt_v_mcp_pt_{s}", 30,0,3000,30,0,3000)
    hists2d[s]["mcp_E_v_mcp_p"] = ROOT.TH2F(f"mcp_E_v_mcp_p_{s}", f"mcp_E_v_mcp_p_{s}", 30,0,1000,30,0,1000)

# Perform matching between two TLVs
def isMatched(tlv1, tlv2):
    if tlv1.DeltaR(tlv2) < 0.1: #and abs(tlv1.Perp()-tlv2.Perp())/tlv2.Perp() < 0.2:       
        return True
    return False

# Create a reader object to use for the rest of the time
reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()
reader.setReadCollectionNames(["MCParticle","SiTracks", "AllTracks", "PandoraPFOs", "SeedTracks","PandoraClusters"]) #changed from SiTracks_refitted for new path
# Loop over the different samples
for s in files:
    print("Working on sample", s)
    i = 0

    # Loop over the files in a sample
    for f in files[s]:
        if max_events > 0 and i >= max_events: break

        reader.open(f)

        # Loop over events in each file
        #for event in reader.get('events'):
        for event in reader:
            if max_events > 0 and i >= max_events: break
            if i%100 == 0: print("\tProcessing event:", i)

            # Make counters and key objects to track
            n_trk = 0
            n_matched_trk = 0
            n_pfo = 0
            n_mcp = 0
            n_pfo_el = 0
            n_mcp_el = 0
            n_matched_clusters = 0
            n_matched_el = 0
            n_clusters = 0
            #my_mcp_el = None
            no_mcp_count = 0
            # Get the collections we care about
            #mcps = event.get("MCParticle")
            #pfos = event.get("PandoraPFOs")
            #trks = event.get("SiTracks_Refitted")
            #lnks = event.get("MCParticle_SiTracks_Refitted")
            try:
                mcps = event.getCollection("MCParticle")
            except:
                mcps = []
                no_mcp_count += 1
                print("No MCP")
            pfos = event.getCollection("PandoraPFOs")
            trks = event.getCollection("SiTracks") #also changed from SiTracks_refitted
            clusters = event.getCollection("PandoraClusters") 
            mcp_electrons = []
            pfo_electrons = []
            trk_electrons = []
            cluster_electrons = []
######## Loop over MCPs
            for mcp in mcps:
                if not mcp.getGeneratorStatus() == 1: continue
                if mcp == 0: continue 
                mcp_tlv = getTLV(mcp)
                if mcp_tlv.E() < 20: continue
                if abs(mcp_tlv.Eta())>2.4: continue
                fillObjHists(hists[s], "mcp", mcp_tlv)
                momentum = math.sqrt(mcp.getMomentum()[0]**2+mcp.getMomentum()[1]**2+mcp.getMomentum()[2]**2)
                hists2d[s]["mcp_E_v_mcp_p"].Fill(mcp.getEnergy(), momentum)
                n_mcp += 1
                # Look at electrons only
                if abs(mcp.getPDG()) == 11:
                    fillObjHists(hists[s], "mcp_el", mcp_tlv)
                    mcp_electrons.append(mcp_tlv)
                    n_mcp_el += 1

            # Only consider events that had at least one electron in our fiducial region
            if n_mcp < 1: continue

            hists[s]["mcp_n"].Fill(n_mcp)
            hists[s]["mcp_el_n"].Fill(n_mcp_el)
    
 
            for pfo in pfos:
                pfo_tlv = getTLV(pfo)
                if pfo_tlv.E() < 10: continue
                fillObjHists(hists[s], "pfo", pfo_tlv)
                n_pfo += 1

                #hists2d[s]["pfo_pt_v_mcp_pt"].Fill(pfo_tlv.Perp(), my_mcp_el.Perp())

                # Look at electrons only
                if abs(pfo.getType()) == 11:
                    fillObjHists(hists[s], "pfo_el", pfo_tlv)
                    pfo_electrons.append(pfo_tlv)
                    n_pfo_el += 1
            hists[s]["pfo_n"].Fill(n_pfo)
            hists[s]["pfo_el_n"].Fill(n_pfo_el)


            ### loop over clusters
            
            for cluster in clusters:
                cluster_position = cluster.getPosition()
                cluster_E = cluster.getEnergy()
                if cluster_E < 10:
                    continue
                cluster_nhits = len(cluster.getCalorimeterHits()) if cluster.getCalorimeterHits() else 0
                cluster_x, cluster_y, cluster_z = cluster_position[0], cluster_position[1], cluster_position[2]
                cluster_r = math.sqrt(cluster_x**2 + cluster_y**2)
                cluster_eta = math.asinh(cluster_z / cluster_r) if cluster_r > 0 else 0
                if i < 3:  # Only for first few events
                    print(f"Cluster: E={cluster_E:.2f}, nhits={cluster_nhits}, eta={cluster_eta:.2f}")
 
                cluster_vec = ROOT.TVector3()
                cluster_vec= ROOT.TVector3(cluster_position[0], cluster_position[1], cluster_position[2])
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

            ######## Loop over tracks
            for trk in trks:
                trk_tlv = getTrackTLV(trk, m=0.0005)
                if trk_tlv.E() < 10: continue
                fillObjHists(hists[s], "trk", trk_tlv)

                hists2d[s]["trk_eta_v_trk_pt"].Fill(trk_tlv.Eta(), trk_tlv.Perp())
                hists2d[s]["trk_eta_v_trk_phi"].Fill(trk_tlv.Eta(), trk_tlv.Phi())
                hists2d[s]["trk_eta_v_trk_n"].Fill(trk_tlv.Eta(), len(trks))
                #hists2d[s]["trk_eta_v_mcp_eta"].Fill(trk_tlv.Eta(), my_mcp_el.Eta())
                #hists2d[s]["trk_eta_v_mcp_pt"].Fill(trk_tlv.Eta(), my_mcp_el.Perp())
                #hists2d[s]["trk_eta_v_mcp_phi"].Fill(trk_tlv.Eta(), my_mcp_el.Phi())
                #hists2d[s]["trk_pt_v_mcp_pt"].Fill(trk_tlv.Perp(), my_mcp_el.Perp())
                trk_electrons.append(trk_tlv)
                n_trk += 1
                
            hists[s]["trk_n"].Fill(n_trk)
            
            #one-to-one matching for each MCP electron
            for mcp_el in mcp_electrons:
                # Find best track match
                best_track = None
                best_track_dr = 20
                for track in trk_electrons:
                    dr = mcp_el.DeltaR(track)
                    if dr < 0.1 and dr < best_track_dr:
                        best_track = track
                        best_track_dr = dr
                
                if best_track:
                    fillObjHists(hists[s], "trk_el_match", mcp_el)  # Fill with MCP!
                    n_matched_trk += 1
                
                # Find best PFO match
                best_pfo = None
                best_pfo_dr = 20
                for pfo in pfo_electrons:
                    dr = mcp_el.DeltaR(pfo)
                    if dr < 0.1 and dr < best_pfo_dr:
                        best_pfo = pfo
                        best_pfo_dr = dr
                
                if best_pfo:
                    fillObjHists(hists[s], "pfo_el_match", mcp_el)  # Fill with MCP!
                    n_matched_el += 1
                
                # Find best cluster match
                best_cluster = None
                best_cluster_dr = 20
                for cluster in cluster_electrons:
                    dr = mcp_el.DeltaR(cluster)
                    if dr < 0.1 and dr < best_cluster_dr:
                        best_cluster = cluster
                        best_cluster_dr = dr
                
                if best_cluster:
                    fillObjHists(hists[s], "clusters_el_match", mcp_el)  # Fill with MCP!
                    n_matched_clusters += 1

            hists[s]["trk_el_match_n"].Fill(n_matched_trk)
            hists[s]["clusters_el_match_n"].Fill(n_matched_clusters)

            i += 1
            #if "1000_5000" in s:
                #mcp_el_counts[s] += n_mcp_el

        reader.close()
#total = sum(mcp_el_counts[s] for s in mcp_el_counts if "1000_5000" in s)
#print(f"\nTotal nymber of MC electrons in 1000-5000 GeV slice: {total}")

from plotHelper import plotHistograms

# Helper for display names
label_map = {
    "electronGun_pT_0_50": "0-50 GeV",
    "electronGun_pT_50_250": "50-250 GeV",
    "electronGun_pT_250_1000": "250-1000 GeV",
    "electronGun_pT_1000_5000": "1000-5000 GeV"
}

# TRACK EFFICIENCY CALCULATION
print("\n=== CALCULATING TRACK EFFICIENCY ===")
eff_eta_trk = {}
for s in hists:
    if hists[s]["mcp_el_eta"].GetEntries() == 0:
        continue
    if hists[s]["trk_el_match_eta"].GetEntries() == 0:
        continue
    
    num = hists[s]["trk_el_match_eta"]
    denom = hists[s]["mcp_el_eta"]
    eff = ROOT.TGraphAsymmErrors()
    eff.BayesDivide(num, denom, "B")
    eff_eta_trk[label_map[s]] = eff

plotEfficiencies(eff_eta_trk,"plots/track_efficiency_as_function_of_eta_allSlices.png",
        xlabel="eta",ylabel="Track Matching Efficiency"
)

print("\nTrack Matching Efficiency Summary:")
for s in hists:
    matched = hists[s]["trk_el_match_eta"].GetEntries() if "trk_el_match_eta" in hists[s] else 0
    total = hists[s]["mcp_el_eta"].GetEntries() if "mcp_el_eta" in hists[s] else 0
    efficiency = matched/total if total > 0 else 0
    print(f"{label_map.get(s, s)}: matched = {int(matched)}, total = {int(total)}, eff = {efficiency:.3f}")

# PFO EFFICIENCY CALCULATION
print("\n=== CALCULATING PFO EFFICIENCY ===")
eff_eta_pfo = {}
for s in hists:
    if "pfo_el_match_eta" not in hists[s] or hists[s]["pfo_el_match_eta"].GetEntries() == 0:
        continue
    if hists[s]["mcp_el_eta"].GetEntries() == 0:
        continue
    num = hists[s]["pfo_el_match_eta"]
    denom = hists[s]["mcp_el_eta"]
    eff = ROOT.TGraphAsymmErrors()
    eff.BayesDivide(num, denom, "B")
    eff_eta_pfo[label_map[s]] = eff

plotEfficiencies(eff_eta_pfo,"plots/pfo_efficiency_as_function_of_eta_allSlices.png",
        xlabel="eta", ylabel="PFO Matching Efficiency"
        )
# PFO EFFICIENCY vs ENERGY
print("\n=== CALCULATING PFO EFFICIENCY vs ENERGY ===")
eff_energy_pfo = {}
for s in hists:
    if "pfo_el_match_E" not in hists[s] or hists[s]["pfo_el_match_E"].GetEntries() == 0:
        continue
    if hists[s]["mcp_el_E"].GetEntries() == 0:
        continue
    num = hists[s]["pfo_el_match_E"]
    denom = hists[s]["mcp_el_E"]
    eff = ROOT.TGraphAsymmErrors()
    eff.BayesDivide(num, denom, "B")
    eff_energy_pfo[label_map[s]] = eff

plotEfficiencies(eff_energy_pfo, "plots/pfo_efficiency_as_function_of_energy_allSlices.png",
    xlabel="Energy [GeV]", ylabel="PFO Matching Efficiency"
)

print("\nPFO Matching Efficiency vs Energy Summary:")
for s in hists:
    matched = hists[s]["pfo_el_match_E"].GetEntries() if "pfo_el_match_E" in hists[s] else 0
    total = hists[s]["mcp_el_E"].GetEntries() if "mcp_el_E" in hists[s] else 0
    efficiency = matched/total if total > 0 else 0
    print(f"{label_map.get(s, s)}: matched = {int(matched)}, total = {int(total)}, eff = {efficiency:.3f}")
print("\nPFO Matching Efficiency Summary:")
for s in hists:
    matched = hists[s]["pfo_el_match_eta"].GetEntries() if "pfo_el_match_eta" in hists[s] else 0
    total = hists[s]["mcp_el_eta"].GetEntries() if "mcp_el_eta" in hists[s] else 0
    efficiency = matched/total if total > 0 else 0
    print(f"{label_map.get(s, s)}: matched = {int(matched)}, total = {int(total)}, eff = {efficiency:.3f}")

# CLUSTER EFFICIENCY CALCULATION
print("\n=== CALCULATING CLUSTER EFFICIENCY ===")
eff_eta_clus = {}
for s in hists:
    if "clusters_el_match_eta" not in hists[s] or hists[s]["clusters_el_match_eta"].GetEntries() == 0:
        print(f"Skipping {s}: no matched clusters")
        continue
    if hists[s]["mcp_el_eta"].GetEntries() == 0:
        print(f"Skipping {s}: no MCP electrons")
        continue
    
    num = hists[s]["clusters_el_match_eta"]
    denom = hists[s]["mcp_el_eta"]
    eff = ROOT.TGraphAsymmErrors()
    eff.BayesDivide(num, denom, "B")
    eff_eta_clus[label_map[s]] = eff
    
plotEfficiencies(eff_eta_clus,"plots/cluster_efficiency_as_function_of_eta_allSlices.png",
        xlabel="eta",ylabel="Cluster Matching Efficiency"
    )

print("\nCluster Matching Efficiency Summary:")
for s in hists:
    matched = hists[s]["clusters_el_match_eta"].GetEntries() if "clusters_el_match_eta" in hists[s] else 0
    total = hists[s]["mcp_el_eta"].GetEntries() if "mcp_el_eta" in hists[s] else 0
    efficiency = matched/total if total > 0 else 0
    print(f"{label_map.get(s, s)}: matched = {int(matched)}, total = {int(total)}, eff = {efficiency:.3f}")


print(no_mcp_count)
for i, h in enumerate(hists[s]):
    if h in ["cluster_nhits", "cluster_r"]:
        continue

    # Collect hists that go on a single plot
    hists_to_plot = {}
    for j, s in enumerate(hists):
        hists_to_plot[s] = hists[s][h]
    var_name = h.split("_")[-1]
    try:
        xlabel = variables["obj"][var_name]["label"]
    except:
        xlabel = variables["evt"][var_name]["label"]

    # Call plotting function
    plotHistograms(hists_to_plot, "plots/"+h+".png", xlabel, "Entries")
    plotHistograms(hists_to_plot, "plots/"+h+".png", xlabel, "Entries")
    #plotHistograms(hists_to_plot, "plots/electrons/"+h+".png", xlabel, "Entries")
    #plotHistograms(hists_to_plot, "plots/electrons/"+h+".root", xlabel, "Entries")
    #plotHistograms(hists_to_plot, "plots/electrons_no_el_req/"+h+".png", xlabel, "Entries")
    #plotHistograms(hists_to_plot, "plots/electrons_no_el_req/"+h+".root", xlabel, "Entries")
cluster_hists_nhits = {}
cluster_hists_r = {}

for s in files:
    if "cluster_nhits" in hists[s]:
        cluster_hists_nhits[s] = hists[s]["cluster_nhits"]
    if "cluster_r" in hists[s]:
        cluster_hists_r[s] = hists[s]["cluster_r"]

if cluster_hists_nhits:
    plotHistograms(cluster_hists_nhits, "plots/cluster_nhits.png", "Number of Hits per Cluster", "Entries")
if cluster_hists_r:
    plotHistograms(cluster_hists_r, "plots/cluster_r.png", "Cluster Radial Position [mm]", "Entries")
for s in hists2d:
    for h in hists2d[s]:
        c = ROOT.TCanvas("can", "can")
        hists2d[s][h].Draw("colz")
        hists2d[s][h].GetXaxis().SetTitle(h.split("_v_")[0])
        hists2d[s][h].GetYaxis().SetTitle(h.split("_v_")[1])
        #c.SaveAs(f"plots/electrons/{hists2d[s][h].GetName()}.png")
        c.SaveAs(f"plots/{hists2d[s][h].GetName()}:.png")
