import ROOT
import pyLCIO
import glob
import math
from plotHelper import plotHistograms

max_events = 10

ecal_collections = [
    "ECalBarrelCollection", "ECalEndcapCollection"
]

collections_to_read = ["PandoraClusters", "MCParticle"] + ecal_collections

def setup_cluster_histograms():
    hists = {}
    hists2d = {}

    samples = ["electronGun_pT_0_50", "electronGun_pT_50_250", "electronGun_pT_250_1000", "electronGun_pT_1000_5000"]

    for sample in samples:
        hists[sample] = {}
        hists2d[sample] = {}

        # Basic cluster properties
        hists[sample]["cluster_energy"] = ROOT.TH1F(f"cluster_energy_{sample}", f"Cluster Energy - {sample}", 100, 0, 100)
        hists[sample]["clusters_per_event"] = ROOT.TH1F(f"clusters_per_event_{sample}", f"Number of Clusters per Event - {sample}", 10, 0, 10)
        
        # Cluster position from cluster.getPosition()
        hists[sample]["cluster_x"] = ROOT.TH1F(f"cluster_x_{sample}", f"Cluster X Position - {sample}", 100, -2000, 2000)
        hists[sample]["cluster_y"] = ROOT.TH1F(f"cluster_y_{sample}", f"Cluster Y Position - {sample}", 100, -2000, 2000) 
        hists[sample]["cluster_z"] = ROOT.TH1F(f"cluster_z_{sample}", f"Cluster Z Position - {sample}", 100, -3000, 3000)
        hists[sample]["cluster_radius"] = ROOT.TH1F(f"cluster_radius_{sample}", f"Cluster Radius from Origin - {sample}", 100, 0, 2000)
        
        # Cluster shape parameters from cluster.getShape()
        hists[sample]["cluster_shape_0"] = ROOT.TH1F(f"cluster_shape_0_{sample}", f"Cluster Shape Parameter 0 - {sample}", 100, 0, 1000)
        hists[sample]["cluster_shape_1"] = ROOT.TH1F(f"cluster_shape_1_{sample}", f"Cluster Shape Parameter 1 - {sample}", 100, 0, 1000)
        hists[sample]["cluster_shape_2"] = ROOT.TH1F(f"cluster_shape_2_{sample}", f"Cluster Shape Parameter 2 - {sample}", 100, 0, 1000)
        hists[sample]["cluster_shape_3"] = ROOT.TH1F(f"cluster_shape_3_{sample}", f"Cluster Shape Parameter 3 - {sample}", 100, 0, 1000)
        hists[sample]["cluster_shape_4"] = ROOT.TH1F(f"cluster_shape_4_{sample}", f"Cluster Shape Parameter 4 - {sample}", 100, 0, 1000)
        hists[sample]["cluster_shape_5"] = ROOT.TH1F(f"cluster_shape_5_{sample}", f"Cluster Shape Parameter 5 - {sample}", 100, 0, 1000)
        
        # Subdetector energies from cluster.getSubdetectorEnergies()
        hists[sample]["subdet_energy_0"] = ROOT.TH1F(f"subdet_energy_0_{sample}", f"Subdetector 0 Energy - {sample}", 100, 0, 100)
        hists[sample]["subdet_energy_1"] = ROOT.TH1F(f"subdet_energy_1_{sample}", f"Subdetector 1 Energy - {sample}", 100, 0, 100)
        hists[sample]["subdet_energy_2"] = ROOT.TH1F(f"subdet_energy_2_{sample}", f"Subdetector 2 Energy - {sample}", 100, 0, 100)
        
        # Cluster type and particle ID related
        hists[sample]["cluster_type"] = ROOT.TH1F(f"cluster_type_{sample}", f"Cluster Type - {sample}", 10, 0, 10)
        
        # Event-level: total energy in clusters vs total energy in hits
        hists[sample]["total_cluster_energy"] = ROOT.TH1F(f"total_cluster_energy_{sample}", f"Total Cluster Energy per Event - {sample}", 100, 0, 100)
        hists[sample]["total_hit_energy"] = ROOT.TH1F(f"total_hit_energy_{sample}", f"Total Hit Energy per Event - {sample}", 100, 0, 100)
        hists[sample]["cluster_efficiency"] = ROOT.TH1F(f"cluster_efficiency_{sample}", f"Cluster Energy / Hit Energy - {sample}", 100, 0, 2)

        # 2D correlations
        hists2d[sample]["cluster_xy"] = ROOT.TH2F(f"cluster_xy_{sample}", f"Cluster Positions (X vs Y) - {sample}", 100, -2000, 2000, 100, -2000, 2000)
        hists2d[sample]["cluster_rz"] = ROOT.TH2F(f"cluster_rz_{sample}", f"Cluster Positions (R vs Z) - {sample}", 100, -3000, 3000, 100, 0, 2000)
        hists2d[sample]["energy_vs_radius"] = ROOT.TH2F(f"energy_vs_radius_{sample}", f"Cluster Energy vs Radius - {sample}", 100, 0, 2000, 100, 0, 100)
        hists2d[sample]["shape_0_vs_energy"] = ROOT.TH2F(f"shape_0_vs_energy_{sample}", f"Shape Parameter 0 vs Energy - {sample}", 100, 0, 100, 100, 0, 1000)

    return hists, hists2d

def analyze_cluster_properties(cluster, hists, hists2d, sample_name, cluster_num, event_num, debug_mode=False):
    try:
        # Basic cluster energy
        cluster_energy = cluster.getEnergy()
        hists[sample_name]["cluster_energy"].Fill(cluster_energy)
        
        if debug_mode and event_num < 3:
            print(f"    Cluster {cluster_num}: Energy = {cluster_energy:.3f} GeV")
        
    except Exception as e:
        if debug_mode:
            print(f"    Cluster {cluster_num}: Error getting energy: {e}")
        return
    
    try:
        # Cluster position
        cluster_pos = cluster.getPosition()
        x, y, z = cluster_pos[0], cluster_pos[1], cluster_pos[2]
        radius = math.sqrt(x*x + y*y)
        
        hists[sample_name]["cluster_x"].Fill(x)
        hists[sample_name]["cluster_y"].Fill(y)
        hists[sample_name]["cluster_z"].Fill(z)
        hists[sample_name]["cluster_radius"].Fill(radius)
        
        # Fill 2D position plots
        hists2d[sample_name]["cluster_xy"].Fill(x, y)
        hists2d[sample_name]["cluster_rz"].Fill(z, radius)
        hists2d[sample_name]["energy_vs_radius"].Fill(radius, cluster_energy)
        
        if debug_mode and event_num < 3:
            print(f"      Position: ({x:.1f}, {y:.1f}, {z:.1f}) mm, R = {radius:.1f} mm")
            
    except Exception as e:
        if debug_mode and event_num < 3:
            print(f"    Cluster {cluster_num}: Error getting position: {e}")
    
    try:
        # Cluster shape parameters
        shape_params = cluster.getShape()
        if len(shape_params) > 0:
            for i, param in enumerate(shape_params[:6]):  # First 6 shape parameters
                hist_name = f"cluster_shape_{i}"
                if hist_name in hists[sample_name]:
                    hists[sample_name][hist_name].Fill(param)
                    
            # Fill 2D correlation for first shape parameter
            if len(shape_params) > 0:
                hists2d[sample_name]["shape_0_vs_energy"].Fill(cluster_energy, shape_params[0])
            
            if debug_mode and event_num < 3:
                print(f"      Shape params: {[f'{p:.2f}' for p in shape_params[:3]]}")
                
    except Exception as e:
        if debug_mode and event_num < 3:
            print(f"    Cluster {cluster_num}: Error getting shape: {e}")
    
    try:
        # Subdetector energies
        subdet_energies = cluster.getSubdetectorEnergies()
        if len(subdet_energies) > 0:
            for i, energy in enumerate(subdet_energies[:3]):  # First 3 subdetectors
                hist_name = f"subdet_energy_{i}"
                if hist_name in hists[sample_name]:
                    hists[sample_name][hist_name].Fill(energy)
            
            if debug_mode and event_num < 3:
                print(f"      Subdet energies: {[f'{e:.3f}' for e in subdet_energies[:3]]}")
                
    except Exception as e:
        if debug_mode and event_num < 3:
            print(f"    Cluster {cluster_num}: Error getting subdetector energies: {e}")
    
    try:
        # Cluster type (if available)
        cluster_type = cluster.getType()
        hists[sample_name]["cluster_type"].Fill(cluster_type)
        
        if debug_mode and event_num < 3:
            print(f"      Type: {cluster_type}")
            
    except Exception as e:
        if debug_mode and event_num < 3:
            print(f"    Cluster {cluster_num}: Error getting type: {e}")

def calculate_total_hit_energy(event, ecal_collections):
    total_energy = 0
    for coll_name in ecal_collections:
        try:
            collection = event.getCollection(coll_name)
            for i in range(collection.getNumberOfElements()):
                hit = collection.getElementAt(i)
                try:
                    total_energy += hit.getEnergy()
                except:
                    continue
        except:
            continue
    return total_energy

def create_cluster_plots(hists, hists2d):
    
    # Plot 1D histograms - let plotHistograms handle everything
    hist_names = ["cluster_energy", "clusters_per_event", "cluster_x", "cluster_y", "cluster_z", 
                  "cluster_radius", "cluster_shape_0", "cluster_shape_1", "cluster_shape_2",
                  "cluster_shape_3", "cluster_shape_4", "cluster_shape_5",
                  "subdet_energy_0", "subdet_energy_1", "subdet_energy_2", "cluster_type",
                  "total_cluster_energy", "total_hit_energy", "cluster_efficiency"]

    xlabel_map = {
        "cluster_energy": "Cluster Energy [GeV]",
        "clusters_per_event": "Number of Clusters per Event",
        "cluster_x": "Cluster X Position [mm]",
        "cluster_y": "Cluster Y Position [mm]",
        "cluster_z": "Cluster Z Position [mm]",
        "cluster_radius": "Cluster Radius [mm]",
        "cluster_shape_0": "Cluster Shape Parameter 0",
        "cluster_shape_1": "Cluster Shape Parameter 1", 
        "cluster_shape_2": "Cluster Shape Parameter 2",
        "cluster_shape_3": "Cluster Shape Parameter 3",
        "cluster_shape_4": "Cluster Shape Parameter 4",
        "cluster_shape_5": "Cluster Shape Parameter 5",
        "subdet_energy_0": "Subdetector 0 Energy [GeV]",
        "subdet_energy_1": "Subdetector 1 Energy [GeV]",
        "subdet_energy_2": "Subdetector 2 Energy [GeV]",
        "cluster_type": "Cluster Type",
        "total_cluster_energy": "Total Cluster Energy per Event [GeV]",
        "total_hit_energy": "Total Hit Energy per Event [GeV]",
        "cluster_efficiency": "Cluster Energy / Hit Energy"
    }

    for hist_name in hist_names:
        hists_to_plot = {s: hists[s][hist_name] for s in hists if hist_name in hists[s]}
        if hists_to_plot:
            plotHistograms(hists_to_plot, f"plots/cluster_{hist_name}.png",
                          xlabel_map.get(hist_name, hist_name), "Entries")

    # Simple 2D histograms - one canvas per plot
    for sample_name in hists2d:
        for hist_name in hists2d[sample_name]:
            c = ROOT.TCanvas("c", "c", 800, 600)
            hists2d[sample_name][hist_name].Draw("colz")
            c.SaveAs(f"plots/cluster_2d_{sample_name}_{hist_name}.png")
            c.Clear()

def main():
    # File discovery
    samples = glob.glob("/data/fmeloni/DataMuC_MuColl10_v0A/v0/reco/electronGun*")
    files = {}
    
    slices = ["0_50", "50_250", "250_1000", "1000_5000"]
    for s in slices:
        files[f"electronGun_pT_{s}"] = []

    for s in samples:
        sname = s.split("/")[-1]
        if sname not in files:
            continue
        files[sname] = glob.glob(f"{s}/*.slcio")

    # Print file counts
    for sample_name, file_list in files.items():
        print(f"Found {len(file_list)} files for {sample_name}")

    hists, hists2d = setup_cluster_histograms()

    # Create reader
    reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()
    reader.setReadCollectionNames(collections_to_read)

    # Loop over samples
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
                if i % 1000 == 0:
                    print(f"\tProcessing event: {i}")

                try:
                    clusters = event.getCollection("PandoraClusters")
                    num_clusters = clusters.getNumberOfElements()
                    
                    # Fill number of clusters per event
                    hists[s]["clusters_per_event"].Fill(num_clusters)
                    
                    # If no clusters, skip but continue (as you requested)
                    if num_clusters == 0:
                        if i < 3:
                            print(f"  Event {i}: No clusters found, continuing...")
                        i += 1
                        continue
                    
                    # Debug info
                    if i < 3:
                        print(f"  Event {i}: Found {num_clusters} clusters")
                    
                    # Calculate total energies for comparison
                    total_cluster_energy = 0
                    total_hit_energy = calculate_total_hit_energy(event, ecal_collections)
                    
                    # Analyze each cluster
                    for j in range(num_clusters):
                        cluster = clusters.getElementAt(j)
                        analyze_cluster_properties(cluster, hists, hists2d, s, j, i, debug_mode=(i < 3))
                        try:
                            total_cluster_energy += cluster.getEnergy()
                        except:
                            pass
                    
                    # Fill event-level histograms
                    hists[s]["total_cluster_energy"].Fill(total_cluster_energy)
                    hists[s]["total_hit_energy"].Fill(total_hit_energy)
                    if total_hit_energy > 0:
                        efficiency = total_cluster_energy / total_hit_energy
                        hists[s]["cluster_efficiency"].Fill(efficiency)
                        
                    if i < 3:
                        print(f"    Total cluster energy: {total_cluster_energy:.3f} GeV")
                        print(f"    Total hit energy: {total_hit_energy:.3f} GeV")
                        print(f"    Efficiency: {efficiency:.3f}" if total_hit_energy > 0 else "    Efficiency: N/A")

                except Exception as e:
                    if i < 3:
                        print(f"  Event {i}: Error accessing clusters: {e}")

                i += 1

            reader.close()

    # Generate plots
    create_cluster_plots(hists, hists2d)

if __name__ == "__main__":
    main()
