import os
import glob
import pyLCIO
import ROOT
from ROOT import TH2F
import plotHelper

# Set up options
max_events = 10

# Get all electronGun samples
samples = glob.glob("/data/fmeloni/DataMuC_MuColl10_v0A/v0/reco/electronGun_pT_*")

# Set up file dictionary
files = {}
for s in samples:
    sname = s.split("/")[-1]  # Get directory name like "electronGun_pT_0_50"
    files[sname] = glob.glob(f"{s}/*.slcio")

print(f"Found samples: {list(files.keys())}")

# Set up 2D histograms for cluster shapes
hists2d = {}
for s in files:
    hists2d[s] = {}
    hists2d[s]["barrel_cluster_shape"] = ROOT.TH2F(f"{s}_barrel_shape", f"Barrel Cluster Shape {s};X [mm];Y [mm]", 
                                                   100, -2000, 2000, 100, -2000, 2000)
    hists2d[s]["endcap_cluster_shape"] = ROOT.TH2F(f"{s}_endcap_shape", f"Endcap Cluster Shape {s};X [mm];Y [mm]", 
                                                   100, -2000, 2000, 100, -2000, 2000)

# Loop over samples
for s in files:
    print(f"Working on sample: {s}")
    
    # Loop over files in sample
    for f in files[s]:
        if max_events > 0 and files[s].index(f) >= max_events: 
            break
            
        # Open LCIO file
        try:
            reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()
            reader.open([f])
            
            # Loop over events in file
            for event in reader:
                
                # Get ECAL barrel collection
                try:
                    barrel_collection = event.getCollection("ECalBarrelCollection")
                    if barrel_collection.getNumberOfElements() > 0:
                        for i in range(barrel_collection.getNumberOfElements()):
                            hit = barrel_collection.getElementAt(i)
                            pos = hit.getPosition()
                            x, y = pos[0], pos[1]
                            hists2d[s]["barrel_cluster_shape"].Fill(x, y)
                except:
                    pass
                
                # Get ECAL endcap collection  
                try:
                    endcap_collection = event.getCollection("ECalEndcapCollection")
                    if endcap_collection.getNumberOfElements() > 0:
                        for i in range(endcap_collection.getNumberOfElements()):
                            hit = endcap_collection.getElementAt(i)
                            pos = hit.getPosition()
                            x, y = pos[0], pos[1]
                            hists2d[s]["endcap_cluster_shape"].Fill(x, y)
                except:
                    pass
                    
            reader.close()
            
        except Exception as e:
            print(f"Error with file {f}: {e}")
            continue

# Save results
output_file = ROOT.TFile("simple_ecal_clusters.root", "RECREATE")
for s in files:
    hists2d[s]["barrel_cluster_shape"].Write()
    hists2d[s]["endcap_cluster_shape"].Write()
output_file.Close()

# Create plots
for s in files:
    # Barrel plots
    if hists2d[s]["barrel_cluster_shape"].GetEntries() > 0:
        canvas = ROOT.TCanvas(f"barrel_{s}", f"Barrel {s}", 800, 600)
        hists2d[s]["barrel_cluster_shape"].Draw("COLZ")
        canvas.SaveAs(f"barrel_cluster_shape_{s}.png")
        print(f"Saved: barrel_cluster_shape_{s}.png")
    
    # Endcap plots
    if hists2d[s]["endcap_cluster_shape"].GetEntries() > 0:
        canvas = ROOT.TCanvas(f"endcap_{s}", f"Endcap {s}", 800, 600)
        hists2d[s]["endcap_cluster_shape"].Draw("COLZ")
        canvas.SaveAs(f"endcap_cluster_shape_{s}.png")
        print(f"Saved: endcap_cluster_shape_{s}.png")

print("Analysis complete!")
