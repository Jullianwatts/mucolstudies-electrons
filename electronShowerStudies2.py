
import math
import glob
import ROOT
import pyLCIO

exec(open("./plotHelper.py").read())
ROOT.gROOT.SetBatch()
# Set up options
max_events = -1  # Process all files

# Get all electronGun samples
samples = glob.glob("/data/fmeloni/DataMuC_MuColl10_v0A/v0/reco/electronGun_pT_*")

# Set up file dictionary
files = {}
slices = ["0_50", "50_250", "250_1000", "1000_5000"]
for s in slices: 
    files[f"electronGun_pT_{s}"] = []

for s in samples:
    sname = s.split("/")[-1]
    if sname not in files: 
        continue
    files[sname] = glob.glob(f"{s}/*.slcio")

print(f"Found samples: {list(files.keys())}")

# Create combined histograms for ALL pT slices together
combined_barrel = ROOT.TH2F("combined_barrel_all_pT", "Combined Barrel All pT;X [mm];Y [mm]", 100, -2000, 2000, 100, -2000, 2000)
combined_endcap = ROOT.TH2F("combined_endcap_all_pT", "Combined Endcap All pT;X [mm];Y [mm]", 100, -2000, 2000, 100, -2000, 2000)

# Loop over samples
for s in files:
    print(f"\nWorking on sample: {s}")
    
    # Loop over ALL files in sample
    for file_count, file_path in enumerate(files[s]):
        
        if file_count % 100 == 0:
            print(f"  Processing file {file_count}")
            
        # Open LCIO file
        try:
            reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()
            reader.open(file_path)
            
            # Loop over events in file
            for event in reader:
                
                # Get ECAL barrel hits
                try:
                    barrel_collection = event.getCollection("ECalBarrelCollection")
                    for i in range(barrel_collection.getNumberOfElements()):
                        hit = barrel_collection.getElementAt(i)
                        pos = hit.getPosition()
                        x, y = pos[0], pos[1]
                        combined_barrel.Fill(x, y)
                except:
                    pass
                
                # Get ECAL endcap hits  
                try:
                    endcap_collection = event.getCollection("ECalEndcapCollection")
                    for i in range(endcap_collection.getNumberOfElements()):
                        hit = endcap_collection.getElementAt(i)
                        pos = hit.getPosition()
                        x, y = pos[0], pos[1]
                        combined_endcap.Fill(x, y)
                except:
                    pass
                    
            reader.close()
            
        except:
            continue
    
    print(f"  Completed {s}: processed {len(files[s])} files")

# Use plotHelper's plotHistograms function like your working script
print("\nCreating plots using plotHelper functions...")

# Create histogram dictionaries like your working script does
barrel_hists = {"All_pT_Combined": combined_barrel}
endcap_hists = {"All_pT_Combined": combined_endcap}

# Use plotHelper's functions - this should work without TCanvas errors
plotHistograms(barrel_hists, "plots/combined_barrel_all_pT.png", "X [mm]", "Y [mm]")
print("Saved: plots/combined_barrel_all_pT.png")

plotHistograms(endcap_hists, "plots/combined_endcap_all_pT.png", "X [mm]", "Y [mm]")
print("Saved: plots/combined_endcap_all_pT.png")

print("Analysis complete!")
