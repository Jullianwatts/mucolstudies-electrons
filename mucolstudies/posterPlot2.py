import math
import pyLCIO
from pyLCIO import EVENT
import ROOT
import os

# PDG Code Map for readable summary
pdg_map = {
    11: "Electrons",
    13: "Muons",
    22: "Photons",
    211: "Charged Pions",
    130: "Neutral Kaons",
    2112: "Neutrons",
    2212: "Protons"
}

file_path = "/scratch/jwatts/mucol/data/reco/electronGun_pT_0_50/electronGun_pT_0_50_reco_7.slcio"
PLOT_DIR = "/scratch/jwatts/mucol/mucolstudies/plots2026"
os.makedirs(PLOT_DIR, exist_ok=True)

# ROOT Setup
ROOT.gROOT.SetBatch()
hists = {}

# Dictionary to store counts
counts = {}
event_count = 0
total_pfos = 0 

reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()

try:
    reader.open(file_path)
    print(f"Reading file: {file_path}")

    # Loop through the events in this single file
    for event in reader:
        event_count += 1
        if "PandoraPFOs" in event.getCollectionNames():
            pfos = event.getCollection("PandoraPFOs")
            total_pfos += pfos.getNumberOfElements()
            
            for pfo in pfos:
                pfo_type = abs(pfo.getType())
                
                # Update counts for summary
                counts[pfo_type] = counts.get(pfo_type, 0) + 1
                
                # Plot for all pfos that are NOT id as pdg id = 11
                if pfo_type != 11:
                    name = pdg_map.get(pfo_type, f"Other_{pfo_type}")
                    if name not in hists:
                        hists[name] = ROOT.TH1F(f"h_{name}", "Non-Electron PFOs; #eta; Counts", 50, -2.5, 2.5)
                    
                    mom = pfo.getMomentum()
                    px, py, pz = mom[0], mom[1], mom[2]
                    p = math.sqrt(px**2 + py**2 + pz**2)
                    
                    if p > 0 and p != abs(pz):
                        eta = 0.5 * math.log((p + pz) / (p - pz))
                        hists[name].Fill(eta)

    reader.close()
except Exception as e:
    print(f"Error: {e}")

# Plotting
canvas = ROOT.TCanvas("c", "c", 800, 600)
stack = ROOT.THStack("stack", "Non-Electron PFOs by Type; #eta; Counts")
legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)

# Standard ROOT colors
colors = [ROOT.kRed, ROOT.kBlue, ROOT.kGreen, ROOT.kOrange, ROOT.kMagenta, ROOT.kCyan, ROOT.kYellow]
for i, (name, hist) in enumerate(hists.items()):
    hist.SetLineColor(colors[i % len(colors)])
    hist.SetLineWidth(2)
    stack.Add(hist)
    legend.AddEntry(hist, name, "l")

stack.Draw("nostack")
legend.Draw()
canvas.SaveAs(os.path.join(PLOT_DIR, "non_electron_pfo_eta.svg"))

# FINAL SUMMARY PRINTING
print("-" * 30)
print(f"SUMMARY FOR {event_count} EVENTS")
print(f"Total PFOs found: {total_pfos}")
print("-" * 30)

if counts:
    sorted_types = sorted(counts.items(), key=lambda item: item[1], reverse=True)
    for pfo_type, count in sorted_types:
        name = pdg_map.get(pfo_type, f"Unknown ({pfo_type})")
        print(f"  - {name:<18}: Total = {count:<8}")
