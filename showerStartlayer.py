#!/usr/bin/env python3

import ROOT
import glob
import os
from pyLCIO import IOIMPL, EVENT, UTIL
from plotHelper import plotHistograms

# Set ROOT to batch mode
ROOT.gROOT.SetBatch()

# Configuration
max_events = 10
plots_dir = "plots"

# Detector geometry parameters (from MAIA detector geometry table)
ECAL_LAYERS = 50
LAYER_THICKNESS = 5.35  # mm
ECAL_START_Z = 2307.0  # mm (230.7 cm from geometry table)
ECAL_END_Z = 2575.0    # mm (257.5 cm from geometry table)

# Create plots directory if it doesn't exist
if not os.path.exists(plots_dir):
    os.makedirs(plots_dir)

def get_shower_start_layer(cluster):
    """
    Calculate the shower start layer for a given cluster.
    Returns the layer number (0-49) where the shower starts depositing energy.
    Uses MAIA detector geometry: ECAL starts at Z=230.7cm, 50 layers of 5.35mm each.
    """
    calorimeter_hits = cluster.getCalorimeterHits()
    
    if len(calorimeter_hits) == 0:
        return -1  # No hits in cluster
    
    # Find the minimum z position of all hits in the cluster
    min_z = float('inf')
    
    for hit in calorimeter_hits:
        pos = hit.getPosition()
        z_pos = abs(pos[2])  # Take absolute value in case of negative z
        if z_pos < min_z:
            min_z = z_pos
    
    # Check if hit is actually in the ECAL
    if min_z < ECAL_START_Z or min_z > ECAL_END_Z:
        return -1  # Hit is outside ECAL range
    
    # Convert z position to layer number
    # Layer 0 starts at ECAL_START_Z
    layer = int((min_z - ECAL_START_Z) / LAYER_THICKNESS)
    
    # Clamp to valid layer range
    if layer < 0:
        layer = 0
    elif layer >= ECAL_LAYERS:
        layer = ECAL_LAYERS - 1
        
    return layer

def process_files():
    """
    Process all SLCIO files and create shower start layer histograms.
    """
    
    # Find sample files
    samples = glob.glob("/data/fmeloni/DataMuC_MAIA_v0/v5/reco/electronGun*")
    files = {}
    slices = ["0_50", "50_250", "250_1000", "1000_5000"]
    
    # Initialize file lists for each slice
    for s in slices:
        files[f"electronGun_pT_{s}"] = []
    
    # Collect SLCIO files for each sample
    for s in samples:
        sname = s.split("/")[-1]
        if sname not in files:
            continue
        files[sname] = glob.glob(f"{s}/*.slcio")
    
    # Create histograms for each pT slice
    histograms = {}
    for slice_name in files.keys():
        if len(files[slice_name]) == 0:
            continue
            
        hist_name = f"h_shower_start_{slice_name}"
        histograms[slice_name] = ROOT.TH1F(
            hist_name, 
            f"Shower Start Layer - {slice_name}",
            21,  # 21 bins to cover layers 0-20 (most showers start early)
            -0.5, 
            20.5
        )
    
    # Process each pT slice
    for slice_name, file_list in files.items():
        if len(file_list) == 0:
            print(f"No files found for {slice_name}")
            continue
            
        print(f"Processing {slice_name} with {len(file_list)} files...")
        
        event_count = 0
        
        # Process each file in the slice
        for file_path in file_list:
            if max_events > 0 and event_count >= max_events:
                break
                
            try:
                # Open LCIO file
                reader = IOIMPL.LCFactory.getInstance().createLCReader()
                reader.open(file_path)
                
                # Process events in the file
                while True:
                    if max_events > 0 and event_count >= max_events:
                        break
                        
                    try:
                        event = reader.readNextEvent()
                        if event is None:
                            break
                            
                        event_count += 1
                        if event_count % 1000 == 0:
                            print(f"  Processed {event_count} events...")
                        
                        # Get PandoraClusters collection
                        try:
                            clusters = event.getCollection("PandoraClusters")
                        except:
                            print(f"  No PandoraClusters collection found in event {event_count}")
                            continue
                        
                        # Process each cluster in the event
                        for cluster in clusters:
                            shower_start_layer = get_shower_start_layer(cluster)
                            
                            if shower_start_layer >= 0:  # Valid layer
                                histograms[slice_name].Fill(shower_start_layer)
                                
                    except Exception as e:
                        print(f"  Error processing event {event_count}: {e}")
                        break
                
                reader.close()
                
            except Exception as e:
                print(f"  Error opening file {file_path}: {e}")
                continue
        
        print(f"  Finished processing {slice_name}: {event_count} events, {histograms[slice_name].GetEntries()} cluster entries")
    
    return histograms

def create_plots(histograms):
    """
    Create and save the shower start layer plots using plotHelper functions.
    """
    
    if len(histograms) == 0:
        print("No histograms to plot!")
        return
    
    # Filter out empty histograms
    non_empty_hists = {k: v for k, v in histograms.items() if v.GetEntries() > 0}
    
    if len(non_empty_hists) == 0:
        print("All histograms are empty!")
        return
    
    print(f"Creating plot with {len(non_empty_hists)} histograms...")
    
    # Create the main shower start layer plot with custom styling
    save_path = os.path.join(plots_dir, "shower_start_layer.png")
    
    # Create custom plot instead of using plotHistograms to get the exact styling
    from plotHelper import colors, getMaximum, colorHists
    
    can = ROOT.TCanvas("can", "can", 800, 600)
    can.SetLeftMargin(0.12)
    can.SetBottomMargin(0.12)
    can.SetTopMargin(0.12)
    can.SetRightMargin(0.05)
    # No grid lines
    
    h_keys = list(non_empty_hists.keys())
    h_values = list(non_empty_hists.values())
    
    # Get maxes/mins of all hists
    maxy = 1.5*getMaximum(h_values)
    miny = 1e-1  # For log scale
    can.SetLogy(1)
    maxy *= 1e4
    
    # Draw histograms
    colorHists(h_values)
    non_empty_hists[h_keys[0]].SetTitle("")
    non_empty_hists[h_keys[0]].GetXaxis().SetTitle("Shower Start Layer")
    non_empty_hists[h_keys[0]].GetYaxis().SetTitle("Entries")
    non_empty_hists[h_keys[0]].GetXaxis().SetTitleSize(0.045)
    non_empty_hists[h_keys[0]].GetYaxis().SetTitleSize(0.045)
    non_empty_hists[h_keys[0]].GetXaxis().SetLabelSize(0.04)
    non_empty_hists[h_keys[0]].GetYaxis().SetLabelSize(0.04)
    h_values[0].SetMinimum(miny)
    h_values[0].SetMaximum(maxy)
    non_empty_hists[h_keys[0]].Draw("hist")
    
    for k in h_keys[1:]:
        non_empty_hists[k].Draw("hist same")
    
    # Create legend
    leg = ROOT.TLegend(0.65, 0.15, 0.93, 0.45)
    leg.SetBorderSize(1)
    leg.SetFillColor(0)
    leg.SetFillStyle(1001)
    leg.SetTextSize(0.035)
    leg.SetHeader("pT slices", "C")
    
    for k in h_keys:
        clean_label = k.replace("electronGun_pT_", "").replace("_", "-") + " GeV"
        leg.AddEntry(non_empty_hists[k], clean_label, "l")
    leg.Draw()
    
    # Add MAIA-style text (like in plotEfficiencies)
    text = ROOT.TLatex()
    text.SetNDC()
    text.SetTextAlign(11)  # Left aligned
    text.SetTextSize(0.04)
    text.SetTextFont(62)  # Bold font
    
    concept_text = ROOT.TLatex()
    concept_text.SetNDC()
    concept_text.SetTextAlign(31)  # Right aligned
    concept_text.SetTextSize(0.035)
    concept_text.SetTextFont(42)
    concept_text.DrawLatex(0.94, 0.92, "MAIA Detector Concept")
    
    text.DrawLatex(0.15, 0.92, "Muon Collider")
    text.SetTextFont(42)  # Regular font
    text.DrawLatex(0.15, 0.88, "Simulation, no BIB")
    text.DrawLatex(0.15, 0.84, "|#eta| < 2.4")
    
    ROOT.gStyle.SetOptStat(0)
    can.SaveAs(save_path)
    can.Close()
    
    print(f"Plot saved to: {save_path}")
    
    # Print some statistics
    print("\nHistogram Statistics:")
    for name, hist in non_empty_hists.items():
        entries = hist.GetEntries()
        mean = hist.GetMean()
        rms = hist.GetRMS()
        print(f"  {name}: {entries:.0f} entries, mean={mean:.2f}, RMS={rms:.2f}")

def main():
    """
    Main function to run the analysis.
    """
    print("Starting Shower Start Layer Analysis...")
    print(f"Looking for files in: /data/fmeloni/DataMuC_MAIA_v0/v5/reco/electronGun*")
    print(f"Plots will be saved to: {plots_dir}/")
    print(f"Detector geometry: {ECAL_LAYERS} layers, {LAYER_THICKNESS} mm each")
    print(f"ECAL range: {ECAL_START_Z/10:.1f} - {ECAL_END_Z/10:.1f} cm")
    
    # Process all files and create histograms
    histograms = process_files()
    
    # Create and save plots
    create_plots(histograms)
    
    print("Analysis complete!")

if __name__ == "__main__":
    main()
