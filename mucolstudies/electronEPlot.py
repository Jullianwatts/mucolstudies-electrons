from math import exp, gamma, log, sqrt 
import glob
import ROOT
import pyLCIO
from pyLCIO import EVENT, UTIL
from ROOT import TH1F, TH2F, TFile, TLorentzVector, TMath, TCanvas
import numpy as np
import math
import os

ROOT.gROOT.SetBatch()
# Set to 1000 events
max_events = 1000 

# Path to the specific directory
samples = ["/scratch/jwatts/mucol/data/reco/electronGun_pT_0_50"]

files = {}
slices = ["0_50"] 

for s in slices:
    files[f"electronGun_pT_{s}"] = []

for s in samples:
    sname = s.split("/")[-1]
    if sname not in files:
        continue
    # TARGET ONLY FILE 4
    files[sname] = glob.glob(f"{s}/electronGun_pT_0_50_reco_4.slcio")

print(f"Found files: {files}")

# Initialize longitudinal profile dictionary
longitudinal_profile = {slice_name: {} for slice_name in files}

# Define histogram variables
variables = {
    "electron": {
        "E": {"nbins": 50, "xmin": 0, "xmax": 100},
        "pt": {"nbins": 50, "xmin": 0, "xmax": 100},
        "eta": {"nbins": 30, "xmin": -3, "xmax": 3},
        "phi": {"nbins": 30, "xmin": -3.14, "xmax": 3.14},
        "theta": {"nbins": 30, "xmin": 0, "xmax": 3.14},
        "shower_start_layer": {"nbins": 20, "xmin": 0, "xmax": 20},
        "profile_discrepancy": {"nbins": 100, "xmin": 0, "xmax": 2}
    }
}

# Set up histograms
hists = {}
for s in files:
    if not files[s]:
        continue
    hists[s] = {}
    for obj in ["mcp", "electron"]:
        for var in variables["electron"]:
            hname = f"{s}_{obj}_{var}"
            hists[s][f"{obj}_{var}"] = ROOT.TH1F(hname, hname,
                                              variables["electron"][var]["nbins"],
                                              variables["electron"][var]["xmin"],
                                              variables["electron"][var]["xmax"])

# Analysis Functions
def get_profile_discrepancy_pandora_style(energy_by_layer, total_energy, energy=None):
    if total_energy == 0 or len(energy_by_layer) < 3:
        return -1, 0.0
    max_layer = max(energy_by_layer.keys(), key=lambda k: energy_by_layer[k])
    max_energy_fraction = energy_by_layer[max_layer] / total_energy
    total_discrepancy = 0.0
    n_layers_compared = 0
    for layer in range(max(0, max_layer - 5), min(50, max_layer + 10)):
        if layer in energy_by_layer:
            f_obs = energy_by_layer[layer] / total_energy
            distance_from_max = abs(layer - max_layer)
            f_exp = max_energy_fraction * exp(-distance_from_max / 3.0) if distance_from_max > 0 else max_energy_fraction
            total_discrepancy += abs(f_obs - f_exp) / (f_exp + 0.01)
            n_layers_compared += 1
    return max_layer, (total_discrepancy / n_layers_compared * 0.5) if n_layers_compared > 0 else 0.0

def find_shower_start_layer_absolute(energy_by_layer, threshold_gev):
    if not energy_by_layer: return -1
    for layer in sorted(energy_by_layer.keys()):
        if energy_by_layer[layer] >= threshold_gev: return layer
    return -1

# Create reader
reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()
reader.setReadCollectionNames(["MCParticle", "PandoraClusters", "EcalBarrelCollectionRec", "EcalEndcapCollectionRec"])

total_events_processed = 0

for slice_name in files:
    for f in files[slice_name]:
        if total_events_processed >= max_events: break
        reader.open(f)
        print(f"Working on file: {f}")

        for ievt, event in enumerate(reader):
            if total_events_processed >= max_events: break
            
            hit_energies_by_layer = {}
            try:
                mcpCollection = event.getCollection('MCParticle')
                truePart = mcpCollection[0]
                if abs(truePart.getPDG()) != 11: continue 
            except: continue

            trueE = truePart.getEnergy()
            dp3 = truePart.getMomentum()
            tlv_truth = TLorentzVector()
            tlv_truth.SetPxPyPzE(dp3[0], dp3[1], dp3[2], trueE)

            hists[slice_name]["mcp_E"].Fill(trueE)
            hists[slice_name]["mcp_pt"].Fill(tlv_truth.Perp())

            try:
                clusterCollection = event.getCollection('PandoraClusters')
                best_cluster = None
                best_dR = 0.1
                for i in range(clusterCollection.getNumberOfElements()):
                    cluster = clusterCollection.getElementAt(i)
                    cp = cluster.getPosition()
                    ctlv = TLorentzVector(); ctlv.SetPxPyPzE(cp[0], cp[1], cp[2], cluster.getEnergy())
                    dR = ctlv.DeltaR(tlv_truth)
                    if dR < best_dR:
                        best_cluster = cluster
                        best_dR = dR
                
                if best_cluster:
                    hits = best_cluster.getCalorimeterHits()
                    for j in range(hits.size()):
                        hit = hits[j]
                        if hit.getEnergy() <= 0: continue
                        # CellID decoding for layer
                        layer = (int(hit.getCellID0()) >> 7) & 0x3F 
                        hit_energies_by_layer[layer] = hit_energies_by_layer.get(layer, 0) + hit.getEnergy()
            except: pass

            total_energy = sum(hit_energies_by_layer.values())
            if total_energy > 0:
                ss_layer = find_shower_start_layer_absolute(hit_energies_by_layer, 0.5)
                _, disc = get_profile_discrepancy_pandora_style(hit_energies_by_layer, total_energy)
                
                hists[slice_name]["electron_shower_start_layer"].Fill(ss_layer)
                hists[slice_name]["electron_profile_discrepancy"].Fill(disc)
                hists[slice_name]["electron_E"].Fill(total_energy)
                hists[slice_name]["electron_pt"].Fill(tlv_truth.Perp())

                for layer, energy in hit_energies_by_layer.items():
                    longitudinal_profile[slice_name][layer] = longitudinal_profile[slice_name].get(layer, 0) + energy

            total_events_processed += 1
            if total_events_processed % 100 == 0: print(f" Processed {total_events_processed} events...")
        reader.close()

print(f"Total events processed: {total_events_processed}")

if not os.path.exists("plots"): os.makedirs("plots")

# Plotting
for s in hists:
    for h_key, hist in hists[s].items():
        if hist.GetEntries() > 0:
            c = ROOT.TCanvas("c","c",800,600)
            hist.Draw("HIST")
            c.SaveAs(f"plots/{h_key}.png")

for s in longitudinal_profile:
    if not longitudinal_profile[s]: continue
    layers = sorted(longitudinal_profile[s].keys())
    hist = ROOT.TH1F(f"long_prof_{s}", "Longitudinal Profile", 60, 0, 60)
    for l in layers: 
        hist.SetBinContent(l+1, longitudinal_profile[s][l]/total_events_processed)
    c = ROOT.TCanvas("c","c",800,600)
    hist.Draw("HIST")
    c.SaveAs(f"plots/longitudinal_profile_{s}.png")
