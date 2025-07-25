import math
from math import exp, gamma, log, sqrt
import glob
import ROOT
import pyLCIO
from pyLCIO import EVENT, UTIL
from ROOT import TH1F, TH2F, TFile, TLorentzVector, TMath, TCanvas
import numpy as np

ROOT.gROOT.SetBatch()

# Set up some options
max_events = -1

import os

#samples = glob.glob("/data/fmeloni/DataMuC_MuColl10_v0A/reco/electronGun*")
#samples = glob.glob("/data/fmeloni/DataMuC_MAIA_v0/v3/reco/electronGun*")
#samples = glob.glob("/data/fmeloni/DataMuC_MAIA_v0/v3/electronGun*")
#samples = glob.glob("/data/fmeloni/DataMuC_MAIA_v0/v4rotated/electronGun*")
#samples = glob.glob("/data/fmeloni/DataMuC_MAIA_v0/v4/electronGun*")
#samples = glob.glob("/data/fmeloni/DataMuC_MuColl10_v0A/v0/reco/electronGun*")
samples = glob.glob("/data/fmeloni/DataMuC_MAIA_v0/v5/reco/electronGun*")
files = {}
slices = ["0_50", "50_250", "250_1000", "1000_5000"]
for s in slices:
    files[f"electronGun_pT_{s}"] = []

for s in samples:
    sname = s.split("/")[-1]
    if sname not in files:
        continue
    files[sname] = glob.glob(f"{s}/*.slcio")

print("Found files:")

# Initialize longitudinal profile dictionary
longitudinal_profile = {slice_name: {} for slice_name in files}

# Define histogram variables for LCElectronId parameters
variables = {
    "electron": {
        "E": {"nbins": 50, "xmin": 0, "xmax": 1000},
        "pt": {"nbins": 50, "xmin": 0, "xmax": 1000},
        "eta": {"nbins": 30, "xmin": -3, "xmax": 3},
        "phi": {"nbins": 30, "xmin": -3.14, "xmax": 3.14},
        "theta": {"nbins": 30, "xmin": 0, "xmax": 3.14},
        "shower_start_layer": {"nbins": 20, "xmin": 0, "xmax": 20},
        "shower_max_layer": {"nbins": 20, "xmin": 0, "xmax": 20},
        "max_cell_energy": {"nbins": 50, "xmin": 0, "xmax": 50},
        "profile_discrepancy": {"nbins": 50, "xmin": 0, "xmax": 1},
        "E_over_p": {"nbins": 50, "xmin": 0, "xmax": 2},
        "residual_E_over_p": {"nbins": 50, "xmin": 0, "xmax": 1},
        "cluster_cone_energy": {"nbins": 50, "xmin": 0, "xmax": 1000},
        "cluster_rms_width": {"nbins": 50, "xmin": 0, "xmax": 200},
        "n_hits_in_cluster": {"nbins": 50, "xmin": 0, "xmax": 100},
        "n_layers_with_energy": {"nbins": 20, "xmin": 0, "xmax": 20},
        "track_momentum": {"nbins": 50, "xmin": 0, "xmax": 1000},
        "track_cluster_dR": {"nbins": 50, "xmin": 0, "xmax": 0.5}
    }
}

# Set up histograms
hists = {}
for s in files:
    if not files[s]:  # Skip empty slices
        continue
    hists[s] = {}

    # Truth and reconstructed electron histograms
    for obj in ["mcp", "mcp_el", "electron", "electron_match"]:
        for var in variables["electron"]:
            hists[s][f"{obj}_{var}"] = ROOT.TH1F(f"{s}_{obj}_{var}", f"{s}_{obj}_{var}",
                                              variables["electron"][var]["nbins"],
                                              variables["electron"][var]["xmin"],
                                              variables["electron"][var]["xmax"])

    # Special cluster histograms
    hists[s]["cluster_nhits"] = ROOT.TH1F(f"{s}_cluster_nhits", f"{s}_cluster_nhits", 50, 0, 100)
    hists[s]["cluster_r"] = ROOT.TH1F(f"{s}_cluster_r", f"{s}_cluster_r", 50, 0, 3000)
    hists[s]["pfo_shower_start_layer"] = ROOT.TH1F(f"{s}_pfo_shower_start_layer", f"{s}_pfo_shower_start_layer", 20, 0, 20)
    # LCElectronId parameter histograms
    hists[s]["lcelectronid_max_inner_layer"] = ROOT.TH1F(f"{s}_lcelectronid_max_inner_layer", f"{s}_lcelectronid_max_inner_layer", 10, 0, 10)
    hists[s]["lcelectronid_max_energy"] = ROOT.TH1F(f"{s}_lcelectronid_max_energy", f"{s}_lcelectronid_max_energy", 50, 0, 20)
    hists[s]["lcelectronid_max_profile_start"] = ROOT.TH1F(f"{s}_lcelectronid_max_profile_start", f"{s}_lcelectronid_max_profile_start", 20, 0, 20)
    hists[s]["lcelectronid_max_profile_discrepancy"] = ROOT.TH1F(f"{s}_lcelectronid_max_profile_discrepancy", f"{s}_lcelectronid_max_profile_discrepancy", 100, 0, 2)
    hists[s]["lcelectronid_max_residual_e_over_p"] = ROOT.TH1F(f"{s}_lcelectronid_max_residual_e_over_p", f"{s}_lcelectronid_max_residual_e_over_p", 50, 0, 1)

# Set up 2D histograms
hists2d = {}
for s in files:
    if not files[s]:  # Skip empty slices
        continue
    hists2d[s] = {}

    # Track vs truth comparisons - COMMENTED OUT
    # hists2d[s]["track_eta_v_mcp_eta"] = ROOT.TH2F(f"track_eta_v_mcp_eta_{s}", f"track_eta_v_mcp_eta_{s}", 30, -3, 3, 30, -3, 3)
    # hists2d[s]["track_pt_v_mcp_pt"] = ROOT.TH2F(f"track_pt_v_mcp_pt_{s}", f"track_pt_v_mcp_pt_{s}", 30, 0, 3000, 30, 0, 3000)
    # hists2d[s]["track_phi_v_mcp_phi"] = ROOT.TH2F(f"track_phi_v_mcp_phi_{s}", f"track_phi_v_mcp_phi_{s}", 30, -3, 3, 30, -3, 3)

    # Cluster vs truth comparisons
    hists2d[s]["cluster_E_v_mcp_E"] = ROOT.TH2F(f"cluster_E_v_mcp_E_{s}", f"cluster_E_v_mcp_E_{s}", 30, 0, 1000, 30, 0, 1000)
    hists2d[s]["cluster_eta_v_mcp_eta"] = ROOT.TH2F(f"cluster_eta_v_mcp_eta_{s}", f"cluster_eta_v_mcp_eta_{s}", 30, -3, 3, 30, -3, 3)

    # LCElectronId parameter correlations
    hists2d[s]["shower_start_layer_v_profile_discrepancy"] = ROOT.TH2F(f"shower_start_layer_v_profile_discrepancy_{s}","Shower Start Layer vs Profile Discrepancy;Shower Start Layer;Profile Discrepancy",
    20, 0, 20,       # 20 bins from layer 0 to 20
    100, 0, 2        # 100 bins from discrepancy 0 to 2 (adjust range if needed)
)

    hists2d[s]["E_over_p_v_profile_discrepancy"] = ROOT.TH2F(f"E_over_p_v_profile_discrepancy_{s}", f"E_over_p_v_profile_discrepancy_{s}", 50, 0, 2, 50, 0, 1)
    max_y = 10  # default
    if "1000_5000" in s:
        max_y = 100  # or 200 depending on what you're seeing
    elif "250_1000" in s:
        max_y = 50
    elif "50_250" in s:
        max_y = 20
    hists2d[s]["shower_start_layer_v_max_cell_energy"] = ROOT.TH2F(f"shower_start_layer_v_max_cell_energy_{s}", f"shower_start_layer_v_max_cell_energy_{s}", 20, 0, 20, 50, 0, 10)
    hists2d[s]["cluster_rms_width_v_n_hits"] = ROOT.TH2F(f"cluster_rms_width_v_n_hits_{s}", f"cluster_rms_width_v_n_hits_{s}", 50, 0, 200, 50, 0, 100)

# Matching function
def isMatched(tlv1, tlv2, dR_cut=0.1):
    """Check if two TLorentzVectors are matched"""
    if tlv1.DeltaR(tlv2) < dR_cut:
        return True
    return False

def calculate_rms_width(hit_positions, hit_energies, center):
    """Calculate RMS width of cluster"""
    if len(hit_positions) == 0:
        return 0.0

    total_energy = sum(hit_energies)
    if total_energy == 0:
        return 0.0

    weighted_r_squared = 0.0
    for i, pos in enumerate(hit_positions):
        r_squared = (pos[0] - center[0])**2 + (pos[1] - center[1])**2
        weighted_r_squared += hit_energies[i] * r_squared

    return math.sqrt(weighted_r_squared / total_energy)

def find_shower_start_layer(energy_by_layer, threshold=0.01):
    """Find first layer with significant energy deposition"""
    if not energy_by_layer:
        return -1

    total_energy = sum(energy_by_layer.values())
    if total_energy == 0:
        return -1

    sorted_layers = sorted(energy_by_layer.keys())
    for layer in sorted_layers:
        if energy_by_layer[layer] / total_energy > threshold:
            return layer
    return sorted_layers[0] if sorted_layers else -1
###NEED TO CHECK
def generate_expected_em_profile(num_layers=50, energy=10.0, X0_scale=0.6286):
    # Shape parameters (approximate) from EM shower physics
    a = 1.0 + 0.5 * math.log(energy / 0.01)  # energy in GeV
    b = 0.5  # shower decay rate

    expected_profile = {}
    total_expected = 0.0
    for layer in range(num_layers):
        t = layer * X0_per_layer
        if t == 0 and a < 1:
            profile_value = 0.0
        else:
            try:
                profile_value = (b**a * t**(a-1) * exp(-b*t)) / gamma(a)
            except (OverflowError, ZeroDivisionError):
                profile_value = 0.0
        expected_profile[layer] = profile_value
        total_expected += profile_value
        if total_expected > 0:
            for layer in expected_profile:
                expected_profile[layer] != total_expexted
        return expected_profile


def get_profile_discrepancy(energy_by_layer, total_energy, energy=10.0,
                             min_fraction_per_layer=0.6,
                             num_layers=50, X0_per_layer=0.628):
    ## .63 fraction/layer comes from barrel ecal tungsten being 2.2mm/3.5mm(each layer) = .6286
    if total_energy == 0 or len(energy_by_layer) < 3:
        return -1, -1.0

    # 1. Generate expected profile using Gamma distribution
    expected_profile = generate_expected_em_profile(num_layers=num_layers, energy=energy, X0_per_layer=X0_per_layer)

    # 2. Normalize observed energy fractions
    observed_fractions = {}
    for later, layer_energy in energy_by_layer.items():
        if layer_energy > 0 and layer < num_layers:
            observed_fractions[layer] = layer_energy / total_energy
    if len(observed_fractions) < 2:
        return -1, -1.0

    
    max_discrepancy = -1.0
    max_discrepancy_layer = -1

    all_layers = set(observed_fractions.keys()) | set(expected_profile.keys())
    for layer in all_layers:
        if layer >= num_layers:
            continue
        f_obs = observed_fractions.get(layer, 0.0)
        f_exp = expected_profile.get(layer, 0.0)
        discrepancy = abs(f_obs-f_exp)
        if discrepancy > max_discrepancy:
            max_discrepancy = discrepancy
            max_discrepancy_later, max_discrepancy
    return max_discrepancy_layer, max_discrepancy
def plotHistograms(hist_dict, output_path, x_title, y_title):
    """Plot multiple histograms on same canvas"""
    if not hist_dict:
        return

    c = ROOT.TCanvas("can", "can", 800, 600)
    c.SetLogy()

    colors = [ROOT.kBlue, ROOT.kRed, ROOT.kGreen, ROOT.kMagenta, ROOT.kCyan]
    legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)

    first = True
    for i, (name, hist) in enumerate(hist_dict.items()):
        if hist.GetEntries() == 0:
            continue

        hist.SetLineColor(colors[i % len(colors)])
        hist.SetLineWidth(2)
        hist.GetXaxis().SetTitle(x_title)
        hist.GetYaxis().SetTitle(y_title)

        if first:
            hist.Draw("HIST")
            first = False
        else:
            hist.Draw("HIST SAME")

        legend.AddEntry(hist, name.replace("electronGun_pT_", ""), "l")

    legend.Draw()
    c.SaveAs(output_path)
    c.Close()

# Create a reader object
reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()
reader.setReadCollectionNames(["MCParticle", "SiTracks", "AllTracks", "PandoraPFOs", "SeedTracks", "PandoraClusters",
                              "EcalBarrelCollectionRec", "EcalEndcapCollectionRec"])  #changed from EcalXCollectionSel to Rec

total_events_processed = 0

# Loop over the different samples
for slice_name in files:
    if not files[slice_name]:  # Skip empty slices
        continue

    print(f"Working on sample {slice_name}")

    file_count = 0
    events_processed_in_slice = 0  # Initialize counter for this slice

    # Loop over the files in a sample
    for f in files[slice_name]:
        if max_events > 0 and file_count >= max_events:
            break

        try:
            reader.open(f)
        except Exception as e:
            print(f"  Error opening file {f}: {e}")
            continue
            
        # Loop over events in this file
        for ievt, event in enumerate(reader):

            # Reset variables for this event
            hit_energies_by_layer = {}
            cluster_hit_positions = []
            cluster_hit_energies = []

            # Get truth particle
            try:
                mcpCollection = event.getCollection('MCParticle')
                trueElectron = mcpCollection[0]
            except:
                continue

            # Check if it's actually an electron
            if abs(trueElectron.getPDG()) != 11:
                continue

            # Get truth energy and momentum
            trueE = trueElectron.getEnergy()
            dp3 = trueElectron.getMomentum()
            tlv_truth = TLorentzVector()
            tlv_truth.SetPxPyPzE(dp3[0], dp3[1], dp3[2], trueE)

            # Apply eta cut - only analyze electrons within |eta| < 2.4
            if abs(tlv_truth.Eta()) > 2.4:
                continue

            # Fill truth histograms
            hists[slice_name]["mcp_E"].Fill(trueE)
            hists[slice_name]["mcp_pt"].Fill(tlv_truth.Perp())
            hists[slice_name]["mcp_eta"].Fill(tlv_truth.Eta())
            hists[slice_name]["mcp_phi"].Fill(tlv_truth.Phi())
            hists[slice_name]["mcp_theta"].Fill(tlv_truth.Theta())

            # Try to get track information - COMMENTED OUT
            track_momentum = -1
            track_cluster_dR = -1

            # ECAL hit analysis
            cluster_energy = 0.0
            max_energy = 0.0  # Initialize max_energy
            ecal_coll = ['EcalBarrelCollectionRec', 'EcalEndcapCollectionRec']   #also changed these from Sel to Rec
            max_cluster_energies = []

            for coll in ecal_coll:
                try:
                    ECALhitCollection = event.getCollection(coll)
                    if ECALhitCollection is None:
                        continue
                        
                    encoding = ECALhitCollection.getParameters().getStringVal(EVENT.LCIO.CellIDEncoding)
                    decoder = UTIL.BitField64(encoding)
                    cluster_hit_energies_local = []
                    
                    for hit in ECALhitCollection:
                        if hit is None:
                            continue
                        hit_energy = hit.getEnergy()
                        if hit_energy > 0:
                            cluster_hit_energies_local.append(hit_energy)
                            
                    if cluster_hit_energies_local:
                        max_energy_local = max(cluster_hit_energies_local)
                        max_cluster_energies.append(max_energy_local)
                        
                except Exception as e:
                    continue  # ECAL collection not available
                    
            # Set the overall max energy
            max_energy = max(max_cluster_energies) if max_cluster_energies else 0.0

            ecal_decoders = {}
            ecal_collections = ['EcalBarrelCollectionRec', 'EcalEndcapCollectionRec']

            for coll_name in ecal_collections:
                try:
                    ECALhitCollection = event.getCollection(coll_name)
                    encoding = ECALhitCollection.getParameters().getStringVal(EVENT.LCIO.CellIDEncoding)
                    decoder = UTIL.BitField64(encoding)
                    ecal_decoders[coll_name] = decoder
                except:
                    continue

            try:
                collection_names = event.getCollectionNames()
                if "PandoraPFOs" in collection_names:
                    pfoCollection = event.getCollection('PandoraPFOs')
                    if pfoCollection is not None:
                        # Process ALL PFOs efficiently
                        for i in range(pfoCollection.getNumberOfElements()):
                            pfo = pfoCollection.getElementAt(i)
                            if pfo is not None and abs(pfo.getCharge()) > 0.5:
                                pfo_clusters = pfo.getClusters()
                                if pfo_clusters is not None:
                                    # Process ALL clusters
                                    for j in range(pfo_clusters.size()):
                                        cluster = pfo_clusters[j]
                                        if cluster is not None:
                                            pfo_hit_energies_by_layer = {}
                                            hits = cluster.getCalorimeterHits()
                                            if hits is not None:
                                                # Process ALL hits efficiently
                                                for k in range(hits.size()):
                                                    hit = hits[k]
                                                    if hit is not None:
                                                        hit_energy = hit.getEnergy()
                                                        if hit_energy > 0:
                                                            cellID = int(hit.getCellID0())
                                                            # Use pre-loaded decoders
                                                            for coll_name, decoder in ecal_decoders.items():
                                                                try:
                                                                    decoder.setValue(cellID)
                                                                    layer = decoder["layer"].value()
                                                                    if layer not in pfo_hit_energies_by_layer:
                                                                        pfo_hit_energies_by_layer[layer] = 0.0
                                                                    pfo_hit_energies_by_layer[layer] += hit_energy
                                                                    break
                                                                except:
                                                                    pass

                                            # Calculate shower start for this cluster
                                            if pfo_hit_energies_by_layer:
                                                pfo_shower_start_layer = find_shower_start_layer(pfo_hit_energies_by_layer, threshold=0.1)
                                                if pfo_shower_start_layer >= 0:
                                                    hists[slice_name]["pfo_shower_start_layer"].Fill(pfo_shower_start_layer)
            except:
                pass  # Ignore all PFO analysis errors and continue
                
            # Store longitudinal profile information
            for layer, energy in hit_energies_by_layer.items():
                if layer not in longitudinal_profile[slice_name]:
                    longitudinal_profile[slice_name][layer] = 0.0
                longitudinal_profile[slice_name][layer] += energy

            # Calculate cluster properties
            n_hits_in_cluster = len(cluster_hit_energies)
            n_layers_with_energy = len(hit_energies_by_layer)

            # Calculate shower start layer
            shower_start_layer = find_shower_start_layer(hit_energies_by_layer, threshold=0.01)

            # Find shower maximum layer
            shower_max_layer = -1
            if hit_energies_by_layer:
                shower_max_layer = max(hit_energies_by_layer.keys(), key=lambda k: hit_energies_by_layer[k])

            # Calculate profile discrepancy
            total_energy = sum(hit_energies_by_layer.values())
            profile_discrepancy_layer, profile_discrepancy = get_profile_discrepancy(hit_energies_by_layer, total_energy, energy=cluster_energy, X0_per_layer=0.6286, num_layers = 50
                    ) 
            # Calculate cluster shape properties
            cluster_rms_width = -1
            if cluster_hit_positions and cluster_hit_energies:
                total_energy = sum(cluster_hit_energies)
                if total_energy > 0:
                    center_x = sum(pos[0] * energy for pos, energy in zip(cluster_hit_positions, cluster_hit_energies)) / total_energy
                    center_y = sum(pos[1] * energy for pos, energy in zip(cluster_hit_positions, cluster_hit_energies)) / total_energy
                    center_z = sum(pos[2] * energy for pos, energy in zip(cluster_hit_positions, cluster_hit_energies)) / total_energy
                    cluster_center = [center_x, center_y, center_z]

                    cluster_rms_width = calculate_rms_width(cluster_hit_positions, cluster_hit_energies, cluster_center)

                    # Fill cluster position histogram
                    cluster_r = math.sqrt(center_x**2 + center_y**2)
                    hists[slice_name]["cluster_r"].Fill(cluster_r)

            # Fill cluster histograms
            hists[slice_name]["cluster_nhits"].Fill(n_hits_in_cluster)
            hists[slice_name]["electron_n_hits_in_cluster"].Fill(n_hits_in_cluster)
            hists[slice_name]["electron_n_layers_with_energy"].Fill(n_layers_with_energy)
            hists[slice_name]["electron_cluster_rms_width"].Fill(cluster_rms_width if cluster_rms_width > 0 else 0)

            # Fill LCElectronId parameter histograms
            hists[slice_name]["electron_shower_start_layer"].Fill(shower_start_layer if shower_start_layer > 0 else 0)
            hists[slice_name]["electron_shower_max_layer"].Fill(shower_max_layer if shower_max_layer > 0 else 0)
            hists[slice_name]["electron_max_cell_energy"].Fill(max_energy)
            hists[slice_name]["electron_profile_discrepancy"].Fill(profile_discrepancy)
            hists[slice_name]["electron_cluster_cone_energy"].Fill(cluster_energy)

            # Fill LCElectronId specific histograms
            hists[slice_name]["lcelectronid_max_profile_start"].Fill(shower_start_layer if shower_start_layer > 0 else 0)
            hists[slice_name]["lcelectronid_max_energy"].Fill(max_energy)
            hists[slice_name]["lcelectronid_max_profile_discrepancy"].Fill(profile_discrepancy)

            # Fill reconstructed quantities if we have a cluster
            if cluster_energy > 0:
                hists[slice_name]["electron_E"].Fill(cluster_energy)
                hists[slice_name]["electron_cluster_cone_energy"].Fill(cluster_energy)

                # Calculate E/p ratio - COMMENTED OUT since we're not using tracks
                # E_over_p = -1
                # residual_E_over_p = -1
                # if track_momentum > 0:
                #     E_over_p = cluster_energy / track_momentum
                #     residual_E_over_p = abs(E_over_p - 1.0)
                #     
                #     hists[slice_name]["electron_E_over_p"].Fill(E_over_p)
                #     hists[slice_name]["electron_residual_E_over_p"].Fill(residual_E_over_p)
                #     hists[slice_name]["lcelectronid_max_residual_e_over_p"].Fill(residual_E_over_p)
                #     
                #     # Fill 2D correlations
                #     hists2d[slice_name]["E_over_p_v_profile_discrepancy"].Fill(E_over_p, profile_discrepancy)

                # Energy-weighted average angles
                if cluster_hit_positions and cluster_hit_energies:
                    total_energy = sum(cluster_hit_energies)
                    avg_theta = sum(np.arctan2(np.sqrt(pos[0]**2 + pos[1]**2), pos[2]) * energy
                                  for pos, energy in zip(cluster_hit_positions, cluster_hit_energies)) / total_energy
                    avg_phi = sum(np.arctan2(pos[1], pos[0]) * energy
                                for pos, energy in zip(cluster_hit_positions, cluster_hit_energies)) / total_energy

                    hists[slice_name]["electron_theta"].Fill(avg_theta)
                    hists[slice_name]["electron_phi"].Fill(avg_phi)
                    hists[slice_name]["electron_eta"].Fill(-np.log(np.tan(avg_theta / 2.0)) if avg_theta > 0 else 0)
                    hists[slice_name]["electron_pt"].Fill(cluster_energy * np.sin(avg_theta))

                    # Fill 2D cluster vs truth
                    hists2d[slice_name]["cluster_E_v_mcp_E"].Fill(cluster_energy, trueE)
                    hists2d[slice_name]["cluster_eta_v_mcp_eta"].Fill(-np.log(np.tan(avg_theta / 2.0)) if avg_theta > 0 else 0, tlv_truth.Eta())

                # Fill 2D parameter correlations
                hists2d[slice_name]["shower_start_layer_v_max_cell_energy"].Fill(shower_start_layer if shower_start_layer > 0 else 0, max_energy)
                hists2d[slice_name]["cluster_rms_width_v_n_hits"].Fill(cluster_rms_width if cluster_rms_width > 0 else 0, n_hits_in_cluster)
                hists2d[slice_name]["shower_start_layer_v_profile_discrepancy"].Fill(shower_start_layer if shower_start_layer > 0 else 0, profile_discrepancy)

                # Check if this is a matched electron
                if isMatched(tlv_truth, TLorentzVector()) or cluster_energy > 0:  # Simple matching criteria
                    hists[slice_name]["electron_match_E"].Fill(cluster_energy)
                    hists[slice_name]["electron_match_pt"].Fill(cluster_energy * np.sin(avg_theta) if cluster_hit_positions else 0)

            events_processed_in_slice += 1
            total_events_processed += 1

            # Print progress every 100 events
            if events_processed_in_slice % 100 == 0:
                print(f"  Processed {events_processed_in_slice} events in {slice_name}")
                
        reader.close()
        file_count += 1

print(f"Total events processed: {total_events_processed}")

# Create longitudinal profile plots
print("Creating longitudinal profile plots...")
for s in longitudinal_profile:
    if not longitudinal_profile[s]:
        continue

    # Sort layers
    layers = sorted(longitudinal_profile[s].keys())
    energy_vals = [longitudinal_profile[s][l] for l in layers]

    # Normalize if you want average per event
    if total_events_processed > 0:
        energy_vals = [e / total_events_processed for e in energy_vals]

    hist = ROOT.TH1F(f"long_profile_{s}", f"Longitudinal Profile {s}",
                    len(layers), 0, len(layers))
    for i, e in enumerate(energy_vals):
        hist.SetBinContent(i+1, e)

    c = ROOT.TCanvas("can", "can", 800, 600)
    hist.GetXaxis().SetTitle("ECAL Layer")
    hist.GetYaxis().SetTitle("Avg Energy Deposited [GeV]")
    hist.Draw("HIST")
    c.SaveAs(f"plots/longitudinal_profile_{s}.png")
    c.Close()

print("Creating plots...")

# Plot cluster histograms
cluster_hists_nhits = {}
cluster_hists_r = {}

for s in hists:
    if "cluster_nhits" in hists[s]:
        cluster_hists_nhits[s] = hists[s]["cluster_nhits"]
    if "cluster_r" in hists[s]:
        cluster_hists_r[s] = hists[s]["cluster_r"]

if cluster_hists_nhits:
    plotHistograms(cluster_hists_nhits, "plots/cluster_nhits.png", "Number of Hits per Cluster", "Entries")
if cluster_hists_r:
    plotHistograms(cluster_hists_r, "plots/cluster_r.png", "Cluster Radial Position [mm]", "Entries")

# Plot LCElectronId parameter histograms
lcelectronid_hists = {}
for param in ["shower_start_layer", "max_cell_energy", "profile_discrepancy", "cluster_cone_energy"]:
    lcelectronid_hists[param] = {}
    for s in hists:
        if f"electron_{param}" in hists[s]:
            lcelectronid_hists[param][s] = hists[s][f"electron_{param}"]

    if lcelectronid_hists[param]:
        plotHistograms(lcelectronid_hists[param], f"plots/lcelectronid_{param}.png",
                      param.replace("_", " ").title(), "Entries")

pfo_shower_hists = {}
for s in hists:
    if "pfo_shower_start_layer" in hists[s]:
        pfo_shower_hists[s] = hists[s]["pfo_shower_start_layer"]

if pfo_shower_hists:
    plotHistograms(pfo_shower_hists, "plots/pfo_shower_start_layer.png", "PFO Shower Start Layer", "Entries")
    for s, hist in pfo_shower_hists.items():
        if hist.GetEntries() > 0:
            c = ROOT.TCanvas("can", "can", 800, 600)
            c.SetLogy()
            hist.SetLineColor(ROOT.kBlue)
            hist.SetLineWidth(2)
            hist.GetXaxis().SetTitle("PFO Shower Start Layer")
            hist.GetYaxis().SetTitle("Entries")
            hist.Draw("HIST")
            c.SaveAs(f"plots/pfo_shower_start_layer_{s}.png")
            c.Close()

# Plot 2D histograms
for s in hists2d:
    for h in hists2d[s]:
        if hists2d[s][h].GetEntries() == 0:
            continue

        c = ROOT.TCanvas("can", "can", 800, 600)
        hists2d[s][h].Draw("colz")
        hists2d[s][h].GetXaxis().SetTitle(h.split("_v_")[0].replace("_", " ").title())
        hists2d[s][h].GetYaxis().SetTitle(h.split("_v_")[1].replace("_", " ").title())
        c.SaveAs(f"plots/{hists2d[s][h].GetName()}.png")
        c.Close()

print("Plots saved in plots/ directory")
print("Key LCElectronId parameter plots:")
print("- plots/lcelectronid_shower_start_layer.png (for m_maxProfileStart)")
print("- plots/lcelectronid_max_cell_energy.png (for m_maxEnergy)")
print("- plots/lcelectronid_profile_discrepancy.png (for m_maxProfileDiscrepancy)")
print("- plots/lcelectronid_cluster_cone_energy.png (cluster energy in cone)")
print("Longitudinal profile plots:")
for s in longitudinal_profile:
    if longitudinal_profile[s]:
        print(f"- plots/longitudinal_profile_{s}.png")
