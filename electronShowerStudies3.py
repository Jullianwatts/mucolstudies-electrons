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
        "profile_discrepancy": {"nbins": 100, "xmin": 0, "xmax": 2},  # Only Pandora method now
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
    
    # LCElectronId parameter histograms - only Pandora method
    hists[s]["lcelectronid_max_inner_layer"] = ROOT.TH1F(f"{s}_lcelectronid_max_inner_layer", f"{s}_lcelectronid_max_inner_layer", 10, 0, 10)
    hists[s]["lcelectronid_max_energy"] = ROOT.TH1F(f"{s}_lcelectronid_max_energy", f"{s}_lcelectronid_max_energy", 50, 0, 20)
    hists[s]["lcelectronid_max_profile_start"] = ROOT.TH1F(f"{s}_lcelectronid_max_profile_start", f"{s}_lcelectronid_max_profile_start", 20, 0, 20)
    hists[s]["lcelectronid_max_profile_discrepancy"] = ROOT.TH1F(f"{s}_lcelectronid_max_profile_discrepancy", f"{s}_lcelectronid_max_profile_discrepancy", 100, 0, 2)
    hists[s]["lcelectronid_max_residual_e_over_p"] = ROOT.TH1F(f"{s}_lcelectronid_max_residual_e_over_p", f"{s}_lcelectronid_max_residual_e_over_p", 50, 0, 1)
    
    # E/p analysis histograms
    hists[s]["electron_E_over_p"] = ROOT.TH1F(f"{s}_electron_E_over_p", f"E/p Distribution {s};E/p;Entries", 100, 0, 3)
    hists[s]["electron_E_minus_p_over_p"] = ROOT.TH1F(f"{s}_electron_E_minus_p_over_p", f"(E-p)/p Distribution {s};(E-p)/p;Entries", 100, -1, 1)
    hists[s]["electron_abs_E_minus_p_over_p"] = ROOT.TH1F(f"{s}_electron_abs_E_minus_p_over_p", f"|E-p|/p Distribution {s};|E-p|/p;Entries", 100, 0, 1)
    hists[s]["track_momentum"] = ROOT.TH1F(f"{s}_track_momentum", f"Track Momentum {s};p [GeV];Entries", 100, 0, 1000)
    hists[s]["matched_cluster_energy"] = ROOT.TH1F(f"{s}_matched_cluster_energy", f"Matched Cluster Energy {s};E [GeV];Entries", 100, 0, 1000)
    hists[s]["track_cluster_dR"] = ROOT.TH1F(f"{s}_track_cluster_dR", f"Track-Cluster dR {s};dR;Entries", 50, 0, 0.5)

# Set up 2D histograms
hists2d = {}
for s in files:
    if not files[s]:  # Skip empty slices
        continue
    hists2d[s] = {}

    # Cluster vs truth comparisons
    hists2d[s]["cluster_E_v_mcp_E"] = ROOT.TH2F(f"cluster_E_v_mcp_E_{s}", f"cluster_E_v_mcp_E_{s}", 30, 0, 1000, 30, 0, 1000)
    hists2d[s]["cluster_eta_v_mcp_eta"] = ROOT.TH2F(f"cluster_eta_v_mcp_eta_{s}", f"cluster_eta_v_mcp_eta_{s}", 30, -3, 3, 30, -3, 3)

    # LCElectronId parameter correlations
    hists2d[s]["shower_start_layer_v_profile_discrepancy"] = ROOT.TH2F(f"shower_start_layer_v_profile_discrepancy_{s}","Shower Start Layer vs Profile Discrepancy;Shower Start Layer;Profile Discrepancy",
    20, 0, 20, 100, 0, 2)

    hists2d[s]["E_over_p_v_profile_discrepancy"] = ROOT.TH2F(f"E_over_p_v_profile_discrepancy_{s}", f"E_over_p_v_profile_discrepancy_{s}", 50, 0, 2, 50, 0, 1)
    
    max_y = 10  # default
    if "1000_5000" in s:
        max_y = 100
    elif "250_1000" in s:
        max_y = 50
    elif "50_250" in s:
        max_y = 20
    hists2d[s]["shower_start_layer_v_max_cell_energy"] = ROOT.TH2F(f"shower_start_layer_v_max_cell_energy_{s}", f"shower_start_layer_v_max_cell_energy_{s}", 20, 0, 20, 50, 0, 10)
    hists2d[s]["cluster_rms_width_v_n_hits"] = ROOT.TH2F(f"cluster_rms_width_v_n_hits_{s}", f"cluster_rms_width_v_n_hits_{s}", 50, 0, 200, 50, 0, 100)
    
    # E/p analysis 2D histograms
    hists2d[s]["E_vs_p"] = ROOT.TH2F(f"E_vs_p_{s}", f"Cluster Energy vs Track Momentum {s};Track p [GeV];Cluster E [GeV]", 50, 0, 1000, 50, 0, 1000)
    hists2d[s]["E_over_p_vs_energy"] = ROOT.TH2F(f"E_over_p_vs_energy_{s}", f"E/p vs Energy {s};True Energy [GeV];E/p", 50, 0, 1000, 50, 0, 3)
    hists2d[s]["E_over_p_vs_profile_discrepancy"] = ROOT.TH2F(f"E_over_p_vs_profile_discrepancy_{s}", f"E/p vs Profile Discrepancy {s};E/p;Profile Discrepancy", 50, 0, 3, 50, 0, 2)

# PANDORA-STYLE PROFILE DISCREPANCY FUNCTION (ONLY)

def get_profile_discrepancy_pandora_style(energy_by_layer, total_energy, energy=None,
                                        num_layers=50, X0_per_layer=0.6286):
    """
    Profile discrepancy calculation using Pandora's approach
    """
    if total_energy == 0 or len(energy_by_layer) < 3:
        return -1, 0.0
    
    if energy is None:
        energy = total_energy
    
    # Find shower maximum layer
    max_layer = max(energy_by_layer.keys(), key=lambda k: energy_by_layer[k])
    max_energy_fraction = energy_by_layer[max_layer] / total_energy
    
    # Pandora-style shower profile analysis
    # Focus on the region around shower maximum
    start_layer = max(0, max_layer - 5)
    end_layer = min(num_layers, max_layer + 10)
    
    # Calculate expected vs observed in shower development region
    total_discrepancy = 0.0
    n_layers_compared = 0
    max_discrepancy = 0.0
    
    for layer in range(start_layer, end_layer):
        if layer in energy_by_layer:
            f_obs = energy_by_layer[layer] / total_energy
            
            # Simple expected profile based on distance from maximum
            distance_from_max = abs(layer - max_layer)
            if distance_from_max == 0:
                f_exp = max_energy_fraction
            else:
                # Exponential decay from maximum
                decay_constant = 3.0  # Tunable parameter
                f_exp = max_energy_fraction * exp(-distance_from_max / decay_constant)
            
            # Relative discrepancy
            discrepancy = abs(f_obs - f_exp) / (f_exp + 0.01)
            total_discrepancy += discrepancy
            n_layers_compared += 1
            
            if discrepancy > max_discrepancy:
                max_discrepancy = discrepancy
    
    if n_layers_compared > 0:
        avg_discrepancy = total_discrepancy / n_layers_compared
    else:
        avg_discrepancy = 0.0
    
    # Scale the discrepancy to make it comparable to Pandora's 0.6 threshold
    scaled_discrepancy = avg_discrepancy * 0.5  # Scaling factor
    
    return max_layer, scaled_discrepancy

# TRACK-CLUSTER MATCHING AND E/P FUNCTIONS

def match_track_to_cluster(track_tlv, cluster_tlv, dR_cut=0.1):
    """
    Match track to cluster using delta R
    Returns True if matched, False otherwise
    """
    if track_tlv.E() <= 0 or cluster_tlv.E() <= 0:
        return False
    
    dR = track_tlv.DeltaR(cluster_tlv)
    return dR < dR_cut

def calculate_e_over_p_metrics(cluster_energy, track_momentum):
    """
    Calculate E/p related metrics
    Returns: E/p, (E-p)/p, |E-p|/p
    """
    if track_momentum <= 0:
        return -1, -999, -1
        
    e_over_p = cluster_energy / track_momentum
    e_minus_p_over_p = (cluster_energy - track_momentum) / track_momentum
    abs_e_minus_p_over_p = abs(e_minus_p_over_p)
    
    return e_over_p, e_minus_p_over_p, abs_e_minus_p_over_p

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

    return sqrt(weighted_r_squared / total_energy)

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

def plotHistograms(hist_dict, output_path, x_title, y_title):
    """Plot multiple histograms on same canvas"""
    if not hist_dict:
        return

    c = ROOT.TCanvas("can", "can", 800, 600)
    c.SetLogy()

    colors = [ROOT.kBlue, ROOT.kRed, ROOT.kGreen+2, ROOT.kMagenta, ROOT.kCyan+1]
    legend = ROOT.TLegend(0.65, 0.65, 0.88, 0.88)

    first = True
    for i, (name, hist) in enumerate(hist_dict.items()):
        if hist.GetEntries() == 0:
            continue

        hist.SetLineColor(colors[i % len(colors)])
        hist.SetLineWidth(3)
        hist.GetXaxis().SetTitle(x_title)
        hist.GetYaxis().SetTitle(y_title)

        if first:
            hist.Draw("HIST")
            first = False
        else:
            hist.Draw("HIST SAME")

        # Clean up legend labels
        clean_name = name.replace("electronGun_pT_", "").replace("_", "-") + " GeV"
        legend.AddEntry(hist, clean_name, "l")

    legend.Draw()
    
    # Add a line at 0.6 to show Pandora's threshold
    if "profile_discrepancy" in output_path.lower():
        line = ROOT.TLine(0.6, c.GetUymin(), 0.6, c.GetUymax())
        line.SetLineColor(ROOT.kBlack)
        line.SetLineStyle(2)  # Dashed line
        line.SetLineWidth(2)
        line.Draw()
        
        # Add text label for the line
        text = ROOT.TText(0.62, c.GetUymax() * 0.7, "Pandora threshold = 0.6")
        text.SetTextSize(0.03)
        text.Draw()
    
    c.SaveAs(output_path)
    c.Close()

# Create a reader object
reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()
reader.setReadCollectionNames(["MCParticle", "SiTracks", "AllTracks", "PandoraPFOs", "SeedTracks", "PandoraClusters",
                              "EcalBarrelCollectionRec", "EcalEndcapCollectionRec"])

total_events_processed = 0

# Loop over the different samples
for slice_name in files:
    if not files[slice_name]:  # Skip empty slices
        continue

    print(f"Working on sample {slice_name}")

    file_count = 0
    events_processed_in_slice = 0

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

            if events_processed_in_slice < 3:  # Only first 3 events
                print(f"\n=== EVENT {events_processed_in_slice} DEBUG ===")

            # Enhanced track analysis with proper track-cluster matching using LCIO Track API
            track_momentum = -1
            track_cluster_dR = -1
            track_tlv = TLorentzVector()
            matched_track_found = False

            # Try track collections in priority order
            track_collections = ["SiTracks_Refitted", "SiTracks", "AllTracks", "SeedTracks"]

            for track_collection_name in track_collections:
                try:
                    # Direct collection access
                    trackCollection = event.getCollection(track_collection_name)
                    
                    if trackCollection is None:
                        continue
                
                    n_tracks = trackCollection.getNumberOfElements()
                    if n_tracks == 0:
                        continue
                
                    if events_processed_in_slice < 3:  # Debug output
                        print(f"DEBUG - Found {n_tracks} tracks in {track_collection_name}")

                    best_track = None
                    best_dR = 999.0

                    # Loop through all tracks
                    for i in range(n_tracks):
                        track = trackCollection.getElementAt(i)
                        if track is None:
                            continue
                        
                        try:
                            # FIRST: Try to find convenience methods for momentum
                            track_momentum_mag = -1
                            track_p = None
                            
                            # Check if track has direct momentum methods (some LCIO versions might)
                            momentum_methods = ['getPx', 'getPy', 'getPz', 'getPt', 'getMomentum', 'getP']
                            for method_name in momentum_methods:
                                if hasattr(track, method_name):
                                    if events_processed_in_slice < 3:
                                        print(f"DEBUG - Track {i}: Found method {method_name}!")
                                    try:
                                        if method_name == 'getMomentum':
                                            track_p = track.getMomentum()
                                            if track_p and len(track_p) >= 3:
                                                track_momentum_mag = sqrt(track_p[0]**2 + track_p[1]**2 + track_p[2]**2)
                                                break
                                        elif method_name == 'getPt':
                                            pT = track.getPt()
                                            # Still need pz...
                                            pass
                                    except:
                                        pass
                            
                            # If no direct momentum methods found, use track parameters
                            if track_momentum_mag < 0:
                                # Use LCIO Track API directly to get track parameters
                                d0 = track.getD0()
                                phi0 = track.getPhi()
                                omega = track.getOmega()
                                z0 = track.getZ0()
                                tanLambda = track.getTanLambda()
                                
                                if events_processed_in_slice < 3:
                                    print(f"DEBUG - Track {i} params: d0={d0:.3f}, phi0={phi0:.3f}, omega={omega:.6f}, z0={z0:.3f}, tanLambda={tanLambda:.3f}")
                                
                                # Check for valid omega (curvature)
                                if abs(omega) < 1e-10:  # Avoid division by zero
                                    if events_processed_in_slice < 3:
                                        print(f"DEBUG - Track {i}: Invalid omega ({omega})")
                                    continue
                                
                                # Calculate momentum from track parameters
                                # IMPORTANT: Need to determine the correct units for omega
                                # Let's try multiple unit conversions and see which gives reasonable results
                                
                                B_field = 3.57  # Tesla (typical for ILD)
                                
                                # Try different unit assumptions for omega:
                                # Option 1: omega in 1/mm (original assumption)
                                pT_v1 = 0.3 * B_field / abs(omega)  # GeV/c
                                
                                # Option 2: omega in 1/m (if omega is actually in 1/m)
                                pT_v2 = 0.0003 * B_field / abs(omega)  # GeV/c
                                
                                # Option 3: omega in mm^-1 but different conversion factor
                                pT_v3 = 2.998e-4 * B_field / abs(omega)  # GeV/c (exact speed of light factor)
                                
                                # Option 4: omega might already be in appropriate units
                                pT_v4 = B_field / abs(omega)  # GeV/c
                                
                                if events_processed_in_slice < 3:
                                    print(f"DEBUG - Track {i} pT options: v1={pT_v1:.1f}, v2={pT_v2:.1f}, v3={pT_v3:.1f}, v4={pT_v4:.1f} GeV")
                                
                                # Choose the most reasonable pT (should be in range 0.1-5000 GeV for your samples)
                                if 0.1 <= pT_v2 <= 5000:
                                    pT = pT_v2
                                    if events_processed_in_slice < 3:
                                        print(f"DEBUG - Track {i}: Using pT_v2 = {pT:.3f} GeV")
                                elif 0.1 <= pT_v3 <= 5000:
                                    pT = pT_v3
                                    if events_processed_in_slice < 3:
                                        print(f"DEBUG - Track {i}: Using pT_v3 = {pT:.3f} GeV")
                                elif 0.1 <= pT_v4 <= 5000:
                                    pT = pT_v4
                                    if events_processed_in_slice < 3:
                                        print(f"DEBUG - Track {i}: Using pT_v4 = {pT:.3f} GeV")
                                elif 0.1 <= pT_v1 <= 5000:
                                    pT = pT_v1
                                    if events_processed_in_slice < 3:
                                        print(f"DEBUG - Track {i}: Using pT_v1 = {pT:.3f} GeV")
                                else:
                                    if events_processed_in_slice < 3:
                                        print(f"DEBUG - Track {i}: All pT calculations unrealistic, skipping")
                                    continue
                                
                                # Calculate longitudinal momentum: pz = pT * tan(lambda)
                                pz = pT * tanLambda
                                
                                # Total momentum magnitude
                                track_momentum_mag = sqrt(pT*pT + pz*pz)
                                
                                # Calculate momentum components
                                px = pT * TMath.Cos(phi0)
                                py = pT * TMath.Sin(phi0)
                                track_p = [px, py, pz]

                            if events_processed_in_slice < 3:
                                print(f"DEBUG - Track {i}: pT={pT:.3f}, pz={pz:.3f}, |p|={track_momentum_mag:.3f}")

                            # Skip tracks with unrealistic momentum (adjust range based on your electron gun energies)
                            expected_energy_range = [0.1, 5000]  # GeV - adjust based on your samples
                            if track_momentum_mag < expected_energy_range[0] or track_momentum_mag > expected_energy_range[1]:
                                if events_processed_in_slice < 3:
                                    print(f"DEBUG - Track {i}: Momentum outside expected range ({track_momentum_mag:.3f} GeV)")
                                continue

                            # Create track TLorentzVector
                            electron_mass = 0.000511
                            track_energy = sqrt(track_momentum_mag**2 + electron_mass**2)
                            temp_track_tlv = TLorentzVector()
                            temp_track_tlv.SetPxPyPzE(track_p[0], track_p[1], track_p[2], track_energy)

                            # Check match to truth electron
                            dR_to_truth = temp_track_tlv.DeltaR(tlv_truth)

                            if dR_to_truth < best_dR and dR_to_truth < 0.2:
                                best_track = track
                                best_dR = dR_to_truth
                                track_tlv = temp_track_tlv
                                track_momentum = track_momentum_mag
                                
                                if events_processed_in_slice < 3:
                                    print(f"DEBUG - Track {i}: GOOD MATCH! p={track_momentum:.2f} GeV, dR={dR_to_truth:.3f}")
                            
                        except Exception as e:
                            if events_processed_in_slice < 3:
                                print(f"DEBUG - Error processing track {i}: {e}")
                            continue
                        
                    if best_track is not None:
                        matched_track_found = True
                        hists[slice_name]["track_momentum"].Fill(track_momentum)
                        if events_processed_in_slice < 3:
                            print(f"DEBUG - Found matching track: p={track_momentum:.2f} GeV, dR={best_dR:.3f}")
                        break  # Found a good track, stop looking at other collections
                    
                except Exception as e:
                    if events_processed_in_slice < 3:
                        print(f"DEBUG - Error accessing {track_collection_name}: {e}")
                    continue  # Try next collection

            if not matched_track_found:
                if events_processed_in_slice < 3:
                    print(f"DEBUG - No matching tracks found in event {events_processed_in_slice}")
            else:
                if events_processed_in_slice < 3:
                    print(f"DEBUG - SUCCESS: Found track with p={track_momentum:.2f} GeV")

            cluster_energy = 0.0
            max_energy = 0.0
            ecal_coll = ['EcalBarrelCollectionRec', 'EcalEndcapCollectionRec']
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
                            
                            # Get layer information
                            cellID = int(hit.getCellID0())
                            decoder.setValue(cellID)
                            layer = decoder["layer"].value()
                            
                            if layer not in hit_energies_by_layer:
                                hit_energies_by_layer[layer] = 0.0
                            hit_energies_by_layer[layer] += hit_energy
                            cluster_energy += hit_energy
                            
                            # Get position for RMS calculation
                            pos = hit.getPosition()
                            cluster_hit_positions.append([pos[0], pos[1], pos[2]])
                            cluster_hit_energies.append(hit_energy)
                            
                    if cluster_hit_energies_local:
                        max_energy_local = max(cluster_hit_energies_local)
                        max_cluster_energies.append(max_energy_local)
                        
                except Exception as e:
                    continue
                    
            # Set the overall max energy
            max_energy = max(max_cluster_energies) if max_cluster_energies else 0.0

            # [PFO analysis section - keeping your existing code]
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
                        for i in range(pfoCollection.getNumberOfElements()):
                            pfo = pfoCollection.getElementAt(i)
                            if pfo is not None and abs(pfo.getCharge()) > 0.5:
                                pfo_clusters = pfo.getClusters()
                                if pfo_clusters is not None:
                                    for j in range(pfo_clusters.size()):
                                        cluster = pfo_clusters[j]
                                        if cluster is not None:
                                            pfo_hit_energies_by_layer = {}
                                            hits = cluster.getCalorimeterHits()
                                            if hits is not None:
                                                for k in range(hits.size()):
                                                    hit = hits[k]
                                                    if hit is not None:
                                                        hit_energy = hit.getEnergy()
                                                        if hit_energy > 0:
                                                            cellID = int(hit.getCellID0())
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

                                            if pfo_hit_energies_by_layer:
                                                pfo_shower_start_layer = find_shower_start_layer(pfo_hit_energies_by_layer, threshold=0.1)
                                                if pfo_shower_start_layer >= 0:
                                                    hists[slice_name]["pfo_shower_start_layer"].Fill(pfo_shower_start_layer)
            except:
                pass
                
            # Store longitudinal profile information
            for layer, energy in hit_energies_by_layer.items():
                if layer not in longitudinal_profile[slice_name]:
                    longitudinal_profile[slice_name][layer] = 0.0
                longitudinal_profile[slice_name][layer] += energy

            # Calculate cluster properties and create cluster TLorentzVector
            n_hits_in_cluster = len(cluster_hit_energies)
            n_layers_with_energy = len(hit_energies_by_layer)
            cluster_rms_width = -1
            tlv_cluster = TLorentzVector()  # Initialize outside the if block
            
            if cluster_hit_positions and cluster_hit_energies:
                total_energy_check = sum(cluster_hit_energies)
                if total_energy_check > 0:
                    center_x = sum(pos[0] * energy for pos, energy in zip(cluster_hit_positions, cluster_hit_energies)) / total_energy_check
                    center_y = sum(pos[1] * energy for pos, energy in zip(cluster_hit_positions, cluster_hit_energies)) / total_energy_check
                    center_z = sum(pos[2] * energy for pos, energy in zip(cluster_hit_positions, cluster_hit_energies)) / total_energy_check
                    cluster_center = [center_x, center_y, center_z]

                    cluster_rms_width = calculate_rms_width(cluster_hit_positions, cluster_hit_energies, cluster_center)

                    cluster_r = sqrt(center_x**2 + center_y**2)
                    hists[slice_name]["cluster_r"].Fill(cluster_r)
                    
                    # Create proper TLorentzVector for cluster
                    tlv_cluster.SetPxPyPzE(center_x, center_y, center_z, cluster_energy)

            # Calculate shower start layer
            shower_start_layer = find_shower_start_layer(hit_energies_by_layer, threshold=0.01)

            # Find shower maximum layer
            shower_max_layer = -1
            if hit_energies_by_layer:
                shower_max_layer = max(hit_energies_by_layer.keys(), key=lambda k: hit_energies_by_layer[k])

            # PANDORA-STYLE PROFILE DISCREPANCY CALCULATION (ONLY)
            total_energy = sum(hit_energies_by_layer.values())
            
            # Pandora-style method only
            profile_discrepancy_layer, profile_discrepancy = get_profile_discrepancy_pandora_style(
                hit_energies_by_layer, 
                total_energy, 
                energy=trueE,
                X0_per_layer=0.6286, 
                num_layers=50
            )

            # Debug output every 100 events (simplified)
            if events_processed_in_slice % 100 == 0 and hit_energies_by_layer:
                print(f"Event {events_processed_in_slice}: Energy={trueE:.2f} GeV, Profile Discrepancy: {profile_discrepancy:.4f}")

            # TRACK-CLUSTER MATCHING AND E/P ANALYSIS
            e_over_p = -1
            e_minus_p_over_p = -999
            abs_e_minus_p_over_p = -1
            
            if matched_track_found and tlv_cluster.E() > 0:
                # Check if track matches cluster
                track_cluster_dR = track_tlv.DeltaR(tlv_cluster)
                
                if match_track_to_cluster(track_tlv, tlv_cluster, dR_cut=0.2):  # Increased from 0.1 to 0.2
                    # We have a matched track-cluster pair!
                    
                    # Fill matching histograms
                    hists[slice_name]["track_cluster_dR"].Fill(track_cluster_dR)
                    hists[slice_name]["matched_cluster_energy"].Fill(cluster_energy)
                    
                    # Calculate E/p metrics
                    e_over_p, e_minus_p_over_p, abs_e_minus_p_over_p = calculate_e_over_p_metrics(cluster_energy, track_momentum)
                    
                    if e_over_p > 0:  # Valid E/p calculation
                        # Fill E/p histograms
                        hists[slice_name]["electron_E_over_p"].Fill(e_over_p)
                        hists[slice_name]["electron_E_minus_p_over_p"].Fill(e_minus_p_over_p)
                        hists[slice_name]["electron_abs_E_minus_p_over_p"].Fill(abs_e_minus_p_over_p)
                        
                        # Fill 2D correlations
                        hists2d[slice_name]["E_vs_p"].Fill(track_momentum, cluster_energy)
                        hists2d[slice_name]["E_over_p_vs_energy"].Fill(trueE, e_over_p)
                        hists2d[slice_name]["E_over_p_vs_profile_discrepancy"].Fill(e_over_p, profile_discrepancy)
                        
                        print(f"DEBUG - Event {events_processed_in_slice}: E/p SUCCESS! Track p={track_momentum:.2f}, Cluster E={cluster_energy:.2f}, E/p={e_over_p:.3f}, |E-p|/p={abs_e_minus_p_over_p:.3f}")
                        
                    # Update debug output to include E/p info
                    if events_processed_in_slice % 100 == 0:
                        if e_over_p > 0:
                            print(f"  Track p={track_momentum:.2f} GeV, Cluster E={cluster_energy:.2f} GeV, E/p={e_over_p:.3f}, |E-p|/p={abs_e_minus_p_over_p:.3f}")
                else:
                    print(f"DEBUG - Event {events_processed_in_slice}: Track-cluster dR too large: {track_cluster_dR:.3f}")
            elif matched_track_found:
                print(f"DEBUG - Event {events_processed_in_slice}: Track found but no valid cluster (cluster E={cluster_energy:.2f})")
            else:
                if events_processed_in_slice % 100 == 0:
                    print(f"DEBUG - Event {events_processed_in_slice}: No track found for this event")

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

                # Create cluster TLorentzVector and fill kinematic variables
                if tlv_cluster.E() > 0:  # Check if we have a valid cluster
                    # Use TLorentzVector methods directly
                    hists[slice_name]["electron_theta"].Fill(tlv_cluster.Theta())
                    hists[slice_name]["electron_phi"].Fill(tlv_cluster.Phi())
                    hists[slice_name]["electron_eta"].Fill(tlv_cluster.Eta())
                    hists[slice_name]["electron_pt"].Fill(tlv_cluster.Perp())

                    # Fill 2D cluster vs truth
                    hists2d[slice_name]["cluster_E_v_mcp_E"].Fill(cluster_energy, trueE)
                    hists2d[slice_name]["cluster_eta_v_mcp_eta"].Fill(tlv_cluster.Eta(), tlv_truth.Eta())

                # Fill 2D parameter correlations
                hists2d[slice_name]["shower_start_layer_v_max_cell_energy"].Fill(shower_start_layer if shower_start_layer > 0 else 0, max_energy)
                hists2d[slice_name]["cluster_rms_width_v_n_hits"].Fill(cluster_rms_width if cluster_rms_width > 0 else 0, n_hits_in_cluster)
                hists2d[slice_name]["shower_start_layer_v_profile_discrepancy"].Fill(shower_start_layer if shower_start_layer > 0 else 0, profile_discrepancy)

                # Check if this is a matched electron (using proper TLorentzVector matching)
                if tlv_cluster.E() > 0 and isMatched(tlv_truth, tlv_cluster):
                    hists[slice_name]["electron_match_E"].Fill(cluster_energy)
                    hists[slice_name]["electron_match_pt"].Fill(tlv_cluster.Perp())

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

# Custom plotting function adapted from plotHelper style
def plotDiscrepancyHistograms(hist_dict, output_path, x_title, y_title, threshold_line=None, with_bib=False):
    """Plot discrepancy histograms in the style of plotHelper efficiency plots"""
    if not hist_dict:
        return

    can = ROOT.TCanvas("can", "can", 800, 600)
    can.SetLeftMargin(0.12)
    can.SetBottomMargin(0.12)
    can.SetTopMargin(0.08)
    can.SetRightMargin(0.05)
    can.SetLogy()

    # Color scheme like your efficiency plots
    colors = [ROOT.kBlue, ROOT.kRed, ROOT.kGreen+2, ROOT.kMagenta, ROOT.kOrange+1, ROOT.kCyan+1]
    markers = [20, 21, 22, 23, 24, 25]  # Different marker styles

    # Style the histograms
    for i, (name, hist) in enumerate(hist_dict.items()):
        if hist.GetEntries() == 0:
            continue
        
        hist.SetMarkerColor(colors[i % len(colors)])
        hist.SetLineColor(colors[i % len(colors)])
        hist.SetMarkerStyle(markers[i % len(markers)])
        hist.SetMarkerSize(1.2)
        hist.SetLineWidth(2)

    first = True
    for i, (name, hist) in enumerate(hist_dict.items()):
        if hist.GetEntries() == 0:
            continue

        if first:
            hist.Draw("HIST")
            hist.SetTitle("")  # Remove title for cleaner look
            hist.GetXaxis().SetTitle(x_title)
            hist.GetYaxis().SetTitle(y_title)
            hist.GetXaxis().SetTitleSize(0.045)
            hist.GetYaxis().SetTitleSize(0.045)
            hist.GetXaxis().SetLabelSize(0.04)
            hist.GetYaxis().SetLabelSize(0.04)
            hist.GetXaxis().SetTickLength(0.02)
            hist.GetYaxis().SetTickLength(0.02)
            first = False
        else:
            hist.Draw("HIST SAME")

    # Add threshold line if specified
    if threshold_line is not None:
        line = ROOT.TLine(threshold_line, can.GetUymin(), threshold_line, can.GetUymax())
        line.SetLineColor(ROOT.kGray+1)
        line.SetLineStyle(2)
        line.SetLineWidth(2)
        line.Draw()

    # Create legend with header like in your desired style
    leg = ROOT.TLegend(0.65, 0.15, 0.93, 0.45)
    leg.SetBorderSize(1)
    leg.SetFillColor(0)
    leg.SetFillStyle(1001)
    leg.SetTextSize(0.035)
    
    for name, hist in hist_dict.items():
        if hist.GetEntries() == 0:
            continue
        # Clean up the legend entries to match your style
        clean_label = name.replace("electronGun_pT_", "").replace("_", "-") + " GeV"
        leg.AddEntry(hist, clean_label, "l")
    
    leg.Draw()

    # Add text labels in your style
    text = ROOT.TLatex()
    text.SetNDC()
    text.SetTextAlign(11)  # Left aligned
    text.SetTextSize(0.04)
    text.SetTextFont(62)  # Bold font like in image
    
    concept_text = ROOT.TLatex()
    concept_text.SetNDC()
    concept_text.SetTextAlign(31)  # Right aligned
    concept_text.SetTextSize(0.035)
    concept_text.SetTextFont(42)
    concept_text.DrawLatex(0.94, 0.92, "MAIA Detector Concept")
    
    text.DrawLatex(0.15, 0.92, "Muon Collider")
    
    bib_text = "Simulation, with BIB" if with_bib else "Simulation, no BIB"
    text.SetTextFont(42)  # Regular font
    text.DrawLatex(0.15, 0.88, bib_text)
    
    # Add eta cut information
    text.DrawLatex(0.15, 0.84, "|#eta| < 2.4")
    
    # Add threshold information if present
    if threshold_line is not None:
        if "profile" in output_path.lower():
            text.DrawLatex(0.15, 0.80, f"Pandora threshold = {threshold_line}")
        elif "E_minus_p" in output_path:
            text.DrawLatex(0.15, 0.80, f"Steering threshold = {threshold_line}")

    ROOT.gStyle.SetOptStat(0)  # Remove stats box
    can.SaveAs(output_path)
    can.Close()
    return can

# Plot E/p analysis histograms
print("Creating E/p analysis plots...")

# E/p distributions for all pT slices
e_over_p_hists = {}
e_minus_p_over_p_hists = {}
abs_e_minus_p_over_p_hists = {}

for s in hists:
    if "electron_E_over_p" in hists[s] and hists[s]["electron_E_over_p"].GetEntries() > 0:
        e_over_p_hists[s] = hists[s]["electron_E_over_p"]
    if "electron_E_minus_p_over_p" in hists[s] and hists[s]["electron_E_minus_p_over_p"].GetEntries() > 0:
        e_minus_p_over_p_hists[s] = hists[s]["electron_E_minus_p_over_p"]
    if "electron_abs_E_minus_p_over_p" in hists[s] and hists[s]["electron_abs_E_minus_p_over_p"].GetEntries() > 0:
        abs_e_minus_p_over_p_hists[s] = hists[s]["electron_abs_E_minus_p_over_p"]

# Create the combined plots using the new style
if e_over_p_hists:
    plotDiscrepancyHistograms(e_over_p_hists, "plots/E_over_p_all_pt_slices.png", "E/p", "Entries", threshold_line=1.0)

if abs_e_minus_p_over_p_hists:
    plotDiscrepancyHistograms(abs_e_minus_p_over_p_hists, "plots/abs_E_minus_p_over_p_all_pt_slices.png", "|E-p|/p", "Entries", threshold_line=0.2)

if e_minus_p_over_p_hists:
    plotDiscrepancyHistograms(e_minus_p_over_p_hists, "plots/E_minus_p_over_p_all_pt_slices.png", "(E-p)/p", "Entries")

# Add this before the plotting section:
print("Creating combined profile discrepancy plot...")
profile_discrepancy_hists = {}

for s in hists:
    if "electron_profile_discrepancy" in hists[s] and hists[s]["electron_profile_discrepancy"].GetEntries() > 0:
        profile_discrepancy_hists[s] = hists[s]["electron_profile_discrepancy"]

# Also update the profile discrepancy plot to use the new style
if profile_discrepancy_hists:
    plotDiscrepancyHistograms(profile_discrepancy_hists, "plots/profile_discrepancy_all_pt_slices.png",
                  "Max Profile Discrepancy", "Entries", threshold_line=0.6)

# Individual plots for each energy range with detailed statistics
for s in abs_e_minus_p_over_p_hists:
    hist = abs_e_minus_p_over_p_hists[s]
    
    c = ROOT.TCanvas("can", "can", 800, 600)
    c.SetLogy()
    
    hist.SetLineColor(ROOT.kBlue)
    hist.SetLineWidth(3)
    hist.GetXaxis().SetTitle("|E-p|/p")
    hist.GetYaxis().SetTitle("Entries")
    hist.SetTitle(f"|E-p|/p Distribution - {s.replace('electronGun_pT_', '').replace('_', '-')} GeV")
    hist.Draw("HIST")
    
    # Add steering file threshold line at 0.2
    y_max = hist.GetMaximum()
    line = ROOT.TLine(0.2, 1, 0.2, y_max)
    line.SetLineColor(ROOT.kRed)
    line.SetLineStyle(2)
    line.SetLineWidth(2)
    line.Draw()
    
    # Add statistics and threshold info
    mean_val = hist.GetMean()
    rms_val = hist.GetRMS()
    entries = hist.GetEntries()
    
    # Calculate percentage above 0.2 threshold
    above_threshold = 0
    for i in range(1, hist.GetNbinsX() + 1):
        bin_center = hist.GetBinCenter(i)
        if bin_center > 0.2:
            above_threshold += hist.GetBinContent(i)
    
    pct_above_threshold = (above_threshold / entries * 100) if entries > 0 else 0
    
    # Add text annotations
    text1 = ROOT.TText(0.5, y_max * 0.8, f"Entries: {int(entries)}")
    text1.SetTextSize(0.03)
    text1.Draw()
    
    text2 = ROOT.TText(0.5, y_max * 0.75, f"Mean: {mean_val:.3f}")
    text2.SetTextSize(0.03)
    text2.Draw()
    
    text3 = ROOT.TText(0.5, y_max * 0.7, f"% > 0.2: {pct_above_threshold:.1f}%")
    text3.SetTextSize(0.03)
    text3.SetTextColor(ROOT.kRed)
    text3.Draw()
    
    text4 = ROOT.TText(0.22, y_max * 0.6, "Steering threshold = 0.2")
    text4.SetTextSize(0.03)
    text4.SetTextColor(ROOT.kRed)
    text4.Draw()
    
    c.SaveAs(f"plots/abs_E_minus_p_over_p_{s}.png")
    c.Close()

# Plot LCElectronId parameter histograms
lcelectronid_hists = {}
for param in ["shower_start_layer", "max_cell_energy", "cluster_cone_energy"]:
    lcelectronid_hists[param] = {}
    for s in hists:
        if f"electron_{param}" in hists[s]:
            lcelectronid_hists[param][s] = hists[s][f"electron_{param}"]

    if lcelectronid_hists[param]:
        plotHistograms(lcelectronid_hists[param], f"plots/lcelectronid_{param}.png",
                      param.replace("_", " ").title(), "Entries")

# Plot Profile Discrepancy for all pT slices on the same plot
print("Creating combined profile discrepancy plot...")
profile_discrepancy_hists = {}

for s in hists:
    if "electron_profile_discrepancy" in hists[s] and hists[s]["electron_profile_discrepancy"].GetEntries() > 0:
        profile_discrepancy_hists[s] = hists[s]["electron_profile_discrepancy"]

if profile_discrepancy_hists:
    plotHistograms(profile_discrepancy_hists, "plots/profile_discrepancy_all_pt_slices.png",
                  "Max Profile Discrepancy", "Entries")

# Also create individual plots for each pT slice
for s, hist in profile_discrepancy_hists.items():
    c = ROOT.TCanvas("can", "can", 800, 600)
    c.SetLogy()
    
    hist.SetLineColor(ROOT.kBlue)
    hist.SetLineWidth(3)
    hist.GetXaxis().SetTitle("Max Profile Discrepancy")
    hist.GetYaxis().SetTitle("Entries")
    hist.SetTitle(f"Profile Discrepancy - {s.replace('electronGun_pT_', '').replace('_', '-')} GeV")
    hist.Draw("HIST")
    
    # Add Pandora threshold line
    y_max = hist.GetMaximum()
    line = ROOT.TLine(0.6, 1, 0.6, y_max)
    line.SetLineColor(ROOT.kRed)
    line.SetLineStyle(2)
    line.SetLineWidth(2)
    line.Draw()
    
    # Add text label
    text = ROOT.TText(0.62, y_max * 0.7, "Pandora threshold = 0.6")
    text.SetTextSize(0.03)
    text.SetTextColor(ROOT.kRed)
    text.Draw()
    
    # Add statistics
    mean_val = hist.GetMean()
    rms_val = hist.GetRMS()
    entries = hist.GetEntries()
    
    stats_text = ROOT.TText(0.15, y_max * 0.8, f"Entries: {int(entries)}")
    stats_text.SetTextSize(0.03)
    stats_text.Draw()
    
    stats_text2 = ROOT.TText(0.15, y_max * 0.75, f"Mean: {mean_val:.3f}")
    stats_text2.Draw()
    
    stats_text3 = ROOT.TText(0.15, y_max * 0.7, f"RMS: {rms_val:.3f}")
    stats_text3.SetTextSize(0.03)
    stats_text3.Draw()
    
    c.SaveAs(f"plots/profile_discrepancy_{s}.png")
    c.Close()

# Create summary statistics table
print("\n" + "="*60)
print("PROFILE DISCREPANCY SUMMARY STATISTICS")
print("="*60)
print(f"{'pT Range (GeV)':<15} {'Entries':<8} {'Mean':<8} {'RMS':<8} {'% > 0.6':<8}")
print("-"*60)

for s, hist in profile_discrepancy_hists.items():
    pt_range = s.replace('electronGun_pT_', '').replace('_', '-')
    entries = int(hist.GetEntries())
    mean_val = hist.GetMean()
    rms_val = hist.GetRMS()
    
    # Calculate percentage above 0.6 threshold
    total_entries = hist.GetEntries()
    above_threshold = 0
    for i in range(1, hist.GetNbinsX() + 1):
        bin_center = hist.GetBinCenter(i)
        if bin_center > 0.6:
            above_threshold += hist.GetBinContent(i)
    
    pct_above_threshold = (above_threshold / total_entries * 100) if total_entries > 0 else 0
    
    print(f"{pt_range:<15} {entries:<8} {mean_val:<8.3f} {rms_val:<8.3f} {pct_above_threshold:<8.1f}")

print("-"*60)

# Print E/p summary statistics
print("\n" + "="*70)
print("E/p ANALYSIS SUMMARY STATISTICS")
print("="*70)
print(f"{'pT Range (GeV)':<15} {'Entries':<8} {'Mean |E-p|/p':<12} {'RMS':<8} {'% > 0.2':<8}")
print("-"*70)

for s in abs_e_minus_p_over_p_hists:
    hist = abs_e_minus_p_over_p_hists[s]
    pt_range = s.replace('electronGun_pT_', '').replace('_', '-')
    entries = int(hist.GetEntries())
    mean_val = hist.GetMean()
    rms_val = hist.GetRMS()
    
    # Calculate percentage above 0.2 threshold
    above_threshold = 0
    for i in range(1, hist.GetNbinsX() + 1):
        bin_center = hist.GetBinCenter(i)
        if bin_center > 0.2:
            above_threshold += hist.GetBinContent(i)
    
    pct_above_threshold = (above_threshold / entries * 100) if entries > 0 else 0
    
    print(f"{pt_range:<15} {entries:<8} {mean_val:<12.3f} {rms_val:<8.3f} {pct_above_threshold:<8.1f}")


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
print("INDIVIDUAL PLOTS:")
for s in profile_discrepancy_hists:
    pt_range = s.replace('electronGun_pT_', '').replace('_', '-')
    print(f"  plots/profile_discrepancy_{s}.png ({pt_range} GeV)")
print("LONGITUDINAL PROFILES:")
for s in longitudinal_profile:
    if longitudinal_profile[s]:
        pt_range = s.replace('electronGun_pT_', '').replace('_', '-')
        print(f"  plots/longitudinal_profile_{s}.png ({pt_range} GeV)")

