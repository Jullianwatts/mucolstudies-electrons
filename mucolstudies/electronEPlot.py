from math import exp, gamma, log, sqrt 
import glob
import ROOT
import pyLCIO
from pyLCIO import EVENT, UTIL
from ROOT import TH1F, TH2F, TFile, TLorentzVector, TMath, TCanvas
import numpy as np
import math

ROOT.gROOT.SetBatch()
max_events = 10
import os
samples = glob.glob("/data/fmeloni/DataMuC_MAIA_v0/v5/reco/pionGun*")
#samples = glob.glob("/data/fmeloni/DataMuC_MAIA_v0/v5/reco/electronGun*")
files = {}
slices = ["0_50", "50_250", "250_1000", "1000_5000"]
for s in slices:
    files[f"pionGun_pT_{s}"] = []

for s in samples:
    sname = s.split("/")[-1]
    if sname not in files:
        continue
    files[sname] = glob.glob(f"{s}/*.slcio")

print("Found files:")

# Initialize longitudinal profile dictionary
longitudinal_profile = {slice_name: {} for slice_name in files}

# Define histogram variables for LCElectronId parameters (using pion key set)
variables = {
    "pion": {
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
    for obj in ["mcp", "mcp_pion", "pion", "pion_match"]:
        for var in variables["pion"]:
            hists[s][f"{obj}_{var}"] = ROOT.TH1F(f"{s}_{obj}_{var}", f"{s}_{obj}_{var}",
                                              variables["pion"][var]["nbins"],
                                              variables["pion"][var]["xmin"],
                                              variables["pion"][var]["xmax"])

    # Special cluster histograms
    hists[s]["cluster_nhits"] = ROOT.TH1F(f"{s}_cluster_nhits", f"{s}_cluster_nhits", 50, 0, 100)
    hists[s]["cluster_r"] = ROOT.TH1F(f"{s}_cluster_r", f"{s}_cluster_r", 50, 0, 3000)
    hists[s]["pfo_shower_start_layer"] = ROOT.TH1F(f"{s}_pfo_shower_start_layer", f"{s}_pfo_shower_start_layer", 20, 0, 20)

    # NEW: book per-pion shower start layer histogram so we actually plot it later
    hists[s]["pion_shower_start_layer"] = ROOT.TH1F(f"{s}_pion_shower_start_layer", f"{s}_pion_shower_start_layer", 20, 0, 20)
    
    # LCElectronId parameter histograms - only Pandora method
    #hists[s]["lcelectronid_max_inner_layer"] = ROOT.TH1F(f"{s}_lcelectronid_max_inner_layer", f"{s}_lcelectronid_max_inner_layer", 10, 0, 10)
    #hists[s]["lcelectronid_max_energy"] = ROOT.TH1F(f"{s}_lcelectronid_max_energy", f"{s}_lcelectronid_max_energy", 50, 0, 20)
    #hists[s]["lcelectronid_max_profile_start"] = ROOT.TH1F(f"{s}_lcelectronid_max_profile_start", f"{s}_lcelectronid_max_profile_start", 20, 0, 20)
    hists[s]["lcpionid_max_profile_discrepancy"] = ROOT.TH1F(f"{s}_lcpionid_max_profile_discrepancy", f"{s}_lcpionid_max_profile_discrepancy", 100, 0, 2)
    hists[s]["lcelectronid_max_residual_e_over_p"] = ROOT.TH1F(f"{s}_lcelectronid_max_residual_e_over_p", f"{s}_lcelectronid_max_residual_e_over_p", 50, 0, 1)

# PANDORA-STYLE PROFILE DISCREPANCY FUNCTION (ONLY)
def get_profile_discrepancy_pandora_style(energy_by_layer, total_energy, energy=None, num_layers=50, X0_per_layer=0.6286):
    if total_energy == 0 or len(energy_by_layer) < 3:
        return -1, 0.0
    
    if energy is None:
        energy = total_energy
    
    # Find shower maximum layer
    max_layer = max(energy_by_layer.keys(), key=lambda k: energy_by_layer[k])
    max_energy_fraction = energy_by_layer[max_layer] / total_energy
    
    # Pandora-style shower profile analysis
    start_layer = max(0, max_layer - 5)
    end_layer = min(num_layers, max_layer + 10)
    
    total_discrepancy = 0.0
    n_layers_compared = 0
    max_discrepancy = 0.0
    
    for layer in range(start_layer, end_layer):
        if layer in energy_by_layer:
            f_obs = energy_by_layer[layer] / total_energy
            distance_from_max = abs(layer - max_layer)
            if distance_from_max == 0:
                f_exp = max_energy_fraction
            else:
                decay_constant = 3.0  # Tunable
                f_exp = max_energy_fraction * exp(-distance_from_max / decay_constant)
            discrepancy = abs(f_obs - f_exp) / (f_exp + 0.01)
            total_discrepancy += discrepancy
            n_layers_compared += 1
            if discrepancy > max_discrepancy:
                max_discrepancy = discrepancy
    
    avg_discrepancy = (total_discrepancy / n_layers_compared) if n_layers_compared > 0 else 0.0
    scaled_discrepancy = avg_discrepancy * 0.5  # Scaling factor to be ~Pandora scale
    return max_layer, scaled_discrepancy

# Matching function
def isMatched(tlv1, tlv2, dR_cut=0.1):
    if tlv1.DeltaR(tlv2) < dR_cut:
        return True
    return False

def find_shower_start_layer_absolute(energy_by_layer, absolute_threshold_gev=1.0):
    if not energy_by_layer:
        return -1
    sorted_layers = sorted(energy_by_layer.keys())
    for layer in sorted_layers:
        energy_gev = energy_by_layer[layer]
        if energy_gev >= absolute_threshold_gev:
            return layer
    for layer in sorted_layers:
        if energy_by_layer[layer] > 0:
            return layer
    return -1

#  old function for comparison -- dont know if correct
def find_shower_start_layer_percentage(energy_by_layer, threshold=0.01):
    total_energy = sum(energy_by_layer.values())
    if total_energy == 0:
        return -1
    sorted_layers = sorted(energy_by_layer.keys())
    for layer in sorted_layers:
        fraction = energy_by_layer[layer] / total_energy
        if fraction >= threshold:
            return layer
    for layer in sorted_layers:
        if energy_by_layer[layer] > 0:
            return layer
    return -1

def generate_expected_em_profile(num_layers=60, energy=1.0, X0_per_layer=0.6286):
    E_c = 0.015  # GeV (critical energy)
    if energy > E_c:
        t_max = log(energy / E_c) - 0.5
    else:
        t_max = 0.0
    a = 1.0 + t_max
    b = 1.0
    expected_profile = {}
    total_expected = 0.0
    for layer in range(num_layers):
        t = layer * X0_per_layer
        if t == 0:
            profile_value = 0.0 if a > 1 else 1.0 / gamma(a)
        else:
            try:
                profile_value = (b**a / gamma(a)) * (t**(a-1)) * exp(-b*t)
            except (OverflowError, ZeroDivisionError, ValueError):
                profile_value = 0.0
        expected_profile[layer] = profile_value
        total_expected += profile_value
    if total_expected > 0:
        for layer in expected_profile:
            expected_profile[layer] /= total_expected
    return expected_profile

def plotHistograms(hist_dict, output_path, x_title, y_title):
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
        clean_name = (name.split("_pT_")[-1]).replace("_", "-") + " GeV"
        legend.AddEntry(hist, clean_name, "l")
    legend.Draw()
    if "profile_discrepancy" in output_path.lower():
        line = ROOT.TLine(0.6, c.GetUymin(), 0.6, c.GetUymax())
        line.SetLineColor(ROOT.kBlack)
        line.SetLineStyle(2)
        line.SetLineWidth(2)
        line.Draw()
        text = ROOT.TText(0.62, c.GetUymax() * 0.7, "Pandora threshold = 0.6")
        text.SetTextSize(0.03)
        text.Draw()
    c.SaveAs(output_path)
    c.Close()

def plotDiscrepancyHistograms(hist_dict, output_path, x_title, y_title, threshold_line=None, with_bib=False):
    can = ROOT.TCanvas("can", "can", 800, 600)
    can.SetLeftMargin(0.12)
    can.SetBottomMargin(0.12)
    can.SetTopMargin(0.08)
    can.SetRightMargin(0.05)
    can.SetLogy()
    colors = [ROOT.kBlue, ROOT.kRed, ROOT.kGreen+2, ROOT.kMagenta, ROOT.kOrange+1, ROOT.kCyan+1]
    markers = [20, 21, 22, 23, 24, 25]
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
            hist.SetTitle("")
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
    if threshold_line is not None:
        line = ROOT.TLine(threshold_line, can.GetUymin(), threshold_line, can.GetUymax())
        line.SetLineColor(ROOT.kGray+1)
        line.SetLineStyle(2)
        line.SetLineWidth(2)
        line.Draw()
    leg = ROOT.TLegend(0.65, 0.15, 0.93, 0.45)
    leg.SetBorderSize(1)
    leg.SetFillColor(0)
    leg.SetFillStyle(1001)
    leg.SetTextSize(0.035)
    for name, hist in hist_dict.items():
        if hist.GetEntries() == 0:
            continue
        clean_label = (name.split("_pT_")[-1]).replace("_", "-") + " GeV"
        leg.AddEntry(hist, clean_label, "l")
    leg.Draw()
    text = ROOT.TLatex()
    text.SetNDC()
    text.SetTextAlign(11)
    text.SetTextSize(0.04)
    text.SetTextFont(62)
    concept_text = ROOT.TLatex()
    concept_text.SetNDC()
    concept_text.SetTextAlign(31)
    concept_text.SetTextSize(0.035)
    concept_text.SetTextFont(42)
    concept_text.DrawLatex(0.94, 0.92, "MAIA Detector Concept")
    text.DrawLatex(0.15, 0.92, "Muon Collider")
    bib_text = "Simulation, with BIB" if with_bib else "Simulation, no BIB"
    text.SetTextFont(42)
    text.DrawLatex(0.15, 0.88, bib_text)
    text.DrawLatex(0.15, 0.84, "|#eta| < 2.4")
    if threshold_line is not None:
        if "profile" in output_path.lower():
            text.DrawLatex(0.15, 0.80, f"Pandora threshold = {threshold_line}")
        elif "E_minus_p" in output_path:
            text.DrawLatex(0.15, 0.80, f"Steering threshold = {threshold_line}")
    ROOT.gStyle.SetOptStat(0)
    can.SaveAs(output_path)
    can.Close()
    return can

# Create a reader object
reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()
reader.setReadCollectionNames(["MCParticle", "SiTracks", "AllTracks", "PandoraPFOs", "SiTracks_refitted", "PandoraClusters",
                              "EcalBarrelCollectionRec", "EcalEndcapCollectionRec"])

total_events_processed = 0

# Loop over the different samples
for slice_name in files:
    if not files[slice_name]:  # Skip empty slices
        print("no files")
        continue

    print(f"Working on sample {slice_name}")

    file_count = 0
    events_processed_in_slice = 0

    # Loop over the files in a sample
    for f in files[slice_name]:
        if max_events > 0 and file_count >= max_events:
            break

        reader.open(f)  # <<< OPEN THE FILE BEFORE ITERATING

        for ievt, event in enumerate(reader):

            # Reset variables for this event
            hit_energies_by_layer = {}
            cluster_hit_positions = []
            cluster_hit_energies = []

            # Get truth particle
            try:
                mcpCollection = event.getCollection('MCParticle')
                truePion = mcpCollection[0]
            except:
                continue

            # Check if it's actually a pion (PDG ±211)
            if abs(truePion.getPDG()) != 211:
                continue

            # Get truth energy and momentum
            trueE = truePion.getEnergy()
            dp3 = truePion.getMomentum()
            tlv_truth = TLorentzVector()
            tlv_truth.SetPxPyPzE(dp3[0], dp3[1], dp3[2], trueE)
            if abs(tlv_truth.Eta()) > 2.4:
                continue

            # Fill truth histograms
            hists[slice_name]["mcp_E"].Fill(trueE)
            hists[slice_name]["mcp_pt"].Fill(tlv_truth.Perp())
            hists[slice_name]["mcp_eta"].Fill(tlv_truth.Eta())
            hists[slice_name]["mcp_phi"].Fill(tlv_truth.Phi())
            hists[slice_name]["mcp_theta"].Fill(tlv_truth.Theta())

            track_momentum = -1
            track_cluster_dR = -1
            track_tlv = TLorentzVector()
            matched_track_found = False

            track_collections = ["SiTracks_Refitted", "SiTracks", "AllTracks", "SeedTracks"]

            for track_collection_name in track_collections:
                try:
                    trackCollection = event.getCollection(track_collection_name)
                    if trackCollection is None:
                        continue
                    n_tracks = trackCollection.getNumberOfElements()
                    if n_tracks == 0:
                        continue
                    best_track = None
                    best_dR = 999.0

                    # Loop through all tracks
                    for i in range(n_tracks):
                        track = trackCollection.getElementAt(i)
                        if track is None:
                            continue
                        try:
                            track_momentum_mag = -1
                            track_p = None
                            pT = None
                            momentum_methods = ['getPx', 'getPy', 'getPz', 'getPt', 'getMomentum', 'getP']
                            for method_name in momentum_methods:
                                if hasattr(track, method_name):
                                    try:
                                        if method_name == 'getMomentum':
                                            track_p = track.getMomentum()
                                            if track_p and len(track_p) >= 3:
                                                track_momentum_mag = sqrt(track_p[0]**2 + track_p[1]**2 + track_p[2]**2)
                                                break
                                        elif method_name == 'getPt':
                                            pT = track.getPt()
                                    except:
                                        pass
                            if track_momentum_mag < 0:
                                d0 = track.getD0()
                                phi0 = track.getPhi()
                                omega = track.getOmega()
                                z0 = track.getZ0()
                                tanLambda = track.getTanLambda()
                                if pT is None:
                                    if abs(omega) > 0:
                                        pT = abs(1.0/omega)
                                    else:
                                        continue
                                pz = pT * tanLambda
                                track_momentum_mag = sqrt(pT*pT + pz*pz)
                                px = pT * TMath.Cos(phi0)
                                py = pT * TMath.Sin(phi0)
                                track_p = [px, py, pz]

                            expected_energy_range = [0.1, 5000]
                            if track_momentum_mag < expected_energy_range[0] or track_momentum_mag > expected_energy_range[1]:
                                continue

                            # Create track TLorentzVector (use pion mass)
                            pion_mass = 0.13957
                            track_energy = sqrt(track_momentum_mag**2 + pion_mass**2)
                            temp_track_tlv = TLorentzVector()
                            temp_track_tlv.SetPxPyPzE(track_p[0], track_p[1], track_p[2], track_energy)

                            # Check match to truth pion
                            dR_to_truth = temp_track_tlv.DeltaR(tlv_truth)
                            if dR_to_truth < best_dR and dR_to_truth < 0.1:
                                best_track = track
                                best_dR = dR_to_truth
                                track_tlv = temp_track_tlv
                                track_momentum = track_momentum_mag
                        except Exception as e:
                            continue
                    if best_track is not None:
                        matched_track_found = True
                        hists[slice_name]["track_momentum"].Fill(track_momentum)
                        break
                except Exception as e:
                    continue  # Try next collection

            cluster_energy = 0.0
            max_energy = 0.0
            ecal_coll = ['EcalBarrelCollectionRec', 'EcalEndcapCollectionRec']
            max_cluster_energies = []

            # Use PandoraClusters instead of raw ECAL hits
            try:
                clusterCollection = event.getCollection('PandoraClusters')
                # Find cluster that best matches truth pion
                best_cluster = None
                best_dR = 999.0
                for i in range(clusterCollection.getNumberOfElements()):
                    cluster = clusterCollection.getElementAt(i)
                    if cluster is None:
                        continue
                    cluster_pos = cluster.getPosition()
                    cluster_E = cluster.getEnergy()
                    if cluster_E > 1.0:  # Minimum energy threshold
                        cluster_tlv = TLorentzVector()
                        cluster_tlv.SetPxPyPzE(cluster_pos[0], cluster_pos[1], cluster_pos[2], cluster_E)
                        dR = cluster_tlv.DeltaR(tlv_truth)
                        if dR < best_dR and dR < 0.1:
                            best_cluster = cluster
                            best_dR = dR
                if best_cluster:
                    cluster_energy = best_cluster.getEnergy()
                    hits = best_cluster.getCalorimeterHits()
                    if hits:
                        ecal_decoders = {}
                        ecal_collections = ['EcalBarrelCollectionRec', 'EcalEndcapCollectionRec']
                        for coll_name in ecal_collections:
                            try:
                                temp_collection = event.getCollection(coll_name)
                                if temp_collection:
                                    encoding = temp_collection.getParameters().getStringVal(EVENT.LCIO.CellIDEncoding)
                                    decoder = UTIL.BitField64(encoding)
                                    ecal_decoders[coll_name] = decoder
                            except:
                                continue
                        for j in range(hits.size()):
                            hit = hits[j]
                            if hit and hit.getEnergy() > 0:
                                hit_energy = hit.getEnergy()
                                cellID = int(hit.getCellID0())
                                layer = -1
                                for coll_name, decoder in ecal_decoders.items():
                                    try:
                                        decoder.setValue(cellID)
                                        layer = decoder["layer"].value()
                                        break
                                    except:
                                        continue
                                if layer >= 0:
                                    if layer not in hit_energies_by_layer:
                                        hit_energies_by_layer[layer] = 0.0
                                    hit_energies_by_layer[layer] += hit_energy
                                    pos = hit.getPosition()
                                    cluster_hit_positions.append([pos[0], pos[1], pos[2]])
                                    cluster_hit_energies.append(hit_energy)
                        if cluster_hit_energies:
                            max_energy = max(cluster_hit_energies)
                            max_cluster_energies = [max_energy]  # Single cluster approach
            except Exception as e:
                if events_processed_in_slice < 3:
                    print(f"DEBUG - PandoraClusters failed: {e}")

            # Set the overall max energy
            max_energy = max(max_cluster_energies) if max_cluster_energies else 0.0

            # Determine absolute threshold based on truth energy
            if trueE > 1000:
                abs_threshold = 2.0
            elif trueE > 100:
                abs_threshold = 1.0
            else:
                abs_threshold = 0.5
            
            shower_start_layer_absolute = find_shower_start_layer_absolute(
                hit_energies_by_layer, 
                absolute_threshold_gev=abs_threshold
            )
            # Also calculate with old percentage method for comparison
            shower_start_layer_percentage = find_shower_start_layer_percentage(hit_energies_by_layer, threshold=0.01)
            # Use the absolute method as the primary result
            shower_start_layer = shower_start_layer_absolute

            # (decoders prepared above if needed later)
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

            # Store longitudinal profile information
            for layer, energy in hit_energies_by_layer.items():
                if layer not in longitudinal_profile[slice_name]:
                    longitudinal_profile[slice_name][layer] = 0.0
                longitudinal_profile[slice_name][layer] += energy

            n_hits_in_cluster = len(cluster_hit_energies)
            n_layers_with_energy = len(hit_energies_by_layer)
            cluster_rms_width = -1
            tlv_cluster = TLorentzVector()  
            if cluster_hit_positions and cluster_hit_energies:
                total_energy_check = sum(cluster_hit_energies)
                if total_energy_check > 0:
                    center_x = sum(pos[0] * energy for pos, energy in zip(cluster_hit_positions, cluster_hit_energies)) / total_energy_check
                    center_y = sum(pos[1] * energy for pos, energy in zip(cluster_hit_positions, cluster_hit_energies)) / total_energy_check
                    center_z = sum(pos[2] * energy for pos, energy in zip(cluster_hit_positions, cluster_hit_energies)) / total_energy_check
                    cluster_center = [center_x, center_y, center_z]
                    cluster_r = sqrt(center_x**2 + center_y**2)
                    hists[slice_name]["cluster_r"].Fill(cluster_r)
                    tlv_cluster.SetPxPyPzE(center_x, center_y, center_z, cluster_energy)

            # Find shower maximum layer
            shower_max_layer = -1
            if hit_energies_by_layer:
                shower_max_layer = max(hit_energies_by_layer.keys(), key=lambda k: hit_energies_by_layer[k])

            total_energy = sum(hit_energies_by_layer.values())
            profile_discrepancy_layer, profile_discrepancy = get_profile_discrepancy_pandora_style(
                hit_energies_by_layer, 
                total_energy, 
                energy=trueE,
                X0_per_layer=0.6286, 
                num_layers=50
            )

            if events_processed_in_slice % 100 == 0 and hit_energies_by_layer:
                print(f"Event {events_processed_in_slice}: Energy={trueE:.2f} GeV, Profile Discrepancy: {profile_discrepancy:.4f}, Shower Start: Layer {shower_start_layer} (abs method)")

            # Fill cluster/PFO-level shower start layer hists (NEW to guarantee plots)
            if shower_start_layer >= 0:
                hists[slice_name]["pion_shower_start_layer"].Fill(shower_start_layer)
                hists[slice_name]["pfo_shower_start_layer"].Fill(shower_start_layer)

            # TRACK-CLUSTER MATCHING AND E/P ANALYSIS (left as-is per your request)
            e_over_p = -1
            e_minus_p_over_p = -999
            abs_e_minus_p_over_p = -1
            if matched_track_found and tlv_cluster.E() > 0:
                track_cluster_dR = track_tlv.DeltaR(tlv_cluster)
                pass  # (keeping your commented/guarded blocks unchanged)

            # Fill LCElectronId-style parameters (pion equivalents)
            hists[slice_name]["pion_max_cell_energy"].Fill(max_energy)

            # >>> Guard discrepancy fill so zeros from empty events don’t dominate
            if total_energy > 0:
                hists[slice_name]["pion_profile_discrepancy"].Fill(profile_discrepancy)

            hists[slice_name]["pion_cluster_cone_energy"].Fill(cluster_energy)

            # Fill reconstructed quantities if we have a cluster
            if cluster_energy > 0:
                hists[slice_name]["pion_E"].Fill(cluster_energy)
                hists[slice_name]["pion_cluster_cone_energy"].Fill(cluster_energy)
                if tlv_cluster.E() > 0:
                    hists[slice_name]["pion_theta"].Fill(tlv_cluster.Theta())
                    hists[slice_name]["pion_phi"].Fill(tlv_cluster.Phi())
                    hists[slice_name]["pion_eta"].Fill(tlv_cluster.Eta())
                    hists[slice_name]["pion_pt"].Fill(tlv_cluster.Perp())

            events_processed_in_slice += 1
            total_events_processed += 1

            if events_processed_in_slice % 100 == 0:
                print(f"  Processed {events_processed_in_slice} events in {slice_name}")
                
        reader.close()  # <<< CLOSE AFTER EACH FILE
        file_count += 1

print(f"Total events processed: {total_events_processed}")

# Create longitudinal profile plots
print("Creating longitudinal profile plots...")
for s in longitudinal_profile:
    if not longitudinal_profile[s]:
        continue
    layers = sorted(longitudinal_profile[s].keys())
    energy_vals = [longitudinal_profile[s][l] for l in layers]
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

# E/p distributions blocks left as-is (they’re electron-oriented and currently not being filled)

profile_discrepancy_hists = {}
for s in hists:
    if "pion_profile_discrepancy" in hists[s] and hists[s]["pion_profile_discrepancy"].GetEntries() > 0:
        profile_discrepancy_hists[s] = hists[s]["pion_profile_discrepancy"]

# (Optional) Per-slice plots for |E-p|/p left unchanged, but titles fixed to pionGun if ever used
for s in []:  # no-op: kept structure but disabled since electron hists aren’t filled
    pass

# Gather pion shower start hists (NEW key)
shower_start_hists = {}
for s in hists:
    if "pion_shower_start_layer" in hists[s] and hists[s]["pion_shower_start_layer"].GetEntries() > 0:
        shower_start_hists[s] = hists[s]["pion_shower_start_layer"]

# Combined pion shower start layer plot
if shower_start_hists:
    plotDiscrepancyHistograms(shower_start_hists, "plots/pion_shower_start_layer_all_pt_slices_ABSOLUTE_THRESHOLD.png",
                  "Shower Start Layer (Absolute Energy Threshold)", "Entries")
    
    # Individual shower start layer plots with MAIA formatting
    for s, hist in shower_start_hists.items():
        pt_range = s.replace('pionGun_pT_', '').replace('_', '-')
        can = ROOT.TCanvas("can", "can", 800, 600)
        can.SetLeftMargin(0.12)
        can.SetBottomMargin(0.12)
        can.SetTopMargin(0.08)
        can.SetRightMargin(0.05)
        can.SetLogy()
        hist.SetLineColor(ROOT.kBlue)
        hist.SetLineWidth(3)
        hist.SetTitle("")
        hist.GetXaxis().SetTitle("Shower Start Layer")
        hist.GetYaxis().SetTitle("Entries")
        hist.GetXaxis().SetTitleSize(0.045)
        hist.GetYaxis().SetTitleSize(0.045)
        hist.GetXaxis().SetLabelSize(0.04)
        hist.GetYaxis().SetLabelSize(0.04)
        hist.Draw("HIST")
        text = ROOT.TLatex()
        text.SetNDC()
        text.SetTextAlign(11)
        text.SetTextSize(0.04)
        text.SetTextFont(62)
        concept_text = ROOT.TLatex()
        concept_text.SetNDC()
        concept_text.SetTextAlign(31)
        concept_text.SetTextSize(0.035)
        concept_text.SetTextFont(42)
        concept_text.DrawLatex(0.94, 0.92, "MAIA Detector Concept")
        text.DrawLatex(0.15, 0.92, "Muon Collider")
        text.SetTextFont(42)
        text.DrawLatex(0.15, 0.88, "Simulation, no BIB")
        text.DrawLatex(0.15, 0.84, "|#eta| < 2.4")
        text.DrawLatex(0.15, 0.80, f"pT: {pt_range} GeV")
        text.DrawLatex(0.15, 0.76, "FIXED: Absolute Energy Threshold!")
        mean_val = hist.GetMean()
        rms_val = hist.GetRMS()
        entries = hist.GetEntries()
        stats_text = ROOT.TLatex()
        stats_text.SetNDC()
        stats_text.SetTextSize(0.03)
        stats_text.DrawLatex(0.65, 0.85, f"Entries: {int(entries)}")
        stats_text.DrawLatex(0.65, 0.82, f"Mean: {mean_val:.2f}")
        stats_text.DrawLatex(0.65, 0.79, f"RMS: {rms_val:.2f}")
        ROOT.gStyle.SetOptStat(0)
        can.SaveAs(f"plots/pion_shower_start_layer_{s}_ABSOLUTE_THRESHOLD.png")
        can.Close()

# Plot Profile Discrepancy for all pT slices on the same plot
if profile_discrepancy_hists:
    plotDiscrepancyHistograms(profile_discrepancy_hists, "plots/profile_discrepancy_all_pt_slices_MAIA_COMBINED.png",
                             "Max Profile Discrepancy", "Entries", threshold_line=0.6, with_bib=False)

# Also create individual plots for each pT slice (titles fixed to pionGun)
for s, hist in profile_discrepancy_hists.items():
    c = ROOT.TCanvas("can", "can", 800, 600)
    c.SetLogy()
    hist.SetLineColor(ROOT.kBlue)
    hist.SetLineWidth(3)
    hist.GetXaxis().SetTitle("Max Profile Discrepancy")
    hist.GetYaxis().SetTitle("Entries")
    hist.SetTitle(f"Profile Discrepancy - {s.replace('pionGun_pT_', '').replace('_', '-')} GeV")
    hist.Draw("HIST")
    y_max = hist.GetMaximum()
    line = ROOT.TLine(0.6, 1, 0.6, y_max)
    line.SetLineColor(ROOT.kRed)
    line.SetLineStyle(2)
    line.SetLineWidth(2)
    line.Draw()
    text = ROOT.TText(0.62, y_max * 0.7, "Pandora threshold = 0.6")
    text.SetTextSize(0.03)
    text.SetTextColor(ROOT.kRed)
    text.Draw()
    mean_val = hist.GetMean()
    rms_val = hist.GetRMS()
    entries = hist.GetEntries()
    stats_text = ROOT.TText(0.15, y_max * 0.8, f"Entries: {int(entries)}")
    stats_text.SetTextSize(0.03)
    stats_text.Draw()
    stats_text2 = ROOT.TText(0.15, y_max * 0.75, f"Mean: {mean_val:.3f}")
    stats_text2.SetTextSize(0.03)
    stats_text2.Draw()
    stats_text3 = ROOT.TText(0.15, y_max * 0.7, f"RMS: {rms_val:.3f}")
    stats_text3.SetTextSize(0.03)
    stats_text3.Draw()
    c.SaveAs(f"plots/profile_discrepancy_{s}.png")
    c.Close()

# Compute simple summaries (kept structure; keys updated to pionGun where visible)
for s, hist in profile_discrepancy_hists.items():
    pt_range = s.replace('pionGun_pT_', '').replace('_', '-')
    entries = int(hist.GetEntries())
    mean_val = hist.GetMean()
    rms_val = hist.GetRMS()
    total_entries = hist.GetEntries()
    above_threshold = 0
    for i in range(1, hist.GetNbinsX() + 1):
        bin_center = hist.GetBinCenter(i)
        if bin_center > 0.6:
            above_threshold += hist.GetBinContent(i)
    pct_above_threshold = (above_threshold / total_entries * 100) if total_entries > 0 else 0

# PFO shower start plots (now filled above)
pfo_shower_hists = {}
for s in hists:
    if "pfo_shower_start_layer" in hists[s] and hists[s]["pfo_shower_start_layer"].GetEntries() > 0:
        pfo_shower_hists[s] = hists[s]["pfo_shower_start_layer"]

if pfo_shower_hists:
    plotDiscrepancyHistograms(pfo_shower_hists, "plots/pfo_shower_start_layer_all_pt_slices_ABSOLUTE_THRESHOLD.png", "PFO Shower Start Layer (Absolute Threshold)", "Entries")
    for s, hist in pfo_shower_hists.items():
        pt_range = s.replace('pionGun_pT_', '').replace('_', '-')
        can = ROOT.TCanvas("can", "can", 800, 600)
        can.SetLeftMargin(0.12)
        can.SetBottomMargin(0.12)
        can.SetTopMargin(0.08)
        can.SetRightMargin(0.05)
        can.SetLogy()
        hist.SetLineColor(ROOT.kBlue)
        hist.SetLineWidth(3)
        hist.SetTitle("")
        hist.GetXaxis().SetTitle("PFO Shower Start Layer")
        hist.GetYaxis().SetTitle("Entries")
        hist.GetXaxis().SetTitleSize(0.045)
        hist.GetYaxis().SetTitleSize(0.045)
        hist.GetXaxis().SetLabelSize(0.04)
        hist.GetYaxis().SetLabelSize(0.04)
        hist.Draw("HIST")
        text = ROOT.TLatex()
        text.SetNDC()
        text.SetTextAlign(11)
        text.SetTextSize(0.04)
        text.SetTextFont(62)
        concept_text = ROOT.TLatex()
        concept_text.SetNDC()
        concept_text.SetTextAlign(31)
        concept_text.SetTextSize(0.035)
        concept_text.SetTextFont(42)
        concept_text.DrawLatex(0.94, 0.92, "MAIA Detector Concept")
        text.DrawLatex(0.15, 0.92, "Muon Collider")
        text.SetTextFont(42)
        text.DrawLatex(0.15, 0.88, "Simulation, no BIB")
        text.DrawLatex(0.15, 0.84, "|#eta| < 2.4")
        text.DrawLatex(0.15, 0.80, f"pT: {pt_range} GeV")
        text.DrawLatex(0.15, 0.76, "FIXED: Absolute Energy Threshold!")
        mean_val = hist.GetMean()
        rms_val = hist.GetRMS()
        entries = hist.GetEntries()
        stats_text = ROOT.TLatex()
        stats_text.SetNDC()
        stats_text.SetTextSize(0.03)
        stats_text.DrawLatex(0.65, 0.85, f"Entries: {int(entries)}")
        stats_text.DrawLatex(0.65, 0.82, f"Mean: {mean_val:.2f}")
        stats_text.DrawLatex(0.65, 0.79, f"RMS: {rms_val:.2f}")
        ROOT.gStyle.SetOptStat(0)
        can.SaveAs(f"plots/pfo_shower_start_layer_{s}_ABSOLUTE_THRESHOLD.png")
        can.Close()

# 2D hist section left intact; it only runs if something booked/fed
if 'hists2d' in locals():
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
