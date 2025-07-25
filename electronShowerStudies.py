import math
import glob
import ROOT
import pyLCIO

exec(open("./plotHelper.py").read())
ROOT.gROOT.SetBatch()
#ROOT.gROOT.Reset()

# Set up some options
max_events = -1

#samples = glob.glob("/data/fmeloni/DataMuC_MAIA_v0/v3/electronGun*")
samples = glob.glob("/data/fmeloni/DataMuC_MuColl10_v0A/v0/reco/electronGun*")
files = {}

slices = ["0_50", "50_250", "250_1000", "1000_5000"]
for s in slices:
    files[f"electronGun_pT_{s}"] = []

for s in samples:
    sname = s.split("/")[-1]
    if sname not in files: continue
    files[sname] = glob.glob(f"{s}/*.slcio")

# Set up histograms for spatial matching analysis
hists = {}
for s in files:
    hists[s] = {}
    # Basic cluster properties
    hists[s]["cluster_energy"] = ROOT.TH1F(s+"_cluster_energy", s, 100, 0, 5000)
    hists[s]["cluster_nhits"] = ROOT.TH1F(s+"_cluster_nhits", s, 50, 0, 100)
    hists[s]["cluster_r"] = ROOT.TH1F(s+"_cluster_r", s, 50, 0, 3000)
    hists[s]["cluster_eta"] = ROOT.TH1F(s+"_cluster_eta", s, 100, -5, 5)  # Add eta histogram
    hists[s]["cluster_theta"] = ROOT.TH1F(s+"_cluster_theta", s, 100, 0, 180)  # Add theta histogram

    # Spatial matching shower analysis (region-independent)
    hists[s]["shower_start_depth"] = ROOT.TH1F(s+"_shower_start_depth", s, 10, 0, 50)  # First hit depth, do 25, 0, 25 if 50,0,50 doesnt work
    hists[s]["shower_start_depth_barrel"] = ROOT.TH1F(s+"_shower_start_depth_barrel", s, 10, 0, 50)
    hists[s]["shower_start_depth_endcap"] = ROOT.TH1F(s+"_shower_start_depth_endcap", s, 10, 0, 50)
    hists[s]["shower_penetration"] = ROOT.TH1F(s+"_shower_penetration", s, 20, 0, 300)  # Last hit depth
    hists[s]["shower_penetration_barrel"] = ROOT.TH1F(s+"_shower_penetration_barrel", s, 20, 0, 300)
    hists[s]["shower_penetration_endcap"] = ROOT.TH1F(s+"_shower_penetration_endcap", s, 20, 0, 300)
    hists[s]["shower_max_position"] = ROOT.TH1F(s+"_shower_max_position", s, 64, 0, 160)  # Energy-weighted center
    hists[s]["shower_max_position_barrel"] = ROOT.TH1F(s+"_shower_max_position_barrel", s, 64, 0, 160)
    hists[s]["shower_max_position_endcap"] = ROOT.TH1F(s+"_shower_max_position_endcap", s, 64, 0, 160)
    hists[s]["shower_width"] = ROOT.TH1F(s+"_shower_width", s, 25, 0, 50)  # RMS width in depth
    hists[s]["shower_width_barrel"] = ROOT.TH1F(s+"_shower_width_barrel", s, 25, 0, 50)
    hists[s]["shower_width_endcap"] = ROOT.TH1F(s+"_shower_width_endcap", s, 25, 0, 50)
    hists[s]["matched_hits_count"] = ROOT.TH1F(s+"_matched_hits_count", s, 20, 0, 100)  # Number of matched hits
    hists[s]["matched_hits_count_barrel"] = ROOT.TH1F(s+"_matched_hits_count_barrel", s, 20, 0, 100)
    hists[s]["matched_hits_count_endcap"] = ROOT.TH1F(s+"_matched_hits_count_endcap", s, 20, 0, 100)
    hists[s]["cluster_hit_ratio"] = ROOT.TH1F(s+"_cluster_hit_ratio", s, 44, 1.0, 12.0)  # Matched hits / cluster energy
    hists[s]["cluster_hit_ratio_barrel"] = ROOT.TH1F(s+"_cluster_hit_ratio_barrel", s, 44, 1.0, 12.0)
    hists[s]["cluster_hit_ratio_endcap"] = ROOT.TH1F(s+"_cluster_hit_ratio_endcap", s, 44, 1.0, 12.0)

def calculate_eta_from_theta(theta_radians):
    #Calculate pseudorapidity eta from theta: eta = -ln(tan(theta/2))
    if theta_radians <= 0 or theta_radians >= math.pi:
        return float('inf')  # Handle edge cases

    tan_half_theta = math.tan(theta_radians / 2.0)
    if tan_half_theta <= 0:
        return float('inf')  # Handle edge cases

    eta = -math.log(tan_half_theta)
    #print(eta)
    return eta
def is_barrel_cluster(cluster):
    try:
        pos = cluster.getPosition()
        z = abs(pos[2])  # mm
        r = math.sqrt(pos[0]**2 + pos[1]**2)  # mm
        # ECAL Barrel: R = 185.7-212.5 cm, |Z| < 230.7 cm
        # ECAL Endcap: R = 31.0-212.5 cm, |Z| = 230.7-257.5 cm
        return z < 2307 and r > 1857  # Barrel region (in mm)
    except:
        return True

def get_available_hits(event, event_num):
    all_hits = []
    hit_sources = {}

    # Try different hit collections
    collections_to_try = [
        ('EcalBarrelCollectionRec', 'barrel_rec'),
        ('EcalEndcapCollectionRec', 'endcap_rec')
    ]

    for collection_name, source_name in collections_to_try:
        try:
            hit_collection = event.getCollection(collection_name)
            if len(hit_collection) > 0:
                all_hits.extend(hit_collection)
                hit_sources[source_name] = len(hit_collection)
                if event_num < 3:
                    print(f"        Found {len(hit_collection)} hits in {collection_name}")
        except:
            hit_sources[source_name] = 0
            continue

    if event_num < 3:
        print(f"        Total available hits: {len(all_hits)}")

    return all_hits, hit_sources

def spatial_match_hits(cluster, all_hits, match_radius_mm=80.0, event_num=0):
    try:
        cluster_pos = cluster.getPosition()
        cluster_x, cluster_y, cluster_z = cluster_pos[0], cluster_pos[1], cluster_pos[2]

        matched_hits = []

        for hit in all_hits:
            if not hit:  # Skip null hits
                continue

            try:
                hit_pos = hit.getPosition()
                hit_energy = hit.getEnergy()

                if hit_energy <= 0:  # Skip zero energy hits
                    continue

                # Calculate 3D distance between cluster center and hit
                dx = hit_pos[0] - cluster_x
                dy = hit_pos[1] - cluster_y
                dz = hit_pos[2] - cluster_z
                distance_3d = math.sqrt(dx*dx + dy*dy + dz*dz)

                if distance_3d <= match_radius_mm:
                    matched_hits.append(hit)

            except:
                continue

        if event_num < 3:
            print(f"        Matched {len(matched_hits)}/{len(all_hits)} hits within {match_radius_mm}mm")

        return matched_hits

    except Exception as e:
        if event_num < 3:
            print(f"        ERROR in spatial matching: {e}")
        return []

def analyze_matched_hits(cluster, matched_hits, event_num):
    try:
        if len(matched_hits) == 0:
            return None, None, None, None, None, None

        cluster_pos = cluster.getPosition()
        is_barrel = is_barrel_cluster(cluster)

        depths = []
        energies = []

        for hit in matched_hits:
            try:
                hit_pos = hit.getPosition()
                hit_energy = hit.getEnergy()

                # Calculate depth from ECAL surface
                if is_barrel:
                    r = math.sqrt(hit_pos[0]**2 + hit_pos[1]**2)
                    depth = r - 1857  # ECAL barrel starts at 185.7 cm = 1857 mm
                else:
                    z = abs(hit_pos[2])
                    depth = z - 2307  # ECAL endcap starts at 230.7 cm = 2307 mm

                if depth >= 0:  # Only consider hits inside ECAL
                    depths.append(depth)
                    energies.append(hit_energy)

            except:
                continue

        if len(depths) == 0:
            return None, None, None, None, None, None

        # Calculate shower properties
        total_energy = sum(energies)

        # 1. Energy-weighted shower center (where most energy is deposited)
        if total_energy > 0:
            weighted_depth = sum(d * e for d, e in zip(depths, energies))
            shower_max_depth = weighted_depth / total_energy
        else:
            shower_max_depth = sum(depths) / len(depths)

        # 2. Shower start (minimum depth)
        shower_start = min(depths)

        # 3. Shower penetration (maximum depth)
        shower_penetration = max(depths)

        # 4. RMS width of shower in depth direction
        mean_depth = sum(depths) / len(depths)
        if len(depths) > 1:
            variance = sum((d - mean_depth)**2 for d in depths) / len(depths)
            shower_width = math.sqrt(variance)
        else:
            shower_width = 0.0

        # 5. Hit efficiency (matched hits per GeV)
        cluster_energy = cluster.getEnergy()
        hit_ratio = len(matched_hits) / cluster_energy if cluster_energy > 0 else 0

        return shower_max_depth, shower_start, shower_penetration, shower_width, len(matched_hits), hit_ratio

    except Exception as e:
        if event_num < 3:
            print(f"        ERROR in analyze_matched_hits: {e}")
        return None, None, None, None, None, None

# Create a reader object
reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()
reader.setReadCollectionNames(["PandoraClusters", "MCParticle", "ECalBarrelCollection", "ECalEndcapCollection", "EcalBarrelCollectionRec", "EcalEndcapCollectionRec","EcalBarrelCollectionDigi","EcalEndcapCollectionDigi", "EcalEndcapCollectionConed", "EcalEndcapCollectionSel"])

# Loop over the different samples
for s in files:
    print("Working on sample", s)
    i = 0

    # Loop over the files in a sample
    for f in files[s]:
        reader.open(f)

        # Loop over events in each file
        for event in reader:
            if max_events > 0 and i >= max_events:
                break
            if i % 100 == 0:
                print(f"\tProcessing event: {i}")

            try:
                clusters = event.getCollection("PandoraClusters")

                # Get all available reconstructed hits for this event
                all_hits, hit_sources = get_available_hits(event, i)

                # Debug: count clusters
                n_clusters_total = len(clusters)
                n_clusters_processed = 0
                n_clusters_with_matches = 0

                if i < 3:
                    print(f"  Event {i}: {n_clusters_total} clusters, {len(all_hits)} total hits available")

                # Loop over clusters
                for cluster in clusters:
                    cluster_position = cluster.getPosition()
                    cluster_E = cluster.getEnergy()
                    cluster_x, cluster_y, cluster_z = cluster_position[0], cluster_position[1], cluster_position[2]
                    cluster_r = math.sqrt(cluster_x**2 + cluster_y**2)

                    # Get theta from LCIO and calculate eta using physics formula
                    cluster_theta = cluster.getITheta()  # theta in radians
                    cluster_eta = calculate_eta_from_theta(cluster_theta)
                    cluster_theta_deg = cluster_theta * 180.0 / math.pi  # For debugging

                    # Apply eta cut to focus on central region (comment out to see full distribution)
                    if abs(cluster_eta) > 1.5: continue

                    is_barrel = is_barrel_cluster(cluster)

                    # Debug: print cluster info for first few events
                    if i < 3:
                        region = "BARREL" if is_barrel else "ENDCAP"
                        print(f"    Cluster: E={cluster_E:.1f}GeV, R={cluster_r:.1f}mm, Z={abs(cluster_z):.1f}mm, theta={cluster_theta_deg:.1f}°, eta={cluster_eta:.2f} -> {region}")

                    n_clusters_processed += 1

                    # Fill basic histograms
                    hists[s]["cluster_energy"].Fill(cluster_E)
                    hists[s]["cluster_r"].Fill(cluster_r)
                    hists[s]["cluster_eta"].Fill(cluster_eta)
                    hists[s]["cluster_theta"].Fill(cluster_theta_deg)

                    # Spatial matching and shower analysis
                    matched_hits = spatial_match_hits(cluster, all_hits, match_radius_mm=80.0, event_num=i)

                    if len(matched_hits) > 0:
                        n_clusters_with_matches += 1

                        # Analyze shower properties from matched hits
                        shower_max, shower_start, shower_penetration, shower_width, n_matched, hit_ratio = analyze_matched_hits(
                            cluster, matched_hits, i
                        )

                        # Fill histograms if analysis succeeded
                        if shower_max is not None:
                            hists[s]["shower_max_position"].Fill(shower_max)
                            if is_barrel:
                                hists[s]["shower_max_position_barrel"].Fill(shower_max)
                            else:
                                hists[s]["shower_max_position_endcap"].Fill(shower_max)
                        if shower_start is not None:
                            hists[s]["shower_start_depth"].Fill(shower_start)
                            if is_barrel:
                                hists[s]["shower_start_depth_barrel"].Fill(shower_start)
                            else:
                                hists[s]["shower_start_depth_endcap"].Fill(shower_start)
                        if shower_penetration is not None:
                            hists[s]["shower_penetration"].Fill(shower_penetration)
                            if is_barrel:
                                hists[s]["shower_penetration_barrel"].Fill(shower_penetration)
                            else:
                                hists[s]["shower_penetration_endcap"].Fill(shower_penetration)
                        if shower_width is not None:
                            hists[s]["shower_width"].Fill(shower_width)
                            if is_barrel:
                                hists[s]["shower_width_barrel"].Fill(shower_width)
                            else:
                                hists[s]["shower_width_endcap"].Fill(shower_width)
                        if n_matched is not None:
                            hists[s]["matched_hits_count"].Fill(n_matched)
                            if is_barrel:
                                hists[s]["matched_hits_count_barrel"].Fill(n_matched)
                            else:
                                hists[s]["matched_hits_count_endcap"].Fill(n_matched)
                        if hit_ratio is not None:
                            hists[s]["cluster_hit_ratio"].Fill(hit_ratio)
                            if is_barrel:
                                hists[s]["cluster_hit_ratio_barrel"].Fill(hit_ratio)
                            else:
                                hists[s]["cluster_hit_ratio_endcap"].Fill(hit_ratio)

                # Debug: print counts for first few events
                if i < 3:
                    print(f"  -> Processed {n_clusters_processed}/{n_clusters_total} clusters, {n_clusters_with_matches} with matched hits")

            except Exception as e:
                # Print errors for first few events to help debug
                if i < 3:
                    print(f"    ERROR in event processing: {e}")
                pass
            i += 1

        reader.close()

# Print final statistics for debugging
print("\n=== FINAL STATISTICS ===")
for s in files:
    energy_entries = hists[s]["cluster_energy"].GetEntries()
    eta_entries = hists[s]["cluster_eta"].GetEntries()
    shower_max_entries = hists[s]["shower_max_position"].GetEntries()
    shower_start_entries = hists[s]["shower_start_depth"].GetEntries()
    matched_hits_entries = hists[s]["matched_hits_count"].GetEntries()
    print(f"{s}:")
    print(f"  Cluster energy entries: {energy_entries}")
    print(f"  Cluster eta entries: {eta_entries}")
    print(f"  Shower max entries: {shower_max_entries}")
    print(f"  Shower start entries: {shower_start_entries}")
    print(f"  Matched hits entries: {matched_hits_entries}")

# Helper for display names
label_map = {
    "electronGun_pT_0_50": "pT 0-50 GeV",
    "electronGun_pT_50_250": "pT 50-250 GeV",
    "electronGun_pT_250_1000": "pT 250-1000 GeV",
    "electronGun_pT_1000_5000": "pT 1000-5000 GeV"
}

# Plot histograms
hist_names_to_plot = [
    "cluster_energy",
    "cluster_eta",
    "cluster_theta",
    "shower_max_position",
    "shower_max_position_barrel",
    "shower_max_position_endcap",
    "shower_start_depth",
    "shower_start_depth_barrel",
    "shower_start_depth_endcap",
    "shower_penetration",
    "shower_penetration_barrel",
    "shower_penetration_endcap",
    "shower_width",
    "shower_width_barrel",
    "shower_width_endcap",
    "matched_hits_count",
    "matched_hits_count_barrel",
    "matched_hits_count_endcap",
    "cluster_hit_ratio",
    "cluster_hit_ratio_barrel",
    "cluster_hit_ratio_endcap"
]

for hist_name in hist_names_to_plot:
    hists_to_plot = {}
    for j, s in enumerate(hists):
        if hist_name in hists[s] and hists[s][hist_name].GetEntries() > 0:
            hists_to_plot[label_map.get(s, s)] = hists[s][hist_name]

    # Determine xlabel
    xlabel_map = {
        "cluster_energy": "Cluster Energy [GeV]",
        "cluster_eta": "Cluster Pseudorapidity η",
        "cluster_theta": "Cluster Polar Angle θ [degrees]",
        "shower_max_position": "Energy-Weighted Shower Center [mm]",
        "shower_max_position_barrel": "Energy-Weighted Shower Center [mm] (Barrel)",
        "shower_max_position_endcap": "Energy-Weighted Shower Center [mm] (Endcap)",
        "shower_start_depth": "Shower Start Depth [mm]",
        "shower_start_depth_barrel": "Shower Start Depth [mm] (Barrel)",
        "shower_start_depth_endcap": "Shower Start Depth [mm] (Endcap)",
        "shower_penetration": "Maximum Shower Penetration [mm]",
        "shower_penetration_barrel": "Maximum Shower Penetration [mm] (Barrel)",
        "shower_penetration_endcap": "Maximum Shower Penetration [mm] (Endcap)",
        "shower_width": "Shower RMS Width [mm]",
        "shower_width_barrel": "Shower RMS Width [mm] (Barrel)",
        "shower_width_endcap": "Shower RMS Width [mm] (Endcap)",
        "matched_hits_count": "Number of Matched Hits per Cluster",
        "matched_hits_count_barrel": "Number of Matched Hits per Cluster (Barrel)",
        "matched_hits_count_endcap": "Number of Matched Hits per Cluster (Endcap)",
        "cluster_hit_ratio": "Matched Hits per GeV",
        "cluster_hit_ratio_barrel": "Matched Hits per GeV (Barrel)",
        "cluster_hit_ratio_endcap": "Matched Hits per GeV (Endcap)"
    }
    xlabel = xlabel_map.get(hist_name, hist_name)

    # Call plotting function
    if hists_to_plot:
        plotHistograms(hists_to_plot, "plots/"+hist_name+".png", xlabel, "Entries")
        print(f"Created plot: plots/{hist_name}.png")

print("Analysis complete!")
