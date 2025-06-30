#!/usr/bin/env python3

import ROOT
import pyLCIO
import glob
import math

max_events = 50  # Look at more events to find clusters

ecal_collections = [
    "ECalBarrelCollection", "ECalEndcapCollection"
]

collections_to_read = ["PandoraClusters", "MCParticle"] + ecal_collections

def debug_cluster_analysis(cluster, event_num, cluster_num):
    """Debug function to understand what's in the clusters"""
    print(f"\n=== EVENT {event_num}, CLUSTER {cluster_num} ===")
    
    try:
        cluster_energy = cluster.getEnergy()
        print(f"Cluster Energy: {cluster_energy}")
    except Exception as e:
        print(f"Error getting cluster energy: {e}")
        return
    
    try:
        cluster_hits = cluster.getCalorimeterHits()
        print(f"Number of calorimeter hits: {len(cluster_hits)}")
        
        if len(cluster_hits) == 0:
            print("  -> NO CALORIMETER HITS IN CLUSTER!")
            # Try to see what other info we can get
            try:
                cluster_pos = cluster.getPosition()
                print(f"  But cluster has position: ({cluster_pos[0]}, {cluster_pos[1]}, {cluster_pos[2]})")
            except:
                pass
                
            try:
                # Check if cluster has subclusters or other associations
                print(f"  Cluster type: {type(cluster)}")
                # Try to access other cluster methods
                try:
                    shape = cluster.getShape()
                    print(f"  Cluster shape parameters: {len(shape)} values")
                except:
                    print("  No shape parameters")
                    
                try:
                    subdet_energies = cluster.getSubdetectorEnergies()
                    print(f"  Subdetector energies: {len(subdet_energies)} values")
                    for i, e in enumerate(subdet_energies):
                        print(f"    Subdet {i}: {e}")
                except:
                    print("  No subdetector energies")
                    
            except Exception as e2:
                print(f"  Error getting additional cluster info: {e2}")
            return
        
        # If we have hits, examine them
        for i, hit in enumerate(cluster_hits[:5]):  # First 5 hits
            try:
                hit_energy = hit.getEnergy()
                hit_pos = hit.getPosition()
                print(f"  Hit {i}: Energy={hit_energy}, Position=({hit_pos[0]}, {hit_pos[1]}, {hit_pos[2]})")
            except Exception as e:
                print(f"  Hit {i}: Error accessing hit data: {e}")
                
    except Exception as e:
        print(f"Error getting calorimeter hits: {e}")
        return

def debug_collections(event, event_num):
    """Debug what collections are available and their contents"""
    print(f"\n=== EVENT {event_num} COLLECTIONS ===")
    
    # Check all available collections
    collection_names = event.getCollectionNames()
    print(f"Available collections: {len(collection_names)}")
    for name in collection_names:
        try:
            collection = event.getCollection(name)
            print(f"  {name}: {collection.getNumberOfElements()} elements")
        except:
            print(f"  {name}: Error accessing")
    
    # Focus on ECAL collections
    print(f"\n=== ECAL COLLECTIONS ===")
    all_hits = 0
    for coll_name in ecal_collections:
        try:
            collection = event.getCollection(coll_name)
            num_hits = collection.getNumberOfElements()
            all_hits += num_hits
            print(f"  {coll_name}: {num_hits} hits")
            
            # Check first hit in each collection
            if num_hits > 0:
                hit = collection.getElementAt(0)
                try:
                    hit_energy = hit.getEnergy()
                    hit_pos = hit.getPosition()
                    print(f"    First hit: E={hit_energy}, pos=({hit_pos[0]}, {hit_pos[1]}, {hit_pos[2]})")
                except Exception as e:
                    print(f"    First hit: Error accessing data: {e}")
                    
        except Exception as e:
            print(f"  {coll_name}: Not found or error: {e}")
    
    print(f"Total ECAL hits across all collections: {all_hits}")

def main():
    """Debug main function"""
    # Just use one sample for debugging
    samples = glob.glob("/data/fmeloni/DataMuC_MuColl10_v0A/v0/reco/electronGun_pT_0_50")
    if not samples:
        print("No samples found!")
        return
    
    files = glob.glob(f"{samples[0]}/*.slcio")[:2]  # Just first 2 files
    print(f"Debugging with {len(files)} files from {samples[0]}")
    
    reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()
    reader.setReadCollectionNames(collections_to_read)
    
    event_count = 0
    
    for f in files:
        print(f"\n{'='*50}")
        print(f"PROCESSING FILE: {f}")
        print(f"{'='*50}")
        
        reader.open(f)
        
        for event in reader:
            if event_count >= max_events:
                break
                
            print(f"\n{'='*30}")
            print(f"EVENT {event_count}")
            print(f"{'='*30}")
            
            # Debug collections first
            debug_collections(event, event_count)
            
            # Now try clusters
            try:
                clusters = event.getCollection("PandoraClusters")
                num_clusters = clusters.getNumberOfElements()
                print(f"\n=== PANDORA CLUSTERS ===")
                print(f"Number of PandoraClusters: {num_clusters}")
                
                if num_clusters > 0:
                    for i in range(min(3, num_clusters)):  # Just first 3 clusters
                        cluster = clusters.getElementAt(i)
                        debug_cluster_analysis(cluster, event_count, i)
                else:
                    print("  -> No clusters to analyze")
                    
            except Exception as e:
                print(f"Error accessing PandoraClusters: {e}")
            
            event_count += 1
            
        reader.close()
        
        if event_count >= max_events:
            break
    
    print(f"\n{'='*50}")
    print("DEBUG COMPLETE")
    print(f"{'='*50}")

if __name__ == "__main__":
    main()
