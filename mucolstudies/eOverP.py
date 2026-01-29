import math
import pyLCIO
import ROOT
from pyLCIO import EVENT, IOIMPL

# --- Configuration ---
file_path = "/scratch/jwatts/mucol/data/reco/electronGun_pT_0_50/electronGun_pT_0_50_reco_4.slcio"
pdg_map = {11: "Electron", 13: "Muon", 22: "Photon", 211: "Pion", 2112: "Neutron"}

# Stats storage: stats[pdg] = [count, total_energy]
stats = {}
total_pfo_count = 0
event_count = 0

reader = IOIMPL.LCFactory.getInstance().createLCReader()
reader.open(file_path)

# Loop directly through the reader to count events accurately
for event in reader:
    event_count += 1
    
    # Access the full PFO collection
    try:
        pfo_coll = event.getCollection("PandoraPFOs")
    except:
        continue

    # Process EVERY PFO in the event
    for pfo in pfo_coll:
        p_type = abs(pfo.getType())
        
        # Calculate energy from all associated clusters
        e_calo = sum([c.getEnergy() for c in pfo.getClusters()])
        
        total_pfo_count += 1
        
        if p_type not in stats:
            stats[p_type] = [0, 0.0] # [count, total_energy]
        
        # Accumulate the count and energy globally
        stats[p_type][0] += 1
        stats[p_type][1] += e_calo

    # Stop exactly after 1000 events
    if event_count == 1000:
        break

reader.close()

# --- FINAL SUMMARY ---
print(f"({event_count} EVENTS)")
# Updated Header
print(f"{'Particle Type':<15} | {'Total Count':<12} | {'Avg/Event':<10} | {'Avg Energy/Part (GeV)':<18}")

for p_type in sorted(stats.keys()):
    count = stats[p_type][0]
    avg_per_event = count / event_count
    total_e = stats[p_type][1]
    
    # Calculate Average Energy per Particle
    avg_e_per_particle = total_e / count if count > 0 else 0
    
    name = pdg_map.get(p_type, f"ID({p_type})")
    print(f"{name:<15} | {count:<12} | {avg_per_event:<10.2f} | {avg_e_per_particle:<20.2f}")

print("-" * 65)
print(f"{'TOTAL PFOS':<15} | {total_pfo_count:<12} | {total_pfo_count/event_count:<10.2f} |")
