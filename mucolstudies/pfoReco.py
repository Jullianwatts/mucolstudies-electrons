import pyLCIO
from pyLCIO import EVENT

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

file_path = "/scratch/jwatts/mucol/data/reco/electronGun_pT_0_50/electronGun_pT_0_50_reco_5.slcio"

# Dictionary to store counts
counts = {}
event_count = 0
total_pfos = 0  # Added to track the grand total

reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()

try:
    reader.open(file_path)
    print(f"Reading file: {file_path}")

    for event in reader:
        event_count += 1
        if "PandoraPFOs" in event.getCollectionNames():
            pfos = event.getCollection("PandoraPFOs")
            
            # Increment total PFO count by the size of the collection
            total_pfos += pfos.getNumberOfElements()
            
            for pfo in pfos:
                pfo_type = abs(pfo.getType())
                counts[pfo_type] = counts.get(pfo_type, 0) + 1

    reader.close()
except Exception as e:
    print(f"Error opening file: {e}")

# FINAL SUMMARY PRINTING
print("-" * 30)
print(f"SUMMARY FOR {event_count} EVENTS")
print(f"Total PFOs found: {total_pfos}")
print("-" * 30)

if not counts:
    print("No PandoraPFOs found in any event.")
else:
    # Sort by count (highest first)
    sorted_types = sorted(counts.items(), key=lambda item: item[1], reverse=True)
    for pfo_type, count in sorted_types:
        name = pdg_map.get(pfo_type, f"Unknown ({pfo_type})")
        print(f"  - {name:<18}: Total = {count:<8}")
