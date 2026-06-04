import math
import pyLCIO
from pyLCIO import EVENT, IOIMPL

file_path = "/scratch/jwatts/mucol/v2.11/reco/electronGun_pT_0_50/electronGun_pT_0_50_reco_0.slcio"
pdg_map = {
    11: "Electron",
    13: "Muon",
    22: "Photon",
    211: "Pion",
    2112: "Neutron"
}

stats = {}
total_pfo_count = 0
event_count = 0

reader = IOIMPL.LCFactory.getInstance().createLCReader()
reader.open(file_path)

for event in reader:

    event_count += 1

    # Print cluster energies for first 50 events
    if event_count <= 50:
        try:
            clusters = event.getCollection("PandoraClusters")
            energies = [c.getEnergy() for c in clusters]
            print(f"Event {event_count} Cluster Energies: {energies}")
        except:
            print(f"Event {event_count} Cluster Energies: []")

    try:
        pfo_coll = event.getCollection("PandoraPFOs")
    except:
        continue

    for pfo in pfo_coll:

        p_type = abs(pfo.getType())

        e_calo = sum(c.getEnergy() for c in pfo.getClusters())

        total_pfo_count += 1

        if p_type not in stats:
            stats[p_type] = [0, 0.0]

        stats[p_type][0] += 1
        stats[p_type][1] += e_calo

    if event_count == 50:
        break

reader.close()

print(f"\n({event_count} EVENTS)\n")
print(f"{'Particle Type':<15} | {'Total Count':<12} | {'Avg/Event':<10} | {'Avg Energy/Part (GeV)':<18}")

for p_type in sorted(stats.keys()):

    count = stats[p_type][0]
    avg_per_event = count / event_count
    total_e = stats[p_type][1]

    avg_e_per_particle = total_e / count if count > 0 else 0

    name = pdg_map.get(p_type, f"ID({p_type})")

    print(
        f"{name:<15} | "
        f"{count:<12} | "
        f"{avg_per_event:<10.2f} | "
        f"{avg_e_per_particle:<20.2f}"
    )

print(
    f"{'TOTAL PFOS':<15} | "
    f"{total_pfo_count:<12} | "
    f"{total_pfo_count/event_count:<10.2f} |"
)
