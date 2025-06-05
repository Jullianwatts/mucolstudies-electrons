from pyLCIO import IOIMPL
import glob
import math

files = glob.glob("electronGun_pT_1000_5000*.slcio")

total_mcps = 0
total_el = 0

def calc_eta(px, py, pz):
    p = math.sqrt(px**2 + py**2 + pz**2)
    if p - abs(pz) == 0:
        return float("inf")  # Avoid division by zero
    return 0.5 * math.log((p + pz) / (p - pz))

for f in files:
    print(f"Processing {f}...")
    reader = IOIMPL.LCFactory.getInstance().createLCReader()
    reader.open(f)

    event = reader.readNextEvent(0)
    while event is not None:
        try:
            if not event.getCollectionNames().contains("MCParticle"):
                print("No MCParticle collection in this event")
                event = reader.readNextEvent(0)
                continue

            mcps = event.getCollection("MCParticle")

            for mcp in mcps:
                total_mcps += 1
                if not mcp.getGeneratorStatus() == 1:
                    continue
                px, py, pz = mcp.getMomentum()
                pdg = mcp.getPDG()
                eta = calc_eta(px, py, pz)
                if abs(pdg) == 11 and abs(eta) <= 2:
                    total_el += 1

        except Exception as e:
            print(f"Error: {e}")

        event = reader.readNextEvent(0)

    reader.close()

print(f"\nTotal MCParticles: {total_mcps}")
print(f"Filtered electrons (PDG=±11, |η|≤2, status=1): {total_el}")

