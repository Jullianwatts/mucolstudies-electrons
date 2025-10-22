import glob
import pyLCIO
from pyLCIO import EVENT

max_events = 5

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

reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()
for sample, filelist in files.items():
    if not filelist:
        continue

    print(f"\n>>> Sample: {sample}")
    for fname in filelist:
        print(f"  Opening {fname}")
        reader.open(fname)

        for ievt, event in enumerate(reader):
            print(f"    Event {event.getEventNumber()}")

            if "PandoraPFOs" in event.getCollectionNames():
                pfos = event.getCollection("PandoraPFOs")
                for i, pfo in enumerate(pfos):
                    print(f"      PFO {i} -> type = {pfo.getType()}")
            else:
                print("    [No PandoraPFOs collection in this event]")

            if ievt + 1 >= max_events:
                break  # stop after max_events per file

        reader.close()

