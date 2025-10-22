import glob
from pyLCIO import IOIMPL, EVENT

dirs = glob.glob("/data/fmeloni/DataMuC_MAIA_v0/v5/reco/electronGun_pT_*")
files = []
for d in dirs:
    files.extend(glob.glob(f"{d}/*.slcio"))
print(f"Found {len(files)} .slcio files across {len(dirs)} directories")

reader = IOIMPL.LCFactory.getInstance().createLCReader()
reader.setReadCollectionNames(["PandoraClusters"])  # only load what you need

max_files = 5      # number of files to open
max_events = 2     # number of events per file

for f_index, fname in enumerate(files):
    if f_index >= max_files:
        break

    print(f"\nOpening {fname}")
    reader.open(fname)

    for evt_num, event in enumerate(reader):
        if evt_num >= max_events:
            print("broke at event number")
            break

        if "PandoraClusters" not in event.getCollectionNames():
            print("broke at collection names")
            continue

        clusters = event.getCollection("PandoraClusters")
        for cl in clusters:
            pids = cl.getParticleIDs()
            if not pids:
                print("get particle id fails")
                continue
            for pid in pids:
                print(f"Type: {pid.getType()}, "
                      f"Likelihood: {pid.getLikelihood():.3f}, "
                      f"Params: {list(pid.getParameters())}")

    reader.close()

