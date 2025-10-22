import glob
from pyLCIO import IOIMPL, EVENT

# ----------------------------------------------
# Collect all .slcio files across pT directories
# ----------------------------------------------
dirs = glob.glob("/data/fmeloni/DataMuC_MAIA_v0/v5/reco/electronGun_pT_*")
files = []
for d in dirs:
    files.extend(glob.glob(f"{d}/*.slcio"))
print(f"Found {len(files)} .slcio files across {len(dirs)} directories")

# ----------------------------------------------
# Set up LCIO reader
# ----------------------------------------------
reader = IOIMPL.LCFactory.getInstance().createLCReader()
reader.setReadCollectionNames(["PandoraClusters"])  # only load what you need

max_files = 5      # number of files to open
max_events = 2     # number of events per file

# ----------------------------------------------
# Loop through files and inspect PID data
# ----------------------------------------------
for f_index, fname in enumerate(files):
    if f_index >= max_files:
        break

    print(f"\nOpening {fname}")
    reader.open(fname)

    for evt_num, event in enumerate(reader):
        if evt_num >= max_events:
            print("Reached max_events limit")
            break

        if "PandoraClusters" not in event.getCollectionNames():
            print("PandoraClusters collection missing")
            continue

        clusters = event.getCollection("PandoraClusters")

        for i, cl in enumerate(clusters):
            try:
                pids = cl.getParticleIDs()
                print(f"Cluster {i}: PID vector length = {len(pids)}")

                if len(pids) == 0:
                    print("   No ParticleID info stored for this cluster")
                    continue

                for pid in pids:
                    print(f"   Type: {pid.getType()}, "
                          f"Likelihood: {pid.getLikelihood():.3f}, "
                          f"Params: {list(pid.getParameters())}")
            except Exception as e:
                print(f"   Cluster {i} failed to getParticleIDs(): {e}")
                continue

    reader.close()

