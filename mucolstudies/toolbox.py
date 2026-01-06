import glob
from pyLCIO import IOIMPL, EVENT, IMPL
dirs = glob.glob("/data/fmeloni/DataMuC_MAIA_v0/v5/reco/electronGun_pT_*")
files = []
for d in dirs:
    files.extend(glob.glob(f"{d}/*.slcio"))

reader = IOIMPL.LCFactory.getInstance().createLCReader()
reader.setReadCollectionNames(["PandoraClusters"])  # only what you need

max_files = 5      # number of files to open
max_events = 2     # number of events per file
for f_index, fname in enumerate(files):
    if f_index >= max_files:
        break

    print(f"\nOpening {fname}")
    reader.open(fname)

    for evt_num, event in enumerate(reader):
        if evt_num >= max_events:
            print("Reached max_events limit")
            break
        print(dir(EVENT))
        break
        clusters = event.getCollection("PandoraClusters")
        if len(clusters) == 0:
            print("No clusters found.")
            continue

        cluster = clusters[0]
        print(f"Testing cluster 0 from event {evt_num}")
        #trying to do toolbox.ParticleID
        pid_event = EVENT.ParticleID
        pid_impl = IMPL.ParticleIDImpl()

        print("Trying IsEmShower on both EVENT and IMPL ParticleID classes...")
        
        #now trying isit?

        try:
            print("EVENT result:", pid_event.IsEmShower(cluster))
        except Exception as e:
            print("EVENT::ParticleID IsEmShower failed:", e)

        try:
            print("IMPL result:", pid_impl.IsEmShower(cluster))
        except Exception as e:
            print("IMPL::ParticleIDImpl IsEmShower failed:", e)

        try:
            print("EVENT result:", pid_event.isEmShower(cluster))
        except Exception as e:
            print("EVENT::ParticleID isEmShower failed:", e)

        try:
            print("IMPL result:", pid_impl.isEmShower(cluster))
        except Exception as e:
            print("IMPL::ParticleIDImpl isEmShower failed:", e)

    reader.close()

