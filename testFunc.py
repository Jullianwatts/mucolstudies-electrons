from pyLCIO import IOIMPL

reader = IOIMPL.LCFactory.getInstance().createLCReader()
reader.open("/data/fmeloni/DataMuC_MAIA_v0/v5/reco/electronGun_pT_0_50/electronGun_pT_0_50_reco_5480.slcio")
event = reader.readNextEvent()

clusters = event.getCollection("PandoraClusters")
cluster = clusters[0]

print("Cluster object type:", type(cluster))
print("Available methods:", dir(cluster))

