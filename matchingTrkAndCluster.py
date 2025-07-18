import numpy as np
import math
import glob
import ROOT
import pyLCIO

ROOT.gROOT.SetBatch()

max_events = 10
samples = glob.glob("/data/fmeloni/....")

files = {}
slices = ["0_50", "50_250", "250_1000", "1000_5000"]
for s in slices:
        files[f"electronGun_pT_{s}] = []
for s in samples:
        sname = s.split("/")[-1]
        if sname not in files: continue
        files[sname] = glob.glob(f"{s}/*.sclio)
print("Found files:")

#creating empty histogram objects
#variables = {}
#variables["E"] = {"nbins": 30, "xmin": 0, "xmax":)



#finding psuedo layer based off of MAIA geometry
def estimatePseudoLayer(hit):
    x, y, z = hit.getPosition()
    r = math.sqrt(x**2 + y**2)
    layer_thickness = 5.35 #mm

    if 185.7 <= r <= 212.5 and abs(z) <= 230.7:     #endcap region
        ecal_front_r = 1857.0 #mm
        depth = r *10.0 - ecal_front_r  # convert cm to mm, then subtract ecal front to find depth
    elif 31.0 <= r <= 212.5 and 230.7 <= abs(z) <= 257.5:
        ecal_front_z = 2307.0
        depth = abs(z) * 10.0 - ecal_front_z   # convert cm to mm

    pseudo_layer = int(depth / layer_thickness)
    return max(0, min(pseudo_layer, 49))  # between 0 and 49

def getInnerPseudoLayer(cluster):
    layers = [
            estimatePsuedoLayer(hit)
            for hit in cluster.getCalorimeterHits()
            if estimatePsudeoLayer(hit) != -1
    ]
    return min(layers) if layers else - 1



def getShowerProfile(hits, ecal_frontface_z=185.7):
    #Returns longitudinal profile and first pseudo-layer over threshold. need to figure out what .6 discrepancy looks like vs .7 or .8
    profile = np.zeros(50)  # One bin per pseudo-layer

    total_energy = 0.0
    for hit in hits:
        energy = hit.getEnergy()
        psuedo_layer = estimatePsuedoLayer(hit)
        if 0 <= pseudo_layer < 50:
            profile[pseudo_layer] += energy
            total_energy += energy
    if total_energy == 0:
        first_layer_over_threshold = -1
    else:
        threshold = THRESHOLD_FRACTION * total_energy
        first_layer_over_threshold = -1
        for i, energy in enumerate(profile):
            if energy > threshold:
                first_layer_over_threshold = i
                break
    return profile, first_layer_over_threshold
