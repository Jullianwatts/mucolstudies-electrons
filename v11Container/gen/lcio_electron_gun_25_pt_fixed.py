#####################################
#
# simple script to create lcio files with single particle
# events - modify as needed
# @author F.Gaede, DESY
# @date 1/07/2014
#
# initialize environment:
# export PYTHONPATH=${LCIO}/src/python:${ROOTSYS}/lib
#
#####################################
import math
import random
from array import array

# --- LCIO dependencies ---
from pyLCIO import UTIL, EVENT, IMPL, IO, IOIMPL

#---- number of events ----------------------
nevt = 10000

outfile = "electronGun_gen25_fixed.slcio"

#--------------------------------------------

wrt = IOIMPL.LCFactory.getInstance().createLCWriter( )

wrt.open( outfile , EVENT.LCIO.WRITE_NEW )

random.seed()

#========== particle properties ===================

# particles per event
npart = 1

genstat = 1

# Fixed transverse momentum at 25 GeV
fixed_pt = 25.0

# Fixed theta at pi/2
fixed_theta = math.pi / 2.

pdg = 11

mass = 0.0005109988851472735

charge = -1.

decayLen = 1.e32

beamspot_sigma = 1.5 #mm

#=================================================

for j in range( 0, nevt ):

    col = IMPL.LCCollectionVec( EVENT.LCIO.MCPARTICLE )

    evt = IMPL.LCEventImpl()

    evt.setEventNumber( j )

    evt.addCollection( col , "MCParticle" )

    print (j, "-----------------------------")

    for ipart in range( 0, npart ):

        pt = fixed_pt

        theta = fixed_theta

        phi = random.random() * math.pi * 2.

        p = pt/math.sin( theta )

        energy = math.sqrt( mass*mass + p * p )

        if energy > 5000:
            continue

        px = pt * math.cos( phi )
        py = pt * math.sin( phi )
        pz = p * math.cos( theta )

        momentum = array('f',[ px, py, pz ] )

        # --- endpoint

        epx = decayLen * math.cos( phi ) * math.sin( theta )
        epy = decayLen * math.sin( phi ) * math.sin( theta )
        epz = decayLen * math.cos( theta )

        endpoint = array('d',[ epx, epy, epz ] )

        # --- production vertex

        vpx = 0.
        vpy = 0.
        vpz = random.gauss(0., beamspot_sigma)

        vertex = array('d',[ vpx, vpy, vpz ] )

        time = 0.

        # --- particle charge

        if ipart % 2 == 1:
            pdg = -pdg
            charge = -charge

        #--------------- create MCParticle -------------------

        mcp = IMPL.MCParticleImpl()

        mcp.setGeneratorStatus( genstat )
        mcp.setMass( mass )
        mcp.setPDG( pdg )
        mcp.setMomentum( momentum )
        mcp.setCharge( charge )
        mcp.setVertex( vertex )
        mcp.setTime( time )

        if( decayLen < 1.e9 ) : # arbitrary ...
            mcp.setEndpoint( endpoint )

        print (" ", ipart, pdg, charge, pt, phi, theta)

        #-------------------------------------------------------

        col.addElement( mcp )

    wrt.writeEvent( evt )

print("Generated fixed pt =", fixed_pt, "GeV, fixed theta = pi/2 sample")

wrt.close()
