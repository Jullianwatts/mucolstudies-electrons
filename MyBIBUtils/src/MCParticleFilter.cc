#include "MCParticleFilter.h"
#include <iostream>

#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <EVENT/MCParticle.h>

#include "TMath.h"
#include "TVector3.h"

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"

using namespace lcio;
using namespace marlin;

MCParticleFilter aMCParticleFilter;

MCParticleFilter::MCParticleFilter() : Processor("MCParticleFilter")
{

    // Modify processor description
    _description = "MCParticleFilter applies space selections to reduce the BIB";

    // Input collection
    registerProcessorParameter("MCParticleCollectionName",
                               "Name of the MCParticle input collection",
                               m_inputCollection,
                               std::string("MCParticle"));

    // Output collection
    registerProcessorParameter("FilteredCollectionName",
                               "Selected MCParticle output collection",
                               m_outputCollection,
                               std::string("FilteredMCParticle"));

    // Angle range, symmetric around pi/2
    registerProcessorParameter("AngleRange",
                               "Cut in radians, around pi/2",
                               m_angle_range,
                               0.2);
}

void MCParticleFilter::init()
{

    streamlog_out(DEBUG) << "   init called  " << std::endl;

    // usually a good idea to
    printParameters();

    _nRun = 0;
    _nEvt = 0;
}

void MCParticleFilter::processRunHeader(LCRunHeader *run)
{

    _nRun++;
}

void MCParticleFilter::processEvent(LCEvent *evt)
{

    streamlog_out(DEBUG8) << "Processing event " << _nEvt << std::endl;

    // Get the collection of MCParticles
    LCCollection *MCParticleCollection = 0;
    getCollection(MCParticleCollection, m_inputCollection, evt);

    // Make the output collections
    LCCollectionVec *FilteredCollection = new LCCollectionVec(MCParticleCollection->getTypeName());
    FilteredCollection->setSubset(true);

    int nParts = MCParticleCollection->getNumberOfElements();

    // Loop over MCParticles    
    for (int itPart = 0; itPart < nParts; itPart++)
    {

        // Get the MCParticle
        MCParticle *part = static_cast<MCParticle *>(MCParticleCollection->getElementAt(itPart));

        // get the position
        TVector3 pos(part->getMomentum()[0], part->getMomentum()[1], part->getMomentum()[2]);

        if (fabs(pos.Phi() - TMath::Pi() / 2.) > m_angle_range)
        {
            streamlog_out(DEBUG0) << " -> rejected particle with dphi " << fabs(pos.Phi() - TMath::Pi() / 2.) << std::endl;
            continue;
        }

        FilteredCollection->addElement(part);
    }

    // Store the filtered hit collections
    evt->addCollection(FilteredCollection, m_outputCollection);

    //-- note: this will not be printed if compiled w/o MARLINDEBUG=1 !
    streamlog_out(DEBUG) << "   done processing event: " << evt->getEventNumber()
                         << "   in run:  " << evt->getRunNumber() << std::endl;
    
    _nEvt++;
}

void MCParticleFilter::check(LCEvent *evt)
{
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}

void MCParticleFilter::end()
{

    //   std::cout << "MCParticleFilter::end()  " << name()
    // 	    << " processed " << _nEvt << " events in " << _nRun << " runs "
    // 	    << std::endl ;
}

void MCParticleFilter::getCollection(LCCollection *&collection, std::string collectionName, LCEvent *evt)
{
    try
    {
        collection = evt->getCollection(collectionName);
    }
    catch (DataNotAvailableException &e)
    {
        streamlog_out(DEBUG5) << "- cannot get collection. Collection " << collectionName.c_str() << " is unavailable" << std::endl;
        return;
    }
    return;
}
