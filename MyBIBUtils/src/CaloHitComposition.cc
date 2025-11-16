#include "CaloHitComposition.h"
#include <iostream>
#include <vector>
#include <map>
#include <math.h>

#include <EVENT/LCCollection.h>
#include <EVENT/CalorimeterHit.h>
#include <EVENT/SimCalorimeterHit.h>
#include <EVENT/LCRelation.h>
#include <EVENT/MCParticle.h>

#include <UTIL/CellIDDecoder.h>

#include <marlin/AIDAProcessor.h>

#include "TVector3.h"
#include "TString.h"

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"

using namespace lcio;
using namespace marlin;

CaloHitComposition aCaloHitComposition;

CaloHitComposition::CaloHitComposition() : Processor("CaloHitComposition")
{

    // Modify processor description
    _description = "CaloHitComposition applies E selections to reduce the BIB";

    // Input collection
    registerProcessorParameter("CaloHitCollectionName",
                               "Name of the CalorimeterHit input collection",
                               m_inputHitCollection,
                               std::string("EcalBarrelCollectionRec"));

    // Input relation collection
    registerProcessorParameter("CaloHitRelationCollectionName",
                               "Name of the CalorimeterHit-SimCalorimeterHit relation input collection",
                               m_inputRelationCollection,
                               std::string("EcalBarrelRelationsSimRec"));

    // Is it a barrel?
    registerProcessorParameter("IsBarrel",
                               "Flag to set if the calo is barrel (true) or endcap (false)",
                               m_is_barrel,
                               true);

    // N calo layers
    registerProcessorParameter("Nlayers",
                               "Number of calorimeter layers",
                               m_nLayers,
                               50);

    // Z min
    registerProcessorParameter("Zmin",
                               "Minimum z coordinate of the calorimeter",
                               m_min_z,
                               -10.0f);
    // Z max
    registerProcessorParameter("Zmax",
                               "Maximum z coordinate of the calorimeter",
                               m_max_z,
                               10.0f);
    // R min
    registerProcessorParameter("Rmin",
                               "Minimum r coordinate of the calorimeter",
                               m_min_r,
                               0.0f);
    // R max
    registerProcessorParameter("Rmax",
                               "Maximum r coordinate of the calorimeter",
                               m_max_r,
                               2.5f);

}

void CaloHitComposition::init()
{

    streamlog_out(DEBUG) << "   init called  " << std::endl;

    // usually a good idea to
    printParameters();

    _nRun = 0;
    _nEvt = 0;

    // --- Initialize the AIDAProcessor and book the diagnostic histograms:

    AIDAProcessor::histogramFactory(this);

    for (int iLayer = 0; iLayer < m_nLayers; iLayer++){
        TString name_photons;
        name_photons.Form("%s%s%i", m_inputHitCollection.c_str(), "_E_photon_L", int(iLayer));
        TString name_neutrons;
        name_neutrons.Form("%s%s%i", m_inputHitCollection.c_str(), "_E_neutron_L", int(iLayer));
        TString name_others;
        name_others.Form("%s%s%i", m_inputHitCollection.c_str(), "_E_other_L", int(iLayer));

        TH1D *histo_photons = new TH1D(name_photons, name_photons, 100, 0, 0.1);
        TH1D *histo_neutrons = new TH1D(name_neutrons, name_neutrons, 100, 0, 0.1);
        TH1D *histo_others = new TH1D(name_others, name_others, 100, 0, 0.1);
        
        m_photon_histograms.push_back(histo_photons);
        m_neutron_histograms.push_back(histo_neutrons);
        m_other_histograms.push_back(histo_others);
    }
}

void CaloHitComposition::processRunHeader(LCRunHeader *run)
{

    _nRun++;
}

void CaloHitComposition::processEvent(LCEvent *evt)
{

    streamlog_out(DEBUG) << "Processing event " << _nEvt << std::endl;

    // Get the collection of calo hits
    LCCollection *caloHitCollection = 0;
    getCollection(caloHitCollection, m_inputHitCollection, evt);

    LCCollection *inputHitRel = 0;
    getCollection(inputHitRel, m_inputRelationCollection, evt);

    std::string encoderString = caloHitCollection->getParameters().getStringVal("CellIDEncoding");
    UTIL::CellIDDecoder<CalorimeterHit> myCellIDEncoding(encoderString);

    // Loop over calo hits
    int nHits = caloHitCollection->getNumberOfElements();
    for (int itHit = 0; itHit < nHits; itHit++)
    {
        std::cout << "in hit loop " << itHit << std::endl;

        // Get the hit
        CalorimeterHit *hit = static_cast<CalorimeterHit *>(caloHitCollection->getElementAt(itHit));
        unsigned int layer = myCellIDEncoding(hit)["layer"];
        double hit_z = hit->getPosition()[2];

        // hit position
        TVector3 hitPos(hit->getPosition()[0], hit->getPosition()[1], hit_z);
        double hit_r = hitPos.Perp();
        
        if (m_is_barrel)
        {
            if (hit_z < m_min_z || hit_z > m_max_z)
                continue;
        }
        else
        {
            if (hit_r < m_min_r || hit_r > m_max_r)
                continue;
        }

        std::cout << "passed cuts " << itHit << std::endl;

        LCRelation *rel = static_cast<LCRelation *>(inputHitRel->getElementAt(itHit));
        SimCalorimeterHit *simhit = static_cast<SimCalorimeterHit *>(rel->getTo());

        size_t n_contributions =  simhit->getNMCContributions();

        float_t photon_frac = 0.;
        float_t neutron_frac = 0.;
        float_t other_frac = 0.;

        streamlog_out(DEBUG) << "n contributions: " << n_contributions << std::endl;

        for (size_t icont = 0; icont < n_contributions; icont++)
        {
            int pdg_shower = static_cast<MCParticle *>(simhit->getParticleCont(icont))->getPDG();

            if (pdg_shower == 2112)
            {
                neutron_frac += simhit->getEnergyCont(icont);
            }
            else if (pdg_shower == 22)
            {
                photon_frac += simhit->getEnergyCont(icont);
            }
            else
            {
                other_frac += simhit->getEnergyCont(icont);
            }

            float_t sumsimE = photon_frac + neutron_frac + other_frac;
            photon_frac = photon_frac / sumsimE;
            neutron_frac = neutron_frac / sumsimE;
            other_frac = other_frac / sumsimE;

            m_photon_histograms[layer]->Fill(hit->getEnergy(), photon_frac);
            m_neutron_histograms[layer]->Fill(hit->getEnergy(), neutron_frac);
            m_other_histograms[layer]->Fill(hit->getEnergy(), other_frac);

        }

    }

    //-- note: this will not be printed if compiled w/o MARLINDEBUG=1 !
    streamlog_out(DEBUG) << "   done processing event: " << evt->getEventNumber()
                         << "   in run:  " << evt->getRunNumber() << std::endl;

    _nEvt++;
}

void CaloHitComposition::check(LCEvent *evt)
{
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}

void CaloHitComposition::end()
{
    //   streamlog_out(DEBUG) << "CaloHitComposition::end()  " << name()
    // 	    << " processed " << _nEvt << " events in " << _nRun << " runs "
    // 	    << std::endl ;
}

void CaloHitComposition::getCollection(LCCollection *&collection, const std::string &collectionName, LCEvent *evt)
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
