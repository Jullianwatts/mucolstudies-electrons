#ifndef CaloHitComposition_h
#define CaloHitComposition_h 1

#include "marlin/Processor.h"
#include "lcio.h"

#include <string>
#include "TH1D.h"
#include <vector>

using namespace lcio;
using namespace marlin;

class CaloHitComposition : public Processor
{

public:
  virtual Processor *newProcessor() { return new CaloHitComposition; }

  CaloHitComposition();

  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual void init();

  /** Called for every run.
   */
  virtual void processRunHeader(LCRunHeader *run);

  /** Called for every event - the working horse.
   */
  virtual void processEvent(LCEvent *evt);

  virtual void check(LCEvent *evt);

  /** Called after data processing for clean up.
   */
  virtual void end();

  // Call to get collections
  void getCollection(LCCollection *&, const std::string &, LCEvent *);

protected:
  // Collection names for (in/out)put
  std::string m_inputHitCollection = "";
  std::string m_inputRelationCollection = "";

  int _nRun{};
  int _nEvt{};

  bool m_is_barrel = true;

  float m_min_z = -10.0;
  float m_max_z = 10.0;
  float m_min_r = 0.0;
  float m_max_r = 2.5;
  int m_nLayers = 50;

  // --- Output threshold histograms:
  std::vector<TH1D *> m_photon_histograms;
  std::vector<TH1D *> m_neutron_histograms;
  std::vector<TH1D *> m_other_histograms;
};

#endif
