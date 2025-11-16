#ifndef MCParticleFilter_h
#define MCParticleFilter_h 1

#include "marlin/Processor.h"

#include "lcio.h"
#include <map>
#include <vector>

#include <EVENT/LCCollection.h>

using namespace lcio;
using namespace marlin;

class MCParticleFilter : public Processor
{

public:
  virtual Processor *newProcessor() { return new MCParticleFilter; }

  MCParticleFilter();

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
  void getCollection(LCCollection *&, std::string, LCEvent *);

protected:
  // Collection names for (in/out)put
  std::string m_inputCollection = "";
  std::string m_outputCollection = "";

  // angle boundary (symmetric around pi/2)
  double m_angle_range = 0.2;

  int _nRun{};
  int _nEvt{};
};

#endif
