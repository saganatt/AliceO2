// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"

#include "Analysis/EventSelection.h"
#include "Analysis/Centrality.h"
#include "Analysis/StepTHn.h"
#include "Analysis/CorrelationContainer.h"

#include <TH1F.h>
#include <cmath>
#include <TDirectory.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};

struct FemtoscopyTask {
  const double PionMass = 0.13956995;
  const double KaonMass = 0.493677;
  const double ProtonMass = 0.938272013;
  const double LambdaMass = 1.115683;

  const int numOfMultBins = 1;
  const int numOfChTypes = 16; //13
  const int numOfkTbins = 1;

  bool performSharedDaughterCut = false;
  bool enablePairMonitors = true;

  int runmults[numOfMultBins] = {1};
  int multbins[numOfMultBins + 1] = {0, 1000};

  int runch[numOfChTypes] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0}; // 1 - wlacza czastki do analizy
  const char* chrgs[numOfChTypes] = {"PP", "aPaP", "PaP", "KpKp", "KmKm", "KpKm", "PIpPIp", "PImPIm", "PIpPIm", "V0LL", "V0ALAL", "V0LAL", "all", "plus", "minus", "mixed"};

  int runqinv = 1;

  int runtype = 0; // Types 0 - global, 1 - ITS only, 2 - TPC Inner	//global tracks ->mfit ITS+TPC
  int owncuts = 0;
  int owndca = 0;

  int gammacut = 1; // cut na ee z gamma

  double shqmax = 2.0;
  int nbinssh = 20;

  // Filters and input definitions
  //#define MYFILTER
  //#ifdef MYFILTER
  //  Filter trackFilter = (aod::etaphi::etam > -0.8f) && (aod::etaphi::etam < 0.8f) && (aod::etaphi::ptm > 1.0f);
  //  using myTracks = soa::Filtered<aod::Tracks>;
  //#else
  using myTracks = aod::Tracks;
  //#endif

  // Configuration
  // TODO: Filterbit 96 == ((aod::track::isGlobalTrack == (uint8_t)1) || (aod::track::isGlobalTrackSDD == (uint8_t)1))
  O2_DEFINE_CONFIGURABLE(cfgFilterBit, int, 96, "Select filterbit");
  // TODO: Corresponding cut in O2?
  O2_DEFINE_CONFIGURABLE(cfgMinPlpContribSPD, int, 3, "SPD pile-up cut");
  O2_DEFINE_CONFIGURABLE(cfgMultBinNo, int, 200, "Number of multiplicity bins");
  O2_DEFINE_CONFIGURABLE(cfgZVertBinNo, int, 10, "Number of z-vertex bins");

  // TODO: Boolean once enabled
  O2_DEFINE_CONFIGURABLE(cfgIfGlobalTracks, int, 0, "Whether to use global tracks: 1 - yes, 0 - no");
  O2_DEFINE_CONFIGURABLE(cfgIfElectronRejection, int, 1, "Whether to reject electrons: 1 - yes, 0 - no");
  O2_DEFINE_CONFIGURABLE(cfgIfIsPileUp, int, 1, "Whether to pile up: 1 - yes, 0 - no");
  // TODO: Are there 'monitors' in O2? What is this?
  O2_DEFINE_CONFIGURABLE(cfgIfMonitors, int, 1, "Whether to use monitors: 1 - yes, 0 - no");
  O2_DEFINE_CONFIGURABLE(cfgIfV0Monitors, int, 0, "Whether to use V0 monitors: 1 - yes, 0 - no");

  O2_DEFINE_CONFIGURABLE(cfgNSigmaVal, float, 2.0f, "nSigmaVal");
  O2_DEFINE_CONFIGURABLE(cfgNSigmaVal2, float, 3.0f, "nSigmaVal2");
  O2_DEFINE_CONFIGURABLE(cfgNEtaMin, float, -0.8f, "nEtaMin");
  O2_DEFINE_CONFIGURABLE(cfgNEtaMax, float, 0.8f, "nEtaMax");
  O2_DEFINE_CONFIGURABLE(cfgMaxPt, float, 2.5f, "MaxPt");

  void init(o2::framework::InitContext&)
  {
  }

  // Use centrality percentile from V0M
  // See AliPhysics/PWGCF/FEMTOSCOPY/AliFemto/AliFemtoEventReaderAOD.cxx:595
  // Reader->SetUseMultiplicity(AliFemtoEventReaderAOD::kCentrality);

  // Version with explicit nested loop
  void process(soa::Join<aod::Collisions, aod::EvSels, aod::Cents>::iterator const& collision, myTracks const& tracks)
  {
    LOGF(info, "Tracks for collision: %d | Trigger mask: %lld | INT7: %d | V0M: %.1f", tracks.size(), collision.bc().triggerMask(), collision.sel7(), collision.centV0M());

    if (!collision.sel7())
      return;

    int bSign = 1; // TODO magnetic field from CCDB
    const float pTCut = 1.0;

    for (auto track1 = tracks.begin(); track1 != tracks.end(); ++track1) {

#ifndef MYFILTER
      if (track1.pt2() < pTCut)
        continue;
      if (track1.eta2() < -0.8 || track1.eta2() > 0.8)
        continue;
#endif

      if (cfgTriggerCharge != 0 && cfgTriggerCharge * track1.charge() < 0)
        continue;

      //LOGF(info, "TRACK %f %f | %f %f | %f %f", track1.eta(), track1.eta2(), track1.phi(), track1.phi2(), track1.pt(), track1.pt2());

      double eventValues[3];
      eventValues[0] = track1.pt();
      eventValues[1] = collision.centV0M();
      eventValues[2] = collision.posZ();

      same->getEventHist()->Fill(eventValues, CorrelationContainer::kCFStepReconstructed);
      //mixed->getEventHist()->Fill(eventValues, CorrelationContainer::kCFStepReconstructed);

      for (auto track2 = track1 + 1; track2 != tracks.end(); ++track2) {
#ifndef MYFILTER
        if (track2.pt2() < pTCut)
          continue;
        if (track2.eta2() < -0.8 || track2.eta2() > 0.8)
          continue;
#endif

        if (cfgAssociatedCharge != 0 && cfgAssociatedCharge * track2.charge() < 0)
          continue;
        if (cfgPairCharge != 0 && cfgPairCharge * track1.charge() * track2.charge() < 0)
          continue;

        if (cfg.mPairCuts && conversionCuts(track1, track2))
          continue;

        if (cfgTwoTrackCut > 0 && twoTrackCut(track1, track2, bSign))
          continue;

        double values[6] = {0};

        values[0] = track1.etam() - track2.etam();
        values[1] = track1.pt();
        values[2] = track2.pt();
        values[3] = collision.centV0M();

        values[4] = track1.phim() - track2.phim();
        if (values[4] > 1.5 * TMath::Pi())
          values[4] -= TMath::TwoPi();
        if (values[4] < -0.5 * TMath::Pi())
          values[4] += TMath::TwoPi();

        values[5] = collision.posZ();

        same->getTrackHist()->Fill(values, CorrelationContainer::kCFStepReconstructed);
        //mixed->getTrackHist()->Fill(values, CorrelationContainer::kCFStepReconstructed);
      }
    }
  }
};

// TODO: PID response task??

WorkflowSpec defineDataProcessing(ConfigContext const&)
{
  return WorkflowSpec{
    adaptAnalysisTask<FemtoscopyTask>("femtoscopy-task")};
}
