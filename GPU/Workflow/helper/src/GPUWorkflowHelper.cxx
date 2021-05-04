// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "GPUWorkflowHelper/GPUWorkflowHelper.h"
#include "DataFormatsTRD/RecoInputContainer.h"
#include "ITStracking/IOUtils.h"
#include "DataFormatsTPC/WorkflowHelper.h"
using namespace o2::globaltracking;
using namespace o2::gpu;

struct GPUWorkflowHelper::tmpDataContainer {
  std::vector<o2::BaseCluster<float>> ITSClustersArray;
};

std::shared_ptr<const GPUWorkflowHelper::tmpDataContainer> GPUWorkflowHelper::fillIOPtr(GPUTrackingInOutPointers& ioPtr, const o2::globaltracking::RecoContainer& recoCont, bool useMC, const GPUCalibObjectsConst* calib, o2::dataformats::GlobalTrackID::mask_t maskCl, o2::dataformats::GlobalTrackID::mask_t maskTrk, o2::dataformats::GlobalTrackID::mask_t maskMatch)
{
  auto retVal = std::make_shared<tmpDataContainer>();

  if (maskCl[GID::ITS] && ioPtr.nItsClusters == 0) {
    const auto& ITSClusterROFRec = recoCont.getITSClustersROFRecords<o2::itsmft::ROFRecord>();
    const auto& clusITS = recoCont.getITSClusters<o2::itsmft::CompClusterExt>();
    if (clusITS.size() && ITSClusterROFRec.size()) {
      if (calib && calib->itsPatternDict) {
        const auto& patterns = recoCont.getITSClustersPatterns();
        auto pattIt = patterns.begin();
        retVal->ITSClustersArray.reserve(clusITS.size());
        o2::its::ioutils::convertCompactClusters(clusITS, pattIt, retVal->ITSClustersArray, *calib->itsPatternDict);
        ioPtr.itsClusters = retVal->ITSClustersArray.data();
      }
      ioPtr.nItsClusters = clusITS.size();
      ioPtr.itsCompClusters = clusITS.data();
      ioPtr.nItsClusterROF = ITSClusterROFRec.size();
      ioPtr.itsClusterROF = ITSClusterROFRec.data();
      if (useMC) {
        const auto& ITSClsLabels = recoCont.mcITSClusters.get();
        ioPtr.itsClusterMC = ITSClsLabels;
      }
    }
  }
  if (maskTrk[GID::ITS] && ioPtr.nItsTracks == 0) {
    const auto& ITSTracksArray = recoCont.getITSTracks<o2::its::TrackITS>();
    const auto& ITSTrackROFRec = recoCont.getITSTracksROFRecords<o2::itsmft::ROFRecord>();
    if (ITSTracksArray.size() && ITSTrackROFRec.size()) {
      const auto& ITSTrackClusIdx = recoCont.getITSTracksClusterRefs();
      ioPtr.nItsTracks = ITSTracksArray.size();
      ioPtr.itsTracks = ITSTracksArray.data();
      ioPtr.itsTrackClusIdx = ITSTrackClusIdx.data();
      ioPtr.nItsTrackROF = ITSTrackROFRec.size();
      ioPtr.itsTrackROF = ITSTrackROFRec.data();
      if (useMC) {
        const auto& ITSTrkLabels = recoCont.getITSTracksMCLabels();
        ioPtr.itsTrackMC = ITSTrkLabels.data();
      }
    }
  }

  if (maskTrk[GID::ITSTPC] && ioPtr.nTracksTPCITSO2 == 0) {
    const auto& trkITSTPC = recoCont.getTPCITSTracks<o2::dataformats::TrackTPCITS>();
    if (trkITSTPC.size()) {
      ioPtr.nTracksTPCITSO2 = trkITSTPC.size();
      ioPtr.tracksTPCITSO2 = trkITSTPC.data();
    }
  }

  if (maskCl[GID::TOF] && ioPtr.nTOFClusters == 0) {
    const auto& tofClusters = recoCont.getTOFClusters<o2::tof::Cluster>();
    if (tofClusters.size()) {
      ioPtr.nTOFClusters = tofClusters.size();
      ioPtr.tofClusters = tofClusters.data();
    }
  }

  if (maskMatch[GID::TOF] && ioPtr.nTOFMatches == 0) {
    const auto& tofMatches = recoCont.getTOFMatches<o2::dataformats::MatchInfoTOF>();
    if (tofMatches.size()) {
      ioPtr.nTOFMatches = tofMatches.size();
      ioPtr.tofMatches = tofMatches.data();
    }
  }

  if (maskMatch[GID::TPCTOF] && ioPtr.nTPCTOFMatches == 0) {
    const auto& tpctofMatches = recoCont.getTPCTOFMatches<o2::dataformats::MatchInfoTOF>();
    if (tpctofMatches.size()) {
      ioPtr.nTPCTOFMatches = tpctofMatches.size();
      ioPtr.tpctofMatches = tpctofMatches.data();
    }
  }

  if (maskCl[GID::TPC] && ioPtr.clustersNative == nullptr) {
    ioPtr.clustersNative = &recoCont.inputsTPCclusters->clusterIndex;
  }

  if (maskTrk[GID::TPC] && ioPtr.nOutputTracksTPCO2 == 0) {
    const auto& tpcTracks = recoCont.getTPCTracks<o2::tpc::TrackTPC>();
    const auto& tpcClusRefs = recoCont.getTPCTracksClusterRefs();
    ioPtr.outputTracksTPCO2 = tpcTracks.data();
    ioPtr.nOutputTracksTPCO2 = tpcTracks.size();
    ioPtr.outputClusRefsTPCO2 = tpcClusRefs.data();
    ioPtr.nOutputClusRefsTPCO2 = tpcClusRefs.size();
    if (useMC) {
      const auto& tpcTracksMC = recoCont.getTPCTracksMCLabels();
      ioPtr.outputTracksTPCO2MC = tpcTracksMC.data();
    }
  }

  if (maskCl[GID::TRD]) {
    // o2::trd::getRecoInputContainer(pc, &ioPtr, &recoCont); // TODO: use this helper here
  }

  return std::move(retVal);
}