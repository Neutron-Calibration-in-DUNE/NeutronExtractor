/**
 * @file NeutronCapture.h
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2022-05-22
 */
#pragma once
#include <string>
#include <vector>
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/Simulation/SimPhotons.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RawData/BeamInfo.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/EndPoint2D.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larsim/Utils/TruthMatchUtils.h"

#include "DetectorGeometry.h"

namespace detinfo 
{
  class DetectorClocksData;
}
namespace neutron
{
    struct NeutronCaptureReco
    {
        std::vector<Double_t> SpacePointX;
        std::vector<Double_t> SpacePointY;
        std::vector<Double_t> SpacePointZ;
        std::vector<Double_t> SpacePointChiSqX;
        std::vector<Double_t> SpacePointChiSqY;
        std::vector<Double_t> SpacePointChiSqZ;
        std::vector<Double_t> SpacePointChiSq;

        std::vector<Int_t>    EdepTrackID;
        std::vector<Int_t>    NeutronTrackID;  
        std::vector<Int_t>    GammaTrackID;
        std::vector<Double_t> GammaEnergy;

        std::vector<Double_t> PeakTime;
        std::vector<Double_t> SigmaPeakTime;
        std::vector<Double_t> RMS;
        std::vector<Double_t> PeakAmplitude;
        std::vector<Double_t> SigmaPeakAmplitude;
        std::vector<Double_t> SummedADC;
    };

    class NeutronCapture
    {
    public:
        NeutronCapture();
        ~NeutronCapture();

        void processEvent(
            detinfo::DetectorClocksData const& clockData,
            const art::ValidHandle<std::vector<simb::MCParticle>>& mcParticles,
            const art::ValidHandle<std::vector<sim::SimChannel>>& mcChannels,
            const art::ValidHandle<std::vector<recob::SpacePoint>>& recoSpacePoints,
            const art::FindManyP<recob::Hit>& hitSpacePointAssn
        );
    private:
        /// ROOT output through art::TFileService
        /** We will save different TTrees to different TFiles specified 
         *  by the directories for each type.
         */ 
        art::ServiceHandle<art::TFileService> fTFileService;
        TTree *mNeutronCaptureRecoTree;
        // geometry information
        DetectorGeometry* fGeometry = DetectorGeometry::getInstance("NeutronCapture");

        NeutronCaptureReco mNeutronCaptureReco;
    };
}