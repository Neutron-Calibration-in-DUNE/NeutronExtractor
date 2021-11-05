/**
 * @file    Gamma.h
 * @brief   A class for holding Gamma information collected from LArSoft
 * @ingroup Gamma
 * @author  Nicholas Carrara (nmcarrara@ucdavis.edu),
**/
#pragma once

// art includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Utilities/Exception.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// special utility includes
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"
#include "art_root_io/TFileService.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

// necessary ROOT libraries
#include <TTree.h>
#include <TH1.h>
#include "TH1F.h"
#include "TGeoMaterial.h"
#include "TGeoElement.h"

// std includes
#include <string>
#include <vector>
#include <memory>

// local includes
#include "DetectorGeometry.h"
#include "Utilities.h"

namespace neutron 
{
    class MCGamma
    {
    public:
        MCGamma();
        ~MCGamma();

        // get index from event_id and track_id
        Int_t getGammaIndex(Int_t eventId, Int_t trackId) {
            return fGammaMap[std::pair(eventId,trackId)];
        }

        // add a new Gamma
        void addGamma(Int_t eventId, simb::MCParticle particle);
        // fill the mc Gamma ttree
        void FillTTree();
    
    private:
        // geometry information
        DetectorGeometry* fGeometry = DetectorGeometry::getInstance("MCGamma");
        // ROOT 
        art::ServiceHandle<art::TFileService> fTFileService;
        TTree *fMCGammaTree;
    // MCParticle information for Gammas
    private:
        // number of Gammas
        size_t fNumberOfGammas;
        // low level information for Gammas
        std::vector<Int_t> fEventId;
        std::vector<Int_t> fGammaTrackId;
        std::vector<Int_t> fGammaParentId;
        std::vector<Int_t> fGammaStatusCode;
        std::vector<Int_t> fGammaNumberOfDaughters;
        std::vector<std::vector<Int_t>> fGammaDaughters;
        std::vector<Int_t> fGammaNumberOfTrajectoryPoints;
        std::vector<std::vector<Double_t>> fGammaT;
        std::vector<std::vector<Double_t>> fGammaX;
        std::vector<std::vector<Double_t>> fGammaY;
        std::vector<std::vector<Double_t>> fGammaZ;
        std::vector<std::vector<Double_t>> fGammaE;
        std::vector<std::vector<Double_t>> fGammaPx;
        std::vector<std::vector<Double_t>> fGammaPy;
        std::vector<std::vector<Double_t>> fGammaPz;
        std::vector<std::string> fGammaProcess;
        std::vector<std::string> fGammaEndProcess;
        // beginning and ending volumes
        std::vector<Int_t> fGammaVolumeTypeBegin;
        std::vector<Int_t> fGammaVolumeTypeEnd;
        std::vector<std::string> fGammaVolumeNameBegin;
        std::vector<std::string> fGammaVolumeNameEnd;
        std::vector<std::string> fGammaMaterialNameBegin;
        std::vector<std::string> fGammaMaterialNameEnd;
        std::vector<Double_t> fGammaMaterialBegin;
        std::vector<Double_t> fGammaMaterialEnd;
        // distances and displacements
        std::vector<Double_t> fGammaTotalDistance;
        std::vector<std::vector<Double_t>> fGammaDistances;
        std::vector<Double_t> fGammaTotalDisplacement;
        // various time stamps for entering/exiting active volume
        std::vector<bool> fGammaEnteredActiveVolume;
        std::vector<bool> fGammaExitActiveVolume;
        std::vector<Int_t> fGammaEnteredActiveVolumeTime;
        std::vector<Int_t> fGammaExitActiveVolumeTime;
        // map from (event_id,track_id) -> index
        std::map<std::pair<Int_t,Int_t>, Int_t> fGammaMap;
        std::vector<std::vector<Int_t>> fGammaMapKeys;
    };
}
