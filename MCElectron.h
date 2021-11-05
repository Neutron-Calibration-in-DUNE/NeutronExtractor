/**
 * @file    Electron.h
 * @brief   A class for holding Electron information collected from LArSoft
 * @ingroup Electron
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
    class MCElectron
    {
    public:
        MCElectron();
        ~MCElectron();

        // get index from event_id and track_id
        Int_t getElectronIndex(Int_t eventId, Int_t trackId) {
            return fElectronMap[std::pair(eventId,trackId)];
        }

        // add a new Electron
        void addElectron(Int_t eventId, simb::MCParticle particle);
        // fill the mc Electron ttree
        void FillTTree();
    
    private:
        // geometry information
        DetectorGeometry* fGeometry = DetectorGeometry::getInstance("MCElectron");
        // ROOT 
        art::ServiceHandle<art::TFileService> fTFileService;
        TTree *fMCElectronTree;
    // MCParticle information for Electrons
    private:
        // number of Electrons
        size_t fNumberOfElectrons;
        // low level information for Electrons
        std::vector<Int_t> fEventId;
        std::vector<Int_t> fElectronTrackId;
        std::vector<Int_t> fElectronParentId;
        std::vector<Int_t> fElectronStatusCode;
        std::vector<Int_t> fElectronNumberOfDaughters;
        std::vector<std::vector<Int_t>> fElectronDaughters;
        std::vector<Int_t> fElectronNumberOfTrajectoryPoints;
        std::vector<std::vector<Double_t>> fElectronT;
        std::vector<std::vector<Double_t>> fElectronX;
        std::vector<std::vector<Double_t>> fElectronY;
        std::vector<std::vector<Double_t>> fElectronZ;
        std::vector<std::vector<Double_t>> fElectronE;
        std::vector<std::vector<Double_t>> fElectronPx;
        std::vector<std::vector<Double_t>> fElectronPy;
        std::vector<std::vector<Double_t>> fElectronPz;
        std::vector<std::string> fElectronProcess;
        std::vector<std::string> fElectronEndProcess;
        // beginning and ending volumes
        std::vector<Int_t> fElectronVolumeTypeBegin;
        std::vector<Int_t> fElectronVolumeTypeEnd;
        std::vector<std::string> fElectronVolumeNameBegin;
        std::vector<std::string> fElectronVolumeNameEnd;
        std::vector<std::string> fElectronMaterialNameBegin;
        std::vector<std::string> fElectronMaterialNameEnd;
        std::vector<Double_t> fElectronMaterialBegin;
        std::vector<Double_t> fElectronMaterialEnd;
        // distances and displacements
        std::vector<Double_t> fElectronTotalDistance;
        std::vector<std::vector<Double_t>> fElectronDistances;
        std::vector<Double_t> fElectronTotalDisplacement;
        // various time stamps for entering/exiting active volume
        std::vector<bool> fElectronEnteredActiveVolume;
        std::vector<bool> fElectronExitActiveVolume;
        std::vector<Int_t> fElectronEnteredActiveVolumeTime;
        std::vector<Int_t> fElectronExitActiveVolumeTime;
        // map from (event_id,track_id) -> index
        std::map<std::pair<Int_t,Int_t>, Int_t> fElectronMap;
        std::vector<std::vector<Int_t>> fElectronMapKeys;
    };
}
