/**
 * @file    Neutron.h
 * @brief   A class for holding neutron information collected from LArSoft
 * @ingroup Neutron
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
    class MCNeutron
    {
    public:
        MCNeutron();
        ~MCNeutron();

        // get index from event_id and track_id
        Int_t getNeutronIndex(Int_t eventId, Int_t trackId) {
            return fNeutronMap[std::pair(eventId,trackId)];
        }

        // add a new neutron
        void addNeutron(Int_t eventId, simb::MCParticle particle);
        // add a gamma
        void addGamma(Int_t eventId, simb::MCParticle gamma);
        // fill the mc neutron ttree
        void FillTTree();
    
    private:
        // geometry information
        DetectorGeometry* fGeometry = DetectorGeometry::getInstance("MCNeutron");
        // ROOT 
        art::ServiceHandle<art::TFileService> fTFileService;
        TTree *fMCNeutronTree;
    // MCParticle information for neutrons
    private:
        // number of neutrons
        size_t fNumberOfNeutrons;
        // low level information for neutrons
        std::vector<Int_t> fEventId;
        std::vector<Int_t> fNeutronTrackId;
        std::vector<Int_t> fNeutronParentId;
        std::vector<Int_t> fNeutronStatusCode;
        std::vector<Int_t> fNeutronNumberOfDaughters;
        std::vector<std::vector<Int_t>> fNeutronDaughters;
        std::vector<Int_t> fNeutronNumberOfTrajectoryPoints;
        std::vector<std::vector<Double_t>> fNeutronT;
        std::vector<std::vector<Double_t>> fNeutronX;
        std::vector<std::vector<Double_t>> fNeutronY;
        std::vector<std::vector<Double_t>> fNeutronZ;
        std::vector<std::vector<Double_t>> fNeutronE;
        std::vector<std::vector<Double_t>> fNeutronPx;
        std::vector<std::vector<Double_t>> fNeutronPy;
        std::vector<std::vector<Double_t>> fNeutronPz;
        std::vector<std::string> fNeutronProcess;
        std::vector<std::string> fNeutronEndProcess;
        // list of inelastics tied to a given neutron
        std::vector<std::vector<Int_t>> fNeutronInelastic;
        // beginning and ending volumes
        std::vector<Int_t> fNeutronVolumeTypeBegin;
        std::vector<Int_t> fNeutronVolumeTypeEnd;
        std::vector<std::string> fNeutronVolumeNameBegin;
        std::vector<std::string> fNeutronVolumeNameEnd;
        std::vector<std::string> fNeutronMaterialNameBegin;
        std::vector<std::string> fNeutronMaterialNameEnd;
        std::vector<Double_t> fNeutronMaterialBegin;
        std::vector<Double_t> fNeutronMaterialEnd;
        // distances and displacements
        std::vector<Double_t> fNeutronTotalDistance;
        std::vector<std::vector<Double_t>> fNeutronDistances;
        std::vector<Double_t> fNeutronTotalDisplacement;
        // various time stamps for entering/exiting active volume
        std::vector<bool> fNeutronEnteredActiveVolume;
        std::vector<bool> fNeutronExitActiveVolume;
        std::vector<Int_t> fNeutronEnteredActiveVolumeTime;
        std::vector<Int_t> fNeutronExitActiveVolumeTime;
        // map from (event_id,track_id) -> index
        std::map<std::pair<Int_t,Int_t>, Int_t> fNeutronMap;
        std::vector<std::vector<Int_t>> fNeutronMapKeys;
    // MCParticle information for gammas that come from captures
    private:
        std::vector<Int_t> fNumberOfCaptureGammas;
        std::vector<std::vector<Double_t>> fCaptureGammaInitialEnergy;
        // initial and final position of the captured gammas
        std::vector<std::vector<Double_t>> fCaptureGammaInitialX;
        std::vector<std::vector<Double_t>> fCaptureGammaFinalX;
        std::vector<std::vector<Double_t>> fCaptureGammaInitialY;
        std::vector<std::vector<Double_t>> fCaptureGammaFinalY;
        std::vector<std::vector<Double_t>> fCaptureGammaInitialZ;
        std::vector<std::vector<Double_t>> fCaptureGammaFinalZ;
        // initial momentum of the captured gammas
        std::vector<std::vector<Double_t>> fCaptureGammaInitialPx;
        std::vector<std::vector<Double_t>> fCaptureGammaInitialPy;
        std::vector<std::vector<Double_t>> fCaptureGammaInitialPz;

    };

}
