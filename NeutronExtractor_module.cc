/**
 * @file    NeutronExtractor_module.cc
 * @brief   A module for extracting truth/reco information about G4 particle trajectories
 *          and conducting some standard analysis tasks. 
 *          Generated at Mon Oct 11 11:21:12 2021 using cetskelgen
 * @ingroup NeutronExtractor
 * @author  Nicholas Carrara (nmcarrara@ucdavis.edu),
 *          Yashwanth Bezawada
**/

// art framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Persistency/Common/PtrMaker.h"
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
#include "larsim/Utils/TruthMatchUtils.h"
#include "lardata/ArtDataHelper/TrackUtils.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

// necessary ROOT libraries
#include <TTree.h>
#include <TH1.h>
#include "TH1F.h"
#include "TGeoMaterial.h"
#include "TGeoElement.h"

// C includes
#include <cmath>

// local includes
#include "DetectorGeometry.h"
#include "Voxelizer.h"

namespace neutron {
    class NeutronExtractor;
}

namespace neutron
{
    struct NeutronTrajectories
    {
        /// struct for storing neutron trajectories for an event.
        /** This struct contains some quantities of interest beyond
         *  the trajectories, such as:
         *      1. When/where the neutron first enters the LAr.
         *      2. When/where the neutron exits the LAr after entering.
         *      3. The total distance traveled inside the LAr.
         *      4. How many different materials the neutron interacts with inside the LAr.
         *      5. The overall dE/dx of the trajectory.
         *      6. The dE/dx just inside the LAr.
         * 
         */ 
        Int_t event_id;
        std::vector<Int_t> track_id;
        std::vector<Int_t> mother;
        std::vector<std::string> process;
        std::vector<std::string> end_process;
        std::vector<std::vector<Double_t>> t;
        std::vector<std::vector<Double_t>> x;
        std::vector<std::vector<Double_t>> y;
        std::vector<std::vector<Double_t>> z;
        std::vector<std::vector<Double_t>> energy;
        std::vector<std::vector<Double_t>> px;
        std::vector<std::vector<Double_t>> py;
        std::vector<std::vector<Double_t>> pz;
        std::vector<std::vector<std::string>> volume_name; 
        std::vector<std::vector<std::string>> material_name;
        std::vector<std::vector<Double_t>> material;

        std::vector<Double_t> enter_lar_t;
        std::vector<Double_t> enter_lar_x;
        std::vector<Double_t> enter_lar_y;
        std::vector<Double_t> enter_lar_z;
        std::vector<Double_t> exit_lar_t;
        std::vector<Double_t> exit_lar_x;
        std::vector<Double_t> exit_lar_y;
        std::vector<Double_t> exit_lar_z;
        
        std::vector<Double_t> total_distance;
        std::vector<Double_t> internal_distance;
        std::vector<Double_t> lar_distance;
        
        std::vector<std::vector<Double_t>> total_dedx;
        std::vector<std::vector<Double_t>> internal_dedx;
        std::vector<std::vector<Double_t>> lar_dedx;

        std::vector<Int_t> num_total_steps;

        NeutronTrajectories(Int_t event) : event_id(event) {}
    };
    struct NeutronList
    {
        // neutron capture related information
        Int_t event_id;
        Int_t primary_neutrons;
        std::vector<Int_t> neutron_ids;
        std::vector<Double_t> neutron_capture_x;
        std::vector<Double_t> neutron_capture_y;
        std::vector<Double_t> neutron_capture_z;
        std::vector<Int_t> gamma_ids;
        std::vector<Int_t> gamma_neutron_ids;
        std::vector<Double_t> gamma_energy;
        std::vector<Double_t> gamma_electron_energy;
        std::vector<Double_t> gamma_edep_energy;

        std::vector<Int_t> electron_ids;
        std::vector<Int_t> electron_parent;
        std::vector<Int_t> electron_gamma_ids;
        std::vector<Int_t> electron_neutron_ids;
        std::vector<Double_t> electron_energy;
        std::vector<Int_t> edep_parent;

        std::vector<Int_t> edep_neutron_ids;
        std::vector<Int_t> edep_gamma_ids;
        std::vector<Double_t> edep_energy; 
        std::vector<Int_t> edep_num_electrons;
        std::vector<Double_t> edep_x;
        std::vector<Double_t> edep_y;
        std::vector<Double_t> edep_z;

        std::vector<Int_t> hit_track_ids;
        std::vector<Int_t> hit_channel_ids;
        std::vector<Int_t> hit_view;

        std::vector<Int_t> space_point_ids;
        std::vector<Double_t> space_point_x;
        std::vector<Double_t> space_point_y;
        std::vector<Double_t> space_point_z;

        NeutronList(Int_t event) : event_id(event) {}
    };

    struct MuonList
    {
        Int_t event_id;
        // muon related information
        Int_t primary_muons;
        std::vector<Int_t> muon_ids;
        std::vector<Int_t> muon_edep_ids;
        std::vector<Double_t> muon_edep_energy;
        std::vector<Double_t> muon_edep_num_electrons;
        std::vector<Double_t> muon_edep_x;
        std::vector<Double_t> muon_edep_y;
        std::vector<Double_t> muon_edep_z;

        MuonList(Int_t event) : event_id(event) {}
    };

    class NeutronExtractor : public art::EDAnalyzer
    {
    public:
        struct Config
        {
            /// Neutron Trajectory Information
            /** We will save the entire trajectory of each neutron if this tag is 
             *  set to true.  This will include information at each step as well
             *  as statistics from the trajectory as a whole.  These trajectories
             *  are stored in the sim::MCParticle class, so little has to be done
             *  on the user side to obtain all the relevant information.
             */ 
            fhicl::Atom<bool> SavePrimaryNeutronTrajectories
            {
                fhicl::Name("SavePrimaryNeutronTrajectories"),
                fhicl::Comment("whether or not to save the neutron trajectories in a TFile.")
            };
            fhicl::Atom<bool> SavePrimaryNeutrondEdx
            {
                fhicl::Name("SavePrimaryNeutrondEdx"),
                fhicl::Comment("whether or not to save the neutron dEdx values.")
            };
            /// larg4 Energy Deposition Information
            /**
             * 
             * 
             */ 
            fhicl::Atom<art::InputTag> LArGeantProducerLabel
            {
                fhicl::Name("LArGeantProducerLabel"),
                fhicl::Comment("tag of the input data product with the largeant side of the simulation")
            };
            fhicl::Atom<art::InputTag> LArGeantEnergyDepositProducerLabel
            {
                fhicl::Name("LArGeantEnergyDepositProducerLabel"),
                fhicl::Comment("tag of the input data product with the largeant side of the simulation")
            };
            fhicl::Atom<art::InputTag> IonAndScintProducerLabel
            {
                fhicl::Name("IonAndScintProducerLabel"),
                fhicl::Comment("tag of the input data product with the ionization and scintillation simulation")
            };
            fhicl::Atom<bool> FindHits
            {
                fhicl::Name("FindHits"),
                fhicl::Comment("whether to look for hits in the data file")
            };
            fhicl::Atom<bool> GenerateNeutronCosmicVoxels
            {
                fhicl::Name("GenerateNeutronCosmicVoxels"),
                fhicl::Comment("whether to generate neutron/cosmic training set")
            };
            fhicl::Atom<Double_t> VoxelSize
            {
                fhicl::Name("VoxelSize"),
                fhicl::Comment("size of the voxels to use in mm")
            };
            fhicl::Atom<bool> DiscretizeVoxelFeatures
            {
                fhicl::Name("DiscretizeVoxelFeatures"),
                fhicl::Comment("whether to discretize voxel features")
            };
            fhicl::Atom<bool> UseMixedVoxelLabeling
            {
                fhicl::Name("UseMixedVoxelLabeling"),
                fhicl::Comment("whether to use mixed labeling for voxels")
            };
            fhicl::Atom<art::InputTag> HitFinderProducerLabel
            {
                fhicl::Name("HitFinderProducerLabel"),
                fhicl::Comment("tag of the data product which contains the hits")
            };
            fhicl::Atom<art::InputTag> SpacePointProducerLabel
            {
                fhicl::Name("SpacePointProducerLabel"),
                fhicl::Comment("tag of the data product which contains the space points")
            };
            fhicl::Atom<art::InputTag> OutputFile
            {
                fhicl::Name("OutputFile"),
                fhicl::Comment("name of the file to output the neutron statistics to")
            };
        };
    public:
        using Parameters = art::EDAnalyzer::Table<Config>;
        explicit NeutronExtractor(Parameters const& config);
        NeutronExtractor(NeutronExtractor const&) = delete;
        NeutronExtractor(NeutronExtractor&&) = delete;
        NeutronExtractor& operator=(NeutronExtractor const&) = delete;
        NeutronExtractor& operator=(NeutronExtractor&&) = delete;

        // required EDAnalyzer functions
        void analyze(art::Event const& event) override;
        void beginJob() override;
        void endJob() override;

        // special functions
        void FillTTree();
        bool checkListOfNeutrons(Int_t eventId, Int_t trackId);
        bool checkListOfGammas(Int_t eventId, Int_t trackId);
        bool checkListOfElectrons(Int_t eventId, Int_t trackId);

    private:
        /// Neutron Trajectory Information
        bool fSavePrimaryNeutronTrajectories;
        bool fSavePrimaryNeutrondEdx;
        /// Energy Deposition Information
        art::InputTag fLArGeantProducerLabel;
        art::InputTag fLArGeantEnergyDepositProducerLabel;
        art::InputTag fIonAndScintProducerLabel;
        bool fFindHits;
        bool fGenerateNeutronCosmicVoxels;
        Double_t fVoxelSize;
        bool fDiscretizeVoxelFeatures;
        bool fUseMixedVoxelLabeling;
        art::InputTag fHitFinderProducerLabel;
        art::InputTag fSpacePointProducerLabel;
        art::InputTag fOutputFileArt;
        
        /// ROOT output through art::TFileService
        /** We will save different TTrees to different TFiles specified 
         *  by the directories for each type.  For example, the neutron
         *  trajectory information will be stored in a TFile. 
         * 
         */ 
        art::ServiceHandle<art::TFileService> fTFileService;
        /// TTrees for different tasks
        TTree *fMetaTree;
        TTree *fNeutronTrajectoryTree;
        TTree *fNeutronTree;
        TTree *fMuonTree;
        TTree *fVoxelTree;
        // event variables
        int fRun;
        int fSubRun;
        int fEvent;

        // geometry information
        DetectorGeometry* fGeometry = DetectorGeometry::getInstance("NeutronExtractor");
        // voxelizer
        Voxelizer fVoxelizer;

        // number of events
        Int_t fNumberOfEvents;

        /// temporary structs for holding event information
        NeutronTrajectories fTempNeutronTrajectories;
        NeutronList         fTempNeutronList;
        MuonList            fTempMuonList;
        Voxels              fTempVoxels;

        std::vector<Int_t> fNumberOfNeutronsPerEvent;
    };

    // constructor
    NeutronExtractor::NeutronExtractor(Parameters const& config)
    : EDAnalyzer(config)
    , fSavePrimaryNeutronTrajectories(config().SavePrimaryNeutronTrajectories())
    , fSavePrimaryNeutrondEdx(config().SavePrimaryNeutrondEdx())
    , fLArGeantProducerLabel(config().LArGeantProducerLabel())
    , fLArGeantEnergyDepositProducerLabel(config().LArGeantEnergyDepositProducerLabel())
    , fIonAndScintProducerLabel(config().IonAndScintProducerLabel())
    , fFindHits(config().FindHits())
    , fGenerateNeutronCosmicVoxels(config().GenerateNeutronCosmicVoxels())
    , fVoxelSize(config().VoxelSize())
    , fDiscretizeVoxelFeatures(config().DiscretizeVoxelFeatures())
    , fUseMixedVoxelLabeling(config().UseMixedVoxelLabeling())
    , fHitFinderProducerLabel(config().HitFinderProducerLabel())
    , fSpacePointProducerLabel(config().SpacePointProducerLabel())
    , fOutputFileArt(config().OutputFile())
    , fTempNeutronTrajectories(0)
    , fTempNeutronList(0)
    , fTempMuonList(0)
    , fTempVoxels(0)
    {
        /// set up TTrees
        fMetaTree    = fTFileService->make<TTree>("meta", "meta");
        fNeutronTrajectoryTree = fTFileService->make<TTree>("neutron_trajectories", "neutron_trajectories");
        fNeutronTree = fTFileService->make<TTree>("neutron", "neutron");
        fMuonTree    = fTFileService->make<TTree>("muon", "muon");
        fVoxelTree   = fTFileService->make<TTree>("voxels", "voxels");

        // construct voxelizer
        fVoxelizer.setBoundingBox(fGeometry->GetTotalActiveTPCBox());

        consumes<std::vector<simb::MCParticle>>(fLArGeantProducerLabel);
        consumes<std::vector<sim::SimEnergyDeposit>>(fLArGeantEnergyDepositProducerLabel);
        if (fFindHits)
        {
            consumes<std::vector<recob::Hit>>(fHitFinderProducerLabel);
            consumes<std::vector<recob::SpacePoint>>(fSpacePointProducerLabel);
        }
        // set up branches
        fNeutronTrajectoryTree->Branch("event_id",  &fTempNeutronTrajectories.event_id);
        fNeutronTrajectoryTree->Branch("track_id",  &fTempNeutronTrajectories.track_id);
        fNeutronTrajectoryTree->Branch("mother",    &fTempNeutronTrajectories.mother);
        fNeutronTrajectoryTree->Branch("process",   &fTempNeutronTrajectories.process);
        fNeutronTrajectoryTree->Branch("end_process", &fTempNeutronTrajectories.end_process);
        fNeutronTrajectoryTree->Branch("t", &fTempNeutronTrajectories.t);
        fNeutronTrajectoryTree->Branch("x", &fTempNeutronTrajectories.x);
        fNeutronTrajectoryTree->Branch("y", &fTempNeutronTrajectories.y);
        fNeutronTrajectoryTree->Branch("z", &fTempNeutronTrajectories.z);
        fNeutronTrajectoryTree->Branch("energy", &fTempNeutronTrajectories.energy);
        fNeutronTrajectoryTree->Branch("px", &fTempNeutronTrajectories.px);
        fNeutronTrajectoryTree->Branch("py", &fTempNeutronTrajectories.py);
        fNeutronTrajectoryTree->Branch("pz", &fTempNeutronTrajectories.pz);
        fNeutronTrajectoryTree->Branch("volume_name", &fTempNeutronTrajectories.volume_name);
        fNeutronTrajectoryTree->Branch("material_name", &fTempNeutronTrajectories.material_name);
        fNeutronTrajectoryTree->Branch("material", &fTempNeutronTrajectories.material);
        fNeutronTrajectoryTree->Branch("enter_lar_t", &fTempNeutronTrajectories.enter_lar_t);
        fNeutronTrajectoryTree->Branch("enter_lar_x", &fTempNeutronTrajectories.enter_lar_x);
        fNeutronTrajectoryTree->Branch("enter_lar_y", &fTempNeutronTrajectories.enter_lar_y);
        fNeutronTrajectoryTree->Branch("enter_lar_z", &fTempNeutronTrajectories.enter_lar_z);
        fNeutronTrajectoryTree->Branch("exit_lar_t", &fTempNeutronTrajectories.exit_lar_t);
        fNeutronTrajectoryTree->Branch("exit_lar_x", &fTempNeutronTrajectories.exit_lar_x);
        fNeutronTrajectoryTree->Branch("exit_lar_y", &fTempNeutronTrajectories.exit_lar_y);
        fNeutronTrajectoryTree->Branch("exit_lar_z", &fTempNeutronTrajectories.exit_lar_z);
        fNeutronTrajectoryTree->Branch("total_distance", &fTempNeutronTrajectories.total_distance);
        fNeutronTrajectoryTree->Branch("internal_distance", &fTempNeutronTrajectories.internal_distance);
        fNeutronTrajectoryTree->Branch("lar_distance", &fTempNeutronTrajectories.lar_distance);
        fNeutronTrajectoryTree->Branch("total_dedx", &fTempNeutronTrajectories.total_dedx);
        fNeutronTrajectoryTree->Branch("internal_dedx", &fTempNeutronTrajectories.internal_dedx);
        fNeutronTrajectoryTree->Branch("lar_dedx", &fTempNeutronTrajectories.lar_dedx);
        fNeutronTrajectoryTree->Branch("num_total_steps", &fTempNeutronTrajectories.num_total_steps);

        fNeutronTree->Branch("event_id", &fTempNeutronList.event_id);
        fNeutronTree->Branch("primary_neutrons", &fTempNeutronList.primary_neutrons);
        fNeutronTree->Branch("neutron_ids", &fTempNeutronList.neutron_ids);
        fNeutronTree->Branch("neutron_capture_x", &fTempNeutronList.neutron_capture_x);
        fNeutronTree->Branch("neutron_capture_y", &fTempNeutronList.neutron_capture_y);
        fNeutronTree->Branch("neutron_capture_z", &fTempNeutronList.neutron_capture_z);
        fNeutronTree->Branch("gamma_ids", &fTempNeutronList.gamma_ids);
        fNeutronTree->Branch("gamma_neutron_ids", &fTempNeutronList.gamma_neutron_ids);
        fNeutronTree->Branch("gamma_energy", &fTempNeutronList.gamma_energy);
        fNeutronTree->Branch("gamma_electron_energy", &fTempNeutronList.gamma_electron_energy);
        fNeutronTree->Branch("gamma_edep_energy", &fTempNeutronList.gamma_edep_energy);
        fNeutronTree->Branch("electron_ids", &fTempNeutronList.electron_ids);
        fNeutronTree->Branch("electron_neutron_ids", &fTempNeutronList.electron_neutron_ids);
        fNeutronTree->Branch("electron_gamma_ids", &fTempNeutronList.electron_gamma_ids);
        fNeutronTree->Branch("electron_parent", &fTempNeutronList.electron_parent);
        fNeutronTree->Branch("electron_energy", &fTempNeutronList.electron_energy);
        fNeutronTree->Branch("edep_parent", &fTempNeutronList.edep_parent);
        fNeutronTree->Branch("edep_neutron_ids", &fTempNeutronList.edep_neutron_ids);
        fNeutronTree->Branch("edep_gamma_ids", &fTempNeutronList.edep_gamma_ids);
        fNeutronTree->Branch("edep_energy", &fTempNeutronList.edep_energy);
        fNeutronTree->Branch("edep_num_electrons", &fTempNeutronList.edep_num_electrons);
        fNeutronTree->Branch("edep_x", &fTempNeutronList.edep_x);
        fNeutronTree->Branch("edep_y", &fTempNeutronList.edep_y);
        fNeutronTree->Branch("edep_z", &fTempNeutronList.edep_z);
        fNeutronTree->Branch("hit_track_ids", &fTempNeutronList.hit_track_ids);
        fNeutronTree->Branch("hit_view", &fTempNeutronList.hit_view);
        fNeutronTree->Branch("space_point_ids", &fTempNeutronList.space_point_ids);
        fNeutronTree->Branch("space_point_x", &fTempNeutronList.space_point_x);
        fNeutronTree->Branch("space_point_y", &fTempNeutronList.space_point_y);
        fNeutronTree->Branch("space_point_z", &fTempNeutronList.space_point_z);

        fMuonTree->Branch("event_id", &fTempMuonList.event_id);
        fMuonTree->Branch("primary_muons", &fTempMuonList.primary_muons);
        fMuonTree->Branch("muon_ids", &fTempMuonList.muon_ids);
        fMuonTree->Branch("muon_edep_ids", &fTempMuonList.muon_edep_ids);
        fMuonTree->Branch("muon_edep_energy", &fTempMuonList.muon_edep_energy);
        fMuonTree->Branch("muon_edep_num_electrons", &fTempMuonList.muon_edep_num_electrons);
        fMuonTree->Branch("muon_edep_x", &fTempMuonList.muon_edep_x);
        fMuonTree->Branch("muon_edep_y", &fTempMuonList.muon_edep_y);
        fMuonTree->Branch("muon_edep_z", &fTempMuonList.muon_edep_z);

        fVoxelTree->Branch("x_min", &fTempVoxels.x_min);
        fVoxelTree->Branch("x_max", &fTempVoxels.x_max);
        fVoxelTree->Branch("y_min", &fTempVoxels.y_min);
        fVoxelTree->Branch("y_max", &fTempVoxels.y_max);
        fVoxelTree->Branch("z_min", &fTempVoxels.z_min);
        fVoxelTree->Branch("z_max", &fTempVoxels.z_max);
        fVoxelTree->Branch("voxel_size", &fTempVoxels.voxel_size);
        fVoxelTree->Branch("num_voxels_x", &fTempVoxels.num_voxels_x);
        fVoxelTree->Branch("num_voxels_y", &fTempVoxels.num_voxels_y);
        fVoxelTree->Branch("num_voxels_z", &fTempVoxels.num_voxels_z);
        fVoxelTree->Branch("x_id", &fTempVoxels.x_id);
        fVoxelTree->Branch("y_id", &fTempVoxels.y_id);
        fVoxelTree->Branch("z_id", &fTempVoxels.z_id);
        fVoxelTree->Branch("values", &fTempVoxels.values);
        fVoxelTree->Branch("labels", &fTempVoxels.labels);
        fVoxelTree->Branch("neutron_edep_ids", &fTempVoxels.neutron_edep_ids);
        fVoxelTree->Branch("muon_edep_ids", &fTempVoxels.muon_edep_ids);
    }

    // analyze function
    void NeutronExtractor::analyze(art::Event const& event)
    {
        if (event.isRealData())
        {
            // If we are looking at real data, then we need to stop the analysis
            // and back out.
            throw cet::exception("NeutronExtractor")
                << " Event contains real data - "
                << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
        }
        // define a "handle" to point to a vector of the objects
        art::Handle<std::vector<simb::MCParticle>> particleHandle;
        if (!event.getByLabel(fLArGeantProducerLabel, particleHandle))
        {
            // if there are no particles for the event truth, then
            // we are in big trouble haha.  throw an exception
            throw cet::exception("NeutronExtractor")
                << " No simb::MCParticle objects in this event - "
                << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
        }

        // get the event meta data 
        fRun    = event.run();
        fSubRun = event.subRun();
        fEvent  = event.id().event();

        fNumberOfEvents++;
        fNumberOfNeutronsPerEvent.emplace_back(0);

        // create a new event list
        NeutronTrajectories neutronTrajectories(fNumberOfEvents-1);
        NeutronList         neutronList(fNumberOfEvents-1);
        MuonList            muonList(fNumberOfEvents-1);
        neutronList.primary_neutrons = 0;
        muonList.primary_muons = 0;

        // get clock information
        auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(event);

        // get the list of MC particles from Geant4
        auto mcParticles = event.getValidHandle<std::vector<simb::MCParticle>>(fLArGeantProducerLabel);
        // iterate over all MC particles and grab all neutrons, all gammas
        // which come from neutron captures regardless of where they happen and
        // all electrons generated from the gammas.
        if (mcParticles.isValid())
        {
            for (auto particle : *mcParticles)
            {
                // check if the particle is a neutron
                if (particle.PdgCode() == 2112)
                {
                    if (particle.Mother() == 0)
                    {
                        neutronList.primary_neutrons += 1;
                        // First grab basic info from trajectories, 
                        // track_id, mother, process and end_process.
                        neutronTrajectories.track_id.emplace_back(particle.TrackId());
                        neutronTrajectories.mother.emplace_back(particle.Mother());
                        neutronTrajectories.process.emplace_back(particle.Process());
                        neutronTrajectories.end_process.emplace_back(particle.EndProcess());

                        // Collect the number of trajectory points
                        auto numTrajectoryPoints = particle.NumberTrajectoryPoints();
                        neutronTrajectories.num_total_steps.emplace_back(numTrajectoryPoints);

                        // // Initialize the full trajectory arrays
                        // std::vector<Double_t> t(numTrajectoryPoints);
                        // std::vector<Double_t> x(numTrajectoryPoints);
                        // std::vector<Double_t> y(numTrajectoryPoints);
                        // std::vector<Double_t> z(numTrajectoryPoints);
                        // std::vector<Double_t> energy(numTrajectoryPoints);
                        // std::vector<Double_t> px(numTrajectoryPoints);
                        // std::vector<Double_t> py(numTrajectoryPoints);
                        // std::vector<Double_t> pz(numTrajectoryPoints);
                        // std::vector<std::string> volume_name(numTrajectoryPoints);
                        // std::vector<std::string> material_name(numTrajectoryPoints);
                        // std::vector<Double_t> material(numTrajectoryPoints);

                        // // keep track of different particle trajectory statistics
                        // std::vector<Double_t> enter_lar = {-1.,-1.,-1.,-1.};
                        // std::vector<Double_t> exit_lar = {-1.,-1.,-1.,-1.};
                        // Double_t total_distance = 0;
                        // Double_t internal_distance = 0;
                        // Double_t lar_distance = 0;          
                        // std::vector<Double_t> total_dedx(numTrajectoryPoints);
                        // std::vector<Double_t> internal_dedx;
                        // std::vector<Double_t> lar_dedx;
                        
                        // // statistics helpers
                        // bool entered_lar = false;
                        // bool exited_lar = false;

                        // // Iterating over the trajectory points
                        // for (size_t i = 0; i < numTrajectoryPoints; i++)
                        // {
                        //     t[i] = particle.T(i);
                        //     x[i] = particle.Vx(i);
                        //     y[i] = particle.Vy(i);
                        //     z[i] = particle.Vz(i);
                        //     energy[i] = particle.E(i);
                        //     px[i] = particle.Px(i);
                        //     py[i] = particle.Py(i);
                        //     pz[i] = particle.Pz(i);
                        //     // some volume information for the current step
                        //     DetectorVolume volume = fGeometry->getVolume(
                        //         particle.Vx(i), particle.Vy(i), particle.Vz(i)
                        //     );
                        //     volume_name[i] = volume.volume_name;
                        //     material_name[i] = volume.material_name;
                        //     material[i] = volume.material;
                        //     // Distance  and dEdx calculations
                        //     Double_t dx = 0;
                        //     Double_t dE = 0;
                        //     if (i != 0)
                        //     {
                        //         dx = sqrt(
                        //             (x[i] - x[i-1])*(x[i] - x[i-1])
                        //         + (y[i] - y[i-1])*(y[i] - y[i-1])
                        //         + (z[i] - z[i-1])*(z[i] - z[i-1])
                        //         ); 
                        //         dE = energy[i] - energy[i-1];
                        //     }
                        //     // calculate dedx
                        //     Double_t dEdx = 0;
                        //     if (i != 0) {
                        //         dEdx = dE/dx;
                        //     }
                        //     // append total distance and dedx values
                        //     total_distance += dx;
                        //     total_dedx[i] = dEdx;
                        //     // append internal and lar values
                        //     if (exited_lar == false)
                        //     {
                        //         // if we are in the active volume
                        //         if (volume.volume_type == 2)
                        //         {
                        //             if (entered_lar == false)
                        //             {
                        //                 entered_lar = true;
                        //                 enter_lar[0] = particle.T(i);
                        //                 enter_lar[1] = particle.Vx(i);
                        //                 enter_lar[2] = particle.Vy(i);
                        //                 enter_lar[3] = particle.Vz(i);
                        //             }
                        //             // accumulate the internal dedx values
                        //             internal_distance += dx;
                        //             internal_dedx.emplace_back(dEdx);
                        //             // accumulate the internal LAr values
                        //             if (volume.material == 18)
                        //             {
                        //                 lar_distance += dx;
                        //                 lar_dedx.emplace_back(dEdx);
                        //             }
                        //         }
                        //         else
                        //         {
                        //             if (entered_lar == true)
                        //             {
                        //                 exited_lar = true;
                        //                 exit_lar[0] = particle.T(i);
                        //                 exit_lar[1] = particle.Vx(i);
                        //                 exit_lar[2] = particle.Vy(i);
                        //                 exit_lar[3] = particle.Vz(i);
                        //             }
                        //         }
                        //     }
                        // }
                        // // save summary statistics
                        // neutronTrajectories.enter_lar_t.emplace_back(enter_lar[0]);
                        // neutronTrajectories.enter_lar_x.emplace_back(enter_lar[1]);
                        // neutronTrajectories.enter_lar_y.emplace_back(enter_lar[2]);
                        // neutronTrajectories.enter_lar_z.emplace_back(enter_lar[3]);
                        // neutronTrajectories.exit_lar_t.emplace_back(exit_lar[0]);
                        // neutronTrajectories.exit_lar_x.emplace_back(exit_lar[1]);
                        // neutronTrajectories.exit_lar_y.emplace_back(exit_lar[2]);
                        // neutronTrajectories.exit_lar_z.emplace_back(exit_lar[3]);
                        // neutronTrajectories.total_distance.emplace_back(total_distance);
                        // neutronTrajectories.internal_distance.emplace_back(internal_distance);
                        // neutronTrajectories.lar_distance.emplace_back(lar_distance);
                        // // if the user wants dEdx, then save it too
                        // if (fSavePrimaryNeutrondEdx)
                        // {
                        //     neutronTrajectories.total_dedx.emplace_back(total_dedx);
                        //     neutronTrajectories.internal_dedx.emplace_back(internal_dedx);
                        //     neutronTrajectories.lar_dedx.emplace_back(lar_dedx);
                        // }
                        // // If the user wants the full trajectory, save it
                        // if (fSavePrimaryNeutronTrajectories)
                        // {
                        //     neutronTrajectories.t.emplace_back(t);
                        //     neutronTrajectories.x.emplace_back(x);
                        //     neutronTrajectories.y.emplace_back(y);
                        //     neutronTrajectories.z.emplace_back(z);
                        //     neutronTrajectories.energy.emplace_back(energy);
                        //     neutronTrajectories.px.emplace_back(px);
                        //     neutronTrajectories.py.emplace_back(py);
                        //     neutronTrajectories.pz.emplace_back(pz);
                        //     neutronTrajectories.volume_name.emplace_back(volume_name);
                        //     neutronTrajectories.material_name.emplace_back(material_name);
                        //     neutronTrajectories.material.emplace_back(material);
                        // }
                        // // otherwise, just save the beginning and end points
                        // else
                        // {
                        //     std::vector<Double_t> t = {particle.T(), particle.EndT()};
                        //     std::vector<Double_t> x = {particle.Vx(), particle.EndX()};
                        //     std::vector<Double_t> y = {particle.Vy(), particle.EndY()};
                        //     std::vector<Double_t> z = {particle.Vz(), particle.EndZ()};
                        //     std::vector<Double_t> energy = {particle.E(), particle.EndE()};
                        //     std::vector<Double_t> px = {particle.Px(), particle.EndPx()};
                        //     std::vector<Double_t> py = {particle.Py(), particle.EndPy()};
                        //     std::vector<Double_t> pz = {particle.Pz(), particle.EndPz()};
                        //     DetectorVolume begin_volume = fGeometry->getVolume(
                        //         particle.Vx(), particle.Vy(), particle.Vz()
                        //     );
                        //     DetectorVolume end_volume = fGeometry->getVolume(
                        //         particle.EndX(), particle.EndY(), particle.EndZ()
                        //     );
                        //     std::vector<std::string> volume_name = {
                        //         begin_volume.volume_name, end_volume.volume_name
                        //     };
                        //     std::vector<std::string> material_name = {
                        //         begin_volume.material_name, end_volume.material_name
                        //     };
                        //     std::vector<Double_t> material = {
                        //         begin_volume.material, end_volume.material
                        //     };
                        //     neutronTrajectories.t.emplace_back(t);
                        //     neutronTrajectories.x.emplace_back(x);
                        //     neutronTrajectories.y.emplace_back(y);
                        //     neutronTrajectories.z.emplace_back(z);
                        //     neutronTrajectories.energy.emplace_back(energy);
                        //     neutronTrajectories.px.emplace_back(px);
                        //     neutronTrajectories.py.emplace_back(py);
                        //     neutronTrajectories.pz.emplace_back(pz);
                        //     neutronTrajectories.volume_name.emplace_back(volume_name);
                        //     neutronTrajectories.material_name.emplace_back(material_name);
                        //     neutronTrajectories.material.emplace_back(material);
                        // }
                    }
                    // Now we look for actual captures.
                    if (particle.EndProcess() == "nCapture")
                    {
                        DetectorVolume endingVolume = fGeometry->getVolume(
                            particle.EndX(), particle.EndY(), particle.EndZ()
                        );
                        if (endingVolume.material_name == "LAr")
                        {
                            neutronList.neutron_ids.emplace_back(particle.TrackId());
                            neutronList.neutron_capture_x.emplace_back(particle.EndX());
                            neutronList.neutron_capture_y.emplace_back(particle.EndY());
                            neutronList.neutron_capture_z.emplace_back(particle.EndZ());
                        }
                    }
                }
                if (particle.PdgCode() == 13)
                {
                    muonList.primary_muons += 1;
                    muonList.muon_ids.emplace_back(particle.TrackId());
                }
                // check if the particle is a gamma
                if (particle.PdgCode() == 22)
                {
                    for(size_t i = 0; i < neutronList.neutron_ids.size(); i++)
                    {
                        if (neutronList.neutron_ids[i] == particle.Mother())
                        {
                            neutronList.gamma_ids.emplace_back(particle.TrackId());
                            neutronList.gamma_neutron_ids.emplace_back(particle.Mother());
                            neutronList.gamma_energy.emplace_back(particle.E());
                            neutronList.gamma_electron_energy.emplace_back(0);
                            neutronList.gamma_edep_energy.emplace_back(0);
                        }
                    }
                }
                // check if the particle is an electron
                if (particle.PdgCode() == 11)
                {
                    // check for gammas first
                    for (size_t i = 0; i < neutronList.gamma_ids.size(); i++)
                    {
                        if (neutronList.gamma_ids[i] == particle.Mother())
                        {
                            neutronList.electron_ids.emplace_back(particle.TrackId());
                            neutronList.electron_parent.emplace_back(particle.Mother());
                            neutronList.electron_gamma_ids.emplace_back(particle.Mother());
                            neutronList.gamma_electron_energy[i] += particle.E();
                            // find the corresponding neutron id
                            for(size_t j = 0; j < neutronList.neutron_ids.size(); j++)
                            {
                                if (neutronList.neutron_ids[j] == neutronList.gamma_neutron_ids[i])
                                {
                                    neutronList.electron_neutron_ids.emplace_back(neutronList.neutron_ids[j]);
                                }
                            }
                            neutronList.electron_energy.emplace_back(particle.E());
                        }
                    }
                    // then check electrons
                    for (size_t i = 0; i < neutronList.electron_ids.size(); i++)
                    {
                        if (neutronList.electron_ids[i] == particle.Mother())
                        {
                            neutronList.electron_ids.emplace_back(particle.TrackId());
                            neutronList.electron_parent.emplace_back(particle.Mother());
                            // find the corresponding gamma
                            neutronList.electron_gamma_ids.emplace_back(neutronList.electron_gamma_ids[i]);
                            neutronList.electron_neutron_ids.emplace_back(neutronList.electron_neutron_ids[i]);
                            neutronList.electron_energy.emplace_back(particle.E());
                        }
                    }
                }
            }
        }
        auto mcEnergyDeposit = event.getValidHandle<std::vector<sim::SimEnergyDeposit>>(fIonAndScintProducerLabel);
        if (mcEnergyDeposit.isValid())
        {
            for (auto energyDeposit : *mcEnergyDeposit)
            {
                // check the list of electrons
                for (size_t i = 0; i < neutronList.electron_ids.size(); i++)
                {
                    if (neutronList.electron_ids[i] == energyDeposit.TrackID())
                    {
                        neutronList.edep_parent.emplace_back(energyDeposit.TrackID());
                        neutronList.edep_neutron_ids.emplace_back(neutronList.electron_neutron_ids[i]);
                        neutronList.edep_gamma_ids.emplace_back(neutronList.electron_gamma_ids[i]);
                        for (size_t j = 0; j < neutronList.gamma_ids.size(); j++)
                        {
                            if (neutronList.gamma_ids[j] == neutronList.electron_gamma_ids[i])
                            {
                                neutronList.gamma_edep_energy[j] += energyDeposit.Energy();
                            }
                        }
                        neutronList.edep_energy.emplace_back(energyDeposit.Energy());
                        neutronList.edep_num_electrons.emplace_back(energyDeposit.NumElectrons());
                        neutronList.edep_x.emplace_back(energyDeposit.StartX());
                        neutronList.edep_y.emplace_back(energyDeposit.StartY());
                        neutronList.edep_z.emplace_back(energyDeposit.StartZ());
                        continue;
                    }
                }
                for (size_t i = 0; i < muonList.muon_ids.size(); i++)
                {
                    if (muonList.muon_ids[i] == energyDeposit.TrackID())
                    {
                        muonList.muon_edep_ids.emplace_back(energyDeposit.TrackID());
                        muonList.muon_edep_energy.emplace_back(energyDeposit.Energy());
                        muonList.muon_edep_num_electrons.emplace_back(energyDeposit.NumElectrons());
                        muonList.muon_edep_x.emplace_back(energyDeposit.StartX());
                        muonList.muon_edep_y.emplace_back(energyDeposit.StartY());
                        muonList.muon_edep_z.emplace_back(energyDeposit.StartZ());
                    }
                }
            }
        }
        // iterate over hits
        if (fFindHits == true)
        {
            std::vector<art::Ptr<recob::Hit>> allHits;
            auto hitHandle = event.getValidHandle<std::vector<recob::Hit>>(fHitFinderProducerLabel);
            if (hitHandle.isValid())
            {
                art::fill_ptr_vector(allHits, hitHandle);
                std::map<int,int> trueParticleHits, trueParticleHitsView0, trueParticleHitsView1, trueParticleHitsView2;
                for (const auto& hit : allHits)
                {
                    TruthMatchUtils::G4ID g4ID(TruthMatchUtils::TrueParticleID(clockData, hit, false));
                    if (TruthMatchUtils::Valid(g4ID))
                    {
                        neutronList.hit_track_ids.emplace_back(g4ID);
                        neutronList.hit_view.emplace_back(hit->View());
                    }
                    
                }
            }
            auto spacePointHandle = event.getValidHandle<std::vector<recob::SpacePoint>>(fSpacePointProducerLabel);
            if (spacePointHandle.isValid())
            {
                for (auto space_point : *spacePointHandle)
                {
                    const Double_t *pos = space_point.XYZ();
                    neutronList.space_point_ids.emplace_back(space_point.ID());
                    neutronList.space_point_x.emplace_back(pos[0]);
                    neutronList.space_point_y.emplace_back(pos[1]);
                    neutronList.space_point_z.emplace_back(pos[2]);
                }
            }
        }
        // save neutron trajectories
        fTempNeutronTrajectories.event_id = neutronTrajectories.event_id;
        fTempNeutronTrajectories.track_id = neutronTrajectories.track_id;
        fTempNeutronTrajectories.mother   = neutronTrajectories.mother;
        fTempNeutronTrajectories.process  = neutronTrajectories.process;
        fTempNeutronTrajectories.end_process = neutronTrajectories.end_process;
        fTempNeutronTrajectories.t = neutronTrajectories.t;
        fTempNeutronTrajectories.x = neutronTrajectories.x;
        fTempNeutronTrajectories.y = neutronTrajectories.y;
        fTempNeutronTrajectories.z = neutronTrajectories.z;
        fTempNeutronTrajectories.energy = neutronTrajectories.energy;
        fTempNeutronTrajectories.px = neutronTrajectories.px;
        fTempNeutronTrajectories.py = neutronTrajectories.py;
        fTempNeutronTrajectories.pz = neutronTrajectories.pz;
        fTempNeutronTrajectories.volume_name = neutronTrajectories.volume_name;
        fTempNeutronTrajectories.material_name = neutronTrajectories.material_name;
        fTempNeutronTrajectories.material = neutronTrajectories.material;
        fTempNeutronTrajectories.enter_lar_t = neutronTrajectories.enter_lar_t;
        fTempNeutronTrajectories.enter_lar_x = neutronTrajectories.enter_lar_x;
        fTempNeutronTrajectories.enter_lar_y = neutronTrajectories.enter_lar_y;
        fTempNeutronTrajectories.enter_lar_z = neutronTrajectories.enter_lar_z;
        fTempNeutronTrajectories.exit_lar_t = neutronTrajectories.exit_lar_t;
        fTempNeutronTrajectories.exit_lar_x = neutronTrajectories.exit_lar_x;
        fTempNeutronTrajectories.exit_lar_y = neutronTrajectories.exit_lar_y;
        fTempNeutronTrajectories.exit_lar_z = neutronTrajectories.exit_lar_z;
        fTempNeutronTrajectories.total_distance = neutronTrajectories.total_distance;
        fTempNeutronTrajectories.internal_distance = neutronTrajectories.internal_distance;
        fTempNeutronTrajectories.lar_distance = neutronTrajectories.lar_distance;
        fTempNeutronTrajectories.total_dedx = neutronTrajectories.total_dedx;
        fTempNeutronTrajectories.internal_dedx = neutronTrajectories.internal_dedx;
        fTempNeutronTrajectories.lar_dedx = neutronTrajectories.lar_dedx;

        fNeutronTrajectoryTree->Fill();
        if (neutronList.edep_x.size() > 0)
        {
            fTempNeutronList.event_id = neutronList.event_id;
            fTempNeutronList.primary_neutrons = neutronList.primary_neutrons;
            fTempNeutronList.neutron_ids = neutronList.neutron_ids;
            fTempNeutronList.neutron_capture_x = neutronList.neutron_capture_x;
            fTempNeutronList.neutron_capture_y = neutronList.neutron_capture_y;
            fTempNeutronList.neutron_capture_z = neutronList.neutron_capture_z;
            fTempNeutronList.gamma_ids = neutronList.gamma_ids;
            fTempNeutronList.gamma_neutron_ids = neutronList.gamma_neutron_ids;
            fTempNeutronList.gamma_energy = neutronList.gamma_energy;
            fTempNeutronList.gamma_electron_energy = neutronList.gamma_electron_energy;
            fTempNeutronList.gamma_edep_energy = neutronList.gamma_edep_energy;
            fTempNeutronList.electron_ids = neutronList.electron_ids;
            fTempNeutronList.electron_neutron_ids = neutronList.electron_neutron_ids;
            fTempNeutronList.electron_gamma_ids = neutronList.electron_gamma_ids;
            fTempNeutronList.electron_parent = neutronList.electron_parent;
            fTempNeutronList.electron_energy = neutronList.electron_energy;
            fTempNeutronList.edep_parent = neutronList.edep_parent;
            fTempNeutronList.edep_neutron_ids = neutronList.edep_neutron_ids;
            fTempNeutronList.edep_gamma_ids = neutronList.edep_gamma_ids;
            fTempNeutronList.edep_energy = neutronList.edep_energy;
            fTempNeutronList.edep_num_electrons = neutronList.edep_num_electrons;
            fTempNeutronList.edep_x = neutronList.edep_x;
            fTempNeutronList.edep_y = neutronList.edep_y;
            fTempNeutronList.edep_z = neutronList.edep_z;
            fTempNeutronList.hit_track_ids = neutronList.hit_track_ids;
            fTempNeutronList.hit_view = neutronList.hit_view;
            fTempNeutronList.space_point_ids = neutronList.space_point_ids;
            fTempNeutronList.space_point_x = neutronList.space_point_x;
            fTempNeutronList.space_point_y = neutronList.space_point_y;
            fTempNeutronList.space_point_z = neutronList.space_point_z;
            fNeutronTree->Fill();
        }
        if (muonList.muon_edep_x.size() > 0)
        {
            fTempMuonList.primary_muons = muonList.primary_muons;
            fTempMuonList.muon_ids = muonList.muon_ids;
            fTempMuonList.muon_edep_ids = muonList.muon_edep_ids;
            fTempMuonList.muon_edep_energy = muonList.muon_edep_energy;
            fTempMuonList.muon_edep_num_electrons = muonList.muon_edep_num_electrons;
            fTempMuonList.muon_edep_x = muonList.muon_edep_x;
            fTempMuonList.muon_edep_y = muonList.muon_edep_y;
            fTempMuonList.muon_edep_z = muonList.muon_edep_z;
            fMuonTree->Fill();
        }
        if (fGenerateNeutronCosmicVoxels)
        {
            fTempVoxels = fVoxelizer.generateLabeledNeutronCosmicVoxels(
                fNumberOfEvents-1,
                fVoxelSize,
                fTempNeutronList.edep_x,
                fTempNeutronList.edep_y,
                fTempNeutronList.edep_z,
                fTempNeutronList.edep_energy,
                fTempNeutronList.edep_neutron_ids,
                fTempMuonList.muon_edep_x,
                fTempMuonList.muon_edep_y,
                fTempMuonList.muon_edep_z,
                fTempMuonList.muon_edep_energy,
                fTempMuonList.muon_edep_ids,
                fDiscretizeVoxelFeatures,
                fUseMixedVoxelLabeling
            );
            fVoxelTree->Fill();
        }
    }
    // begin job
    void NeutronExtractor::beginJob()
    {
        fGeometry->FillTTree();
        fNumberOfEvents = 0;
        fNumberOfNeutronsPerEvent.clear();
    }
    // end job
    void NeutronExtractor::endJob()
    {
        // global neutron info
        fMetaTree->Branch("number_of_events", &fNumberOfEvents);
        fMetaTree->Branch("number_of_neutrons_per_event", &fNumberOfNeutronsPerEvent);
        fMetaTree->Fill();
    }
}

DEFINE_ART_MODULE(neutron::NeutronExtractor)
