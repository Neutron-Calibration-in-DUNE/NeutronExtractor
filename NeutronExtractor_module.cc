/**
 * @file    NeutronExtractor_module.cc
 * @brief   A module for extracting truth/reco information about G4 particle trajectories
 *          and conducting some standard analysis tasks. 
 *          Generated at Mon Oct 11 11:21:12 2021 using cetskelgen
 * @ingroup NeutronExtractor
 * @author  Nicholas Carrara (nmcarrara@ucdavis.edu),
 *          Jingbo Wang
 *          Junying Haung
 *          Yashwanth Bezawada
 *          Walker Johnson
**/

// art framework includes
#include "art/Framework/Core/EDAnalyzer.h"
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
#include "larsim/Utils/TruthMatchUtils.h"
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

namespace neutron {
    class NeutronExtractor;
}

namespace neutron
{
    struct EventList
    {
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

        EventList(Int_t event) : event_id(event){}
    };

    class NeutronExtractor : public art::EDAnalyzer
    {
    public:
        struct Config
        {
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
            fhicl::Atom<art::InputTag> HitFinderProducerLabel
            {
                fhicl::Name("HitFinderProducerLabel"),
                fhicl::Comment("tag of the data product which contains the hits")
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
        art::InputTag fLArGeantProducerLabel;
        art::InputTag fLArGeantEnergyDepositProducerLabel;
        art::InputTag fIonAndScintProducerLabel;
        bool fFindHits;
        art::InputTag fHitFinderProducerLabel;
        art::InputTag fOutputFileArt;
        // geometry information
        DetectorGeometry* fGeometry = DetectorGeometry::getInstance("NeutronExtractor");
        // ROOT
        art::ServiceHandle<art::TFileService> fTFileService;
        TTree *fMetaTree;
        TTree *fNeutronTree;
        // event variables
        int fRun;
        int fSubRun;
        int fEvent;

        // number of events
        Int_t fNumberOfEvents;

        EventList fTempEventList;

        std::vector<Int_t> fNumberOfNeutronsPerEvent;
    };

    // constructor
    NeutronExtractor::NeutronExtractor(Parameters const& config)
    : EDAnalyzer(config)
    , fLArGeantProducerLabel(config().LArGeantProducerLabel())
    , fLArGeantEnergyDepositProducerLabel(config().LArGeantEnergyDepositProducerLabel())
    , fIonAndScintProducerLabel(config().IonAndScintProducerLabel())
    , fFindHits(config().FindHits())
    , fHitFinderProducerLabel(config().HitFinderProducerLabel())
    , fOutputFileArt(config().OutputFile())
    , fTempEventList(0)
    {
        fMetaTree = fTFileService->make<TTree>("meta", "meta");
        fNeutronTree = fTFileService->make<TTree>("neutron", "neutron");
        consumes<std::vector<simb::MCParticle>>(fLArGeantProducerLabel);
        consumes<std::vector<sim::SimEnergyDeposit>>(fLArGeantEnergyDepositProducerLabel);

        fNeutronTree->Branch("event_id", &fTempEventList.event_id);
        fNeutronTree->Branch("primary_neutrons", &fTempEventList.primary_neutrons);
        fNeutronTree->Branch("neutron_ids", &fTempEventList.neutron_ids);
        fNeutronTree->Branch("neutron_capture_x", &fTempEventList.neutron_capture_x);
        fNeutronTree->Branch("neutron_capture_y", &fTempEventList.neutron_capture_y);
        fNeutronTree->Branch("neutron_capture_z", &fTempEventList.neutron_capture_z);
        fNeutronTree->Branch("gamma_ids", &fTempEventList.gamma_ids);
        fNeutronTree->Branch("gamma_neutron_ids", &fTempEventList.gamma_neutron_ids);
        fNeutronTree->Branch("gamma_energy", &fTempEventList.gamma_energy);
        fNeutronTree->Branch("gamma_electron_energy", &fTempEventList.gamma_electron_energy);
        fNeutronTree->Branch("gamma_edep_energy", &fTempEventList.gamma_edep_energy);
        fNeutronTree->Branch("electron_ids", &fTempEventList.electron_ids);
        fNeutronTree->Branch("electron_neutron_ids", &fTempEventList.electron_neutron_ids);
        fNeutronTree->Branch("electron_gamma_ids", &fTempEventList.electron_gamma_ids);
        fNeutronTree->Branch("electron_parent", &fTempEventList.electron_parent);
        fNeutronTree->Branch("electron_energy", &fTempEventList.electron_energy);
        fNeutronTree->Branch("edep_parent", &fTempEventList.edep_parent);
        fNeutronTree->Branch("edep_neutron_ids", &fTempEventList.edep_neutron_ids);
        fNeutronTree->Branch("edep_gamma_ids", &fTempEventList.edep_gamma_ids);
        fNeutronTree->Branch("edep_energy", &fTempEventList.edep_energy);
        fNeutronTree->Branch("edep_num_electrons", &fTempEventList.edep_num_electrons);
        fNeutronTree->Branch("edep_x", &fTempEventList.edep_x);
        fNeutronTree->Branch("edep_y", &fTempEventList.edep_y);
        fNeutronTree->Branch("edep_z", &fTempEventList.edep_z);
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
        EventList eventList(fNumberOfEvents-1);
        eventList.primary_neutrons = 0;
        // get clock information
        detinfo::DetectorClocksData const detClocks
            = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob();
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
                        eventList.primary_neutrons += 1;
                    }
                    DetectorVolume endingVolume = fGeometry->getVolume(
                        particle.EndX(), particle.EndY(), particle.EndZ()
                    );
                    if (particle.EndProcess() == "nCapture" and endingVolume.material_name == "LAr")
                    {
                        eventList.neutron_ids.emplace_back(particle.TrackId());
                        eventList.neutron_capture_x.emplace_back(particle.EndX());
                        eventList.neutron_capture_y.emplace_back(particle.EndY());
                        eventList.neutron_capture_z.emplace_back(particle.EndZ());
                    }

                }
                // check if the particle is a gamma
                if (particle.PdgCode() == 22)
                {
                    for(size_t i = 0; i < eventList.neutron_ids.size(); i++)
                    {
                        if (eventList.neutron_ids[i] == particle.Mother())
                        {
                            eventList.gamma_ids.emplace_back(particle.TrackId());
                            eventList.gamma_neutron_ids.emplace_back(particle.Mother());
                            eventList.gamma_energy.emplace_back(particle.E());
                            eventList.gamma_electron_energy.emplace_back(0);
                            eventList.gamma_edep_energy.emplace_back(0);
                        }
                    }
                }
                // check if the particle is an electron
                if (particle.PdgCode() == 11)
                {
                    // check for gammas first
                    for (size_t i = 0; i < eventList.gamma_ids.size(); i++)
                    {
                        if (eventList.gamma_ids[i] == particle.Mother())
                        {
                            eventList.electron_ids.emplace_back(particle.TrackId());
                            eventList.electron_parent.emplace_back(particle.Mother());
                            eventList.electron_gamma_ids.emplace_back(particle.Mother());
                            eventList.gamma_electron_energy[i] += particle.E();
                            // find the corresponding neutron id
                            for(size_t j = 0; j < eventList.neutron_ids.size(); j++)
                            {
                                if (eventList.neutron_ids[j] == eventList.gamma_neutron_ids[i])
                                {
                                    eventList.electron_neutron_ids.emplace_back(eventList.neutron_ids[j]);
                                }
                            }
                            eventList.electron_energy.emplace_back(particle.E());
                        }
                    }
                    // then check electrons
                    for (size_t i = 0; i < eventList.electron_ids.size(); i++)
                    {
                        if (eventList.electron_ids[i] == particle.Mother())
                        {
                            eventList.electron_ids.emplace_back(particle.TrackId());
                            eventList.electron_parent.emplace_back(particle.Mother());
                            // find the corresponding gamma
                            eventList.electron_gamma_ids.emplace_back(eventList.electron_gamma_ids[i]);
                            eventList.electron_neutron_ids.emplace_back(eventList.electron_neutron_ids[i]);
                            eventList.electron_energy.emplace_back(particle.E());
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
                for (size_t i = 0; i < eventList.electron_ids.size(); i++)
                {
                    if (eventList.electron_ids[i] == energyDeposit.TrackID())
                    {
                        eventList.edep_parent.emplace_back(energyDeposit.TrackID());
                        eventList.edep_neutron_ids.emplace_back(eventList.electron_neutron_ids[i]);
                        eventList.edep_gamma_ids.emplace_back(eventList.electron_gamma_ids[i]);
                        for (size_t j = 0; j < eventList.gamma_ids.size(); j++)
                        {
                            if (eventList.gamma_ids[j] == eventList.electron_gamma_ids[i])
                            {
                                eventList.gamma_edep_energy[j] += energyDeposit.Energy();
                            }
                        }
                        eventList.edep_energy.emplace_back(energyDeposit.Energy());
                        eventList.edep_num_electrons.emplace_back(energyDeposit.NumElectrons());
                        eventList.edep_x.emplace_back(energyDeposit.StartX());
                        eventList.edep_y.emplace_back(energyDeposit.StartY());
                        eventList.edep_z.emplace_back(energyDeposit.StartZ());
                    }
                }
            }
        }
        // iterate over hits
        if (fFindHits == true)
        {
            auto recoHits = event.getValidHandle<std::vector<recob::Hit>>(fHitFinderProducerLabel);
            if (recoHits.isValid())
            {
                for (auto hit : *recoHits)
                {
                    art::Ptr<recob::Hit> hit_ptr = static_cast<art::Ptr<recob::Hit>>(hit);
                    int trackId = TruthMatchUtils::TrueParticleID(detClocks, hit_ptr, false);
                    std::cout << "event: " << fEvent << ", hit track id: " << trackId << std::endl;
                }
            }
        }
        if (eventList.edep_x.size() > 0)
        {
            fTempEventList.event_id = eventList.event_id;
            fTempEventList.primary_neutrons = eventList.primary_neutrons;
            fTempEventList.neutron_ids = eventList.neutron_ids;
            fTempEventList.neutron_capture_x = eventList.neutron_capture_x;
            fTempEventList.neutron_capture_y = eventList.neutron_capture_y;
            fTempEventList.neutron_capture_z = eventList.neutron_capture_z;
            fTempEventList.gamma_ids = eventList.gamma_ids;
            fTempEventList.gamma_neutron_ids = eventList.gamma_neutron_ids;
            fTempEventList.gamma_energy = eventList.gamma_energy;
            fTempEventList.gamma_electron_energy = eventList.gamma_electron_energy;
            fTempEventList.gamma_edep_energy = eventList.gamma_edep_energy;
            fTempEventList.electron_ids = eventList.electron_ids;
            fTempEventList.electron_neutron_ids = eventList.electron_neutron_ids;
            fTempEventList.electron_gamma_ids = eventList.electron_gamma_ids;
            fTempEventList.electron_parent = eventList.electron_parent;
            fTempEventList.electron_energy = eventList.electron_energy;
            fTempEventList.edep_parent = eventList.edep_parent;
            fTempEventList.edep_neutron_ids = eventList.edep_neutron_ids;
            fTempEventList.edep_gamma_ids = eventList.edep_gamma_ids;
            fTempEventList.edep_energy = eventList.edep_energy;
            fTempEventList.edep_num_electrons = eventList.edep_num_electrons;
            fTempEventList.edep_x = eventList.edep_x;
            fTempEventList.edep_y = eventList.edep_y;
            fTempEventList.edep_z = eventList.edep_z;
            fNeutronTree->Fill();
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
