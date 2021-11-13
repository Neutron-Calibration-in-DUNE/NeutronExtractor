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
        std::vector<Int_t> neutron_ids;
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

        std::vector<EventList> fEventList;

        std::vector<Int_t> fNumberOfNeutronsPerEvent;
    };

    // constructor
    NeutronExtractor::NeutronExtractor(Parameters const& config)
    : EDAnalyzer(config)
    , fLArGeantProducerLabel(config().LArGeantProducerLabel())
    , fLArGeantEnergyDepositProducerLabel(config().LArGeantEnergyDepositProducerLabel())
    , fIonAndScintProducerLabel(config().IonAndScintProducerLabel())
    , fOutputFileArt(config().OutputFile())
    {
        fMetaTree = fTFileService->make<TTree>("meta", "meta");
        fNeutronTree = fTFileService->make<TTree>("neutron", "neutron");
        consumes<std::vector<simb::MCParticle>>(fLArGeantProducerLabel);
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
                    DetectorVolume currentVolume = fGeometry->getVolume(
                        particle.Vx(), particle.Vy(), particle.Vz()
                    );
                    if (particle.EndProcess() == "nCapture" and currentVolume.material_name == "LAr")
                    {
                        eventList.neutron_ids.emplace_back(particle.TrackId());
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
                                    eventList.electron_neutron_ids.emplace_back(eventList.neutrons_ids[j]);
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
                        for (size_t j = 0; j < eventList.gamma_ids; j++)
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
        fEventList.emplace_back(eventList);
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

        EventList event_list(0);
        fNeutronTree->Branch("event_id", &event_list.event_id);
        fNeutronTree->Branch("neutron_ids", &event_list.neutron_ids);
        fNeutronTree->Branch("gamma_ids", &event_list.gamma_ids);
        fNeutronTree->Branch("gamma_neutron_ids", &event_list.gamma_neutron_ids);
        fNeutronTree->Branch("gamma_energy", &event_list.gamma_energy);
        fNeutronTree->Branch("gamma_electron_energy", &event_list.gamma_electron_energy);
        fNeutronTree->Branch("gamma_edep_energy", &event_list.gamma_edep_energy);
        fNeutronTree->Branch("electron_ids", &event_list.electron_ids);
        fNeutronTree->Branch("electron_neutron_ids", &event_list.electron_neutron_ids);
        fNeutronTree->Branch("electron_gamma_ids", &event_list.electron_gamma_ids);
        fNeutronTree->Branch("electron_parent", &event_list.electron_parent);
        fNeutronTree->Branch("electron_energy", &event_list.electron_energy);
        fNeutronTree->Branch("edep_parent", &event_list.edep_parent);
        fNeutronTree->Branch("edep_neutron_ids", &event_list.edep_neutron_ids);
        fNeutronTree->Branch("edep_gamma_ids", &event_list.edep_gamma_ids);
        fNeutronTree->Branch("edep_energy", &event_list.edep_energy);
        fNeutronTree->Branch("edep_num_electrons", &event_list.edep_num_electrons);
        fNeutronTree->Branch("edep_x", &event_list.edep_x);
        fNeutronTree->Branch("edep_y", &event_list.edep_y);
        fNeutronTree->Branch("edep_z", &event_list.edep_z);
        for (size_t i = 0; i < fEventList.size(); i++) 
        {
            if (fEventList[i].edep_x.size() > 0)
            {
                event_list.event_id = fEventList[i].event_id;
                event_list.neutron_ids = fEventList[i].neutron_ids;
                event_list.gamma_ids = fEventList[i].gamma_ids;
                event_list.gamma_neutron_ids = fEventList[i].gamma_neutron_ids;
                event_list.gamma_energy = fEventList[i].gamma_energy;
                event_list.gamma_electron_energy = fEventList[i].gamma_electron_energy;
                event_list.gamma_edep_energy = fEventList[i].gamma_edep_energy;
                event_list.electron_ids = fEventList[i].electron_ids;
                event_list.electron_neutron_ids = fEventList[i].electron_neutron_ids;
                event_list.electron_gamma_ids = fEventList[i].electron_gamma_ids;
                event_list.electron_parent = fEventList[i].electron_parent;
                event_list.electron_energy = fEventList[i].electron_energy;
                event_list.edep_parent = fEventList[i].edep_parent;
                event_list.edep_neutron_ids = fEventList[i].edep_neutron_ids;
                event_list.edep_gamma_ids = fEventList[i].edep_gamma_ids;
                event_list.edep_energy = fEventList[i].edep_energy;
                event_list.edep_num_electrons = fEventList[i].edep_num_electrons;
                event_list.edep_x = fEventList[i].edep_x;
                event_list.edep_y = fEventList[i].edep_y;
                event_list.edep_z = fEventList[i].edep_z;
                fNeutronTree->Fill();
            }
        }

        fMetaTree->Fill();
    }
}

DEFINE_ART_MODULE(neutron::NeutronExtractor)
