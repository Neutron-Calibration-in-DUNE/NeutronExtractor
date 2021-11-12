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
#include "MCNeutron.h"
#include "MCGamma.h"
#include "MCElectron.h"

namespace neutron {
    class NeutronExtractor;
}

namespace neutron
{
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
        // event variables
        int fRun;
        int fSubRun;
        int fEvent;

        // number of events
        Int_t fNumberOfEvents;

        // mc neutrons
        MCNeutron fMCNeutrons;
        // number of neutrons per event
        std::vector<Int_t> fNumberOfNeutronsPerEvent;
        // list of neutrons for each event
        std::vector<std::vector<Int_t>> fListOfNeutrons;

        // mc gammas
        MCGamma fMCGammas;
        // list of gammas for each event
        std::vector<std::vector<Int_t>> fListOfGammas;

        // mc electrons
        MCElectron fMCElectrons;
        // list of electrons for each event
        std::vector<std::vector<Int_t>> fListOfElectrons;
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
        consumes<std::vector<simb::MCParticle>>(fLArGeantProducerLabel);
    }

    // check list of neutrons
    bool NeutronExtractor::checkListOfNeutrons(Int_t eventId, Int_t trackId)
    {
        if (eventId-1 >= fNumberOfEvents) { 
            return false; 
        }
        else {
            for (size_t i = 0; i < fListOfNeutrons[eventId-1].size(); i++) {
                if (fListOfNeutrons[eventId-1][i] == trackId) {
                    return true;
                }
            }
        }
        return false;
    }

    // check list of gammas
    bool NeutronExtractor::checkListOfGammas(Int_t eventId, Int_t trackId)
    {
        if (eventId-1 >= fNumberOfEvents) { 
            return false; 
        }
        else {
            for (size_t i = 0; i < fListOfGammas[eventId-1].size(); i++) {
                if (fListOfGammas[eventId-1][i] == trackId) {
                    return true;
                }
            }
        }
        return false;
    }

    // check list of electrons
    bool NeutronExtractor::checkListOfElectrons(Int_t eventId, Int_t trackId)
    {
        if (eventId-1 >= fNumberOfEvents) { 
            return false; 
        }
        else {
            for (size_t i = 0; i < fListOfElectrons[eventId-1].size(); i++) {
                if (fListOfElectrons[eventId-1][i] == trackId) {
                    return true;
                }
            }
        }
        return false;
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
        fListOfNeutrons.emplace_back(std::vector<Int_t>());
        fListOfGammas.emplace_back(std::vector<Int_t>());
        fListOfElectrons.emplace_back(std::vector<Int_t>());

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
                    fMCNeutrons.addNeutron(fEvent, particle);
                    fNumberOfNeutronsPerEvent[fEvent-1]++;
                    fListOfNeutrons[fEvent-1].emplace_back(particle.TrackId());
                    std::cout << "Neutron: " << particle.TrackId() << "," << particle.Mother() << std::endl;
                }
                // check if the particle is a gamma
                if (particle.PdgCode() == 22)
                {
                    // check that the gamma has a parent which is a neutron in the list
                    if (checkListOfNeutrons(fEvent, particle.Mother()))
                    {
                        // check if the gamma comes from a capture
                        if (particle.Process() == "nCapture") {
                            fMCNeutrons.addGamma(fEvent, particle);
                            fMCGammas.addGamma(fEvent, particle);
                            // add gamma to the list
                            fListOfGammas[fEvent-1].emplace_back(particle.TrackId());
                            std::cout << "Gamma: " << particle.TrackId() << "," << particle.Mother() << std::endl;
                        }
                    }
                }
                // check if the particle is an electron
                if (particle.PdgCode() == 11)
                {
                    // check that the electron has a parent which is a gamma 
                    // or another electron in the list
                    if (checkListOfGammas(fEvent, particle.Mother()) || 
                        checkListOfElectrons(fEvent, particle.Mother())
                    )
                    {
                        fMCElectrons.addElectron(fEvent, particle);
                        // add electron to the list
                        fListOfElectrons[fEvent-1].emplace_back(particle.TrackId());
                        DetectorVolume currentVolume = fGeometry->getVolume(
                            particle.Vx(), particle.Vy(), particle.Vz()
                        );
                        std::cout << "Electron: " << particle.TrackId() << "," << particle.Mother() <<  "," <<  currentVolume.material_name << "," << particle.E() << "," << particle.EndE() << "," << particle.T() << std::endl;
                    }
                }
            }
        }
        // get the energy depositions from IonAndScint for the event
        auto mcGeantEnergyDeposit = event.getValidHandle<std::vector<sim::SimEnergyDeposit>>(fLArGeantEnergyDepositProducerLabel);
        if (mcGeantEnergyDeposit.isValid())
        {

            for (auto GeantEnergyDeposit : *mcGeantEnergyDeposit)
            {
                DetectorVolume currentVolume = fGeometry->getVolume(
                    GeantEnergyDeposit.StartX(), GeantEnergyDeposit.StartY(), GeantEnergyDeposit.StartZ()
                );
                std::cout << "ED: " << GeantEnergyDeposit.TrackID() << "," << currentVolume.material_name << "," << GeantEnergyDeposit.Energy() << "," << GeantEnergyDeposit.StartT() << std::endl;
                // std::cout << fEvent << "," << GeantEnergyDeposit.TrackID();
                // std::cout << "," << GeantEnergyDeposit.NumElectrons();
                // std::cout << ",(" << GeantEnergyDeposit.StartX() << "," << GeantEnergyDeposit.StartY() << "," << GeantEnergyDeposit.StartZ();
                // std::cout << "),(" << GeantEnergyDeposit.T0() << "," << GeantEnergyDeposit.T1() << ")"  << std::endl;
                // for (auto particle : *mcParticles)
                // {   
                //     if (particle.TrackId() == GeantEnergyDeposit.TrackID())
                //     {
                //         std::cout << "\t" << particle.TrackId() << ',' << particle.PdgCode() << "," << particle.Mother();
                //         std::cout << ",(" << particle.Vx() << "," << particle.Vy() << "," << particle.Vz();
                //         std::cout << "),(" << particle.EndX() << "," << particle.EndY() << "," << particle.EndZ();
                //         std::cout << "),(" << particle.T() << "," << particle.EndT() << ")" << std::endl;
                //         for (auto parent : *mcParticles)
                //         {
                //             if (parent.TrackId() == particle.Mother())
                //             {
                //                 std::cout << "\t\t" << parent.TrackId() << ',' << parent.PdgCode() << "," << parent.Mother();
                //                 std::cout << ",(" << parent.Vx() << "," << parent.Vy() << "," << parent.Vz();
                //                 std::cout << "),(" << parent.EndX() << "," << parent.EndY() << "," << parent.EndZ();
                //                 std::cout << "),(" << parent.T() << "," << parent.EndT() << ")" << std::endl;
                //             }
                //         }
                //     }
                // }

                // check if the deposit has a parent in the list of electrons
                if (checkListOfElectrons(fEvent, GeantEnergyDeposit.TrackID()))
                {
                    std::cout << "Edep: " << GeantEnergyDeposit.TrackID() << std::endl;
                }
            }
        }
        auto mcEnergyDeposit = event.getValidHandle<std::vector<sim::SimEnergyDeposit>>(fIonAndScintProducerLabel);
        if (mcEnergyDeposit.isValid())
        {
            for (auto energyDeposit : *mcEnergyDeposit)
            {
                DetectorVolume currentVolume = fGeometry->getVolume(
                    energyDeposit.StartX(), energyDeposit.StartY(), energyDeposit.StartZ()
                );
                std::cout << "NEST: " << energyDeposit.TrackID() << "," << currentVolume.material_name << "," << energyDeposit.Energy() << "," << energyDeposit.StartT() << std::endl;
                //std::cout << fEvent << "," << energyDeposit.TrackID() << "," << energyDeposit.NumElectrons() << std::endl;

                // check if the deposit has a parent in the list of electrons
                if (checkListOfElectrons(fEvent, energyDeposit.TrackID()))
                {
                    std::cout << "Edep: " << energyDeposit.TrackID() << std::endl;
                }
            }
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
        fMCNeutrons.FillTTree();
        fMCGammas.FillTTree();
        fMCElectrons.FillTTree();
        // global neutron info
        fMetaTree->Branch("number_of_events", &fNumberOfEvents);
        fMetaTree->Branch("number_of_neutrons_per_event", &fNumberOfNeutronsPerEvent);
        fMetaTree->Fill();
    }
}

DEFINE_ART_MODULE(neutron::NeutronExtractor)
