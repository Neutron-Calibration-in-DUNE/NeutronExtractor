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
            fhicl::Atom<art::InputTag> SimulationLabel
            {
                fhicl::Name("SimulationLabel"),
                fhicl::Comment("tag of the input data product with the detector simulation")
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

    private:
        art::InputTag fSimulationProducerLabel;
        art::InputTag fOutputFileArt;
        // geometry information
        DetectorGeometry fGeometry;
        // ROOT
        art::ServiceHandle<art::TFileService> fTFileService;
        TTree *fMetaTree;
        // event variables
        int fRun;
        int fSubRun;
        int fEvent;

        // mc neutrons
        MCNeutron fMCNeutrons;

    };

    // constructor
    NeutronExtractor::NeutronExtractor(Parameters const& config)
    : EDAnalyzer(config)
    , fSimulationProducerLabel(config().SimulationLabel())
    , fOutputFileArt(config().OutputFile())
    {
        fMetaTree = fTFileService->make<TTree>("meta", "meta");
        consumes<std::vector<simb::MCParticle>>(fSimulationProducerLabel);
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
        if (!event.getByLabel(fSimulationProducerLabel, particleHandle))
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

        // tell MCNeutrons that you have a new event
        fMCNeutrons.initializeNewEvent();

        // get the list of MC particles from Geant4
        auto mcParticles = event.getValidHandle<std::vector<simb::MCParticle>>(fSimulationProducerLabel);
        if (mcParticles.isValid())
        {
            for (auto particle : *mcParticles)
            {
                // check if the particle is a neutron
                if (particle.PdgCode() == 2112)
                {
                    fMCNeutrons.addNeutron(fEvent, particle);
                }
            }
        }
    }
    // begin job
    void NeutronExtractor::beginJob()
    {
        fGeometry.FillTTree();
    }
    // end job
    void NeutronExtractor::endJob()
    {
        fMCNeutrons.FillTTree();
        // global neutron info
        fMetaTree->Branch("number_of_events", &fMCNeutrons.getNumberOfEvents());
        fMetaTree->Branch("number_of_neutrons_per_event", &fMCNeutrons.getNumberOfNeutronsPerEvent());
        fMetaTree->Fill();
    }
}

DEFINE_ART_MODULE(neutron::NeutronExtractor)
