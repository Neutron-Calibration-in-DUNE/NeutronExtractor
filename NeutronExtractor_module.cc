/**
 * @file    NeutronExtractor_module.cc
 * @brief   A module for extracting truth/reco information about G4 particle trajectories
 *          and conducting some standard analysis tasks. 
 *          Generated at Mon Oct 11 11:21:12 2021 using cetskelgen
 * @ingroup NeutronExtractor
 * @author  Nicholas Carrara (nmcarrara@ucdavis.edu),
 *          Yashwanth Bezawada
**/
#include "Core.h"

#include "Configuration.h"
#include "DetectorGeometry.h"
#include "ParticleMap.h"
#include "SingleNeutronCaptures.h"

namespace neutron
{
    class NeutronExtractor : public art::EDAnalyzer
    {
    public:
        explicit NeutronExtractor(const Parameters& config);
        NeutronExtractor(const NeutronExtractor&) = delete;
        NeutronExtractor(NeutronExtractor&&) = delete;
        NeutronExtractor& operator=(const NeutronExtractor&) = delete;
        NeutronExtractor& operator=(NeutronExtractor&&) = delete;

        // required EDAnalyzer functions
        void analyze(const art::Event& event) override;
        void beginJob() override;
        void endJob() override;

    private:
        Parameters mParameters;
        DetectorGeometry* mDetectorGeometry = DetectorGeometry::getInstance("NeutronExtractor");
        ParticleMap mParticleMap;
        SingleNeutronCaptures mSingleNeutronCaptures;
    };

    // constructor
    NeutronExtractor::NeutronExtractor(const Parameters& config)
    : EDAnalyzer(config)
    , mParameters(config)
    {
    }

    // begin job
    void NeutronExtractor::beginJob()
    {
        mDetectorGeometry->FillTTree();
    }

    // analyze function
    void NeutronExtractor::analyze(art::Event const& event)
    {
        Int_t event_id = event.id().event();
        /**
         * @details For each event, we will look through the various
         * available data products and send event info to the 
         * corresponding submodules that process them, starting with MCParticles
         * 
         */
        art::Handle<std::vector<simb::MCParticle>> particleHandle;
        if (!event.getByLabel(mParameters().LArGeantProducerLabel(), particleHandle))
        {
            // if there are no particles for the event truth, then
            // we are in big trouble haha.  throw an exception
            throw cet::exception("NeutronExtractor")
                << " No simb::MCParticle objects in this event - "
                << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
        }

        // get the list of MC particles from Geant4
        std::cout << "Collecting MC Particles..." << std::endl;
        auto mcParticles = event.getValidHandle<std::vector<simb::MCParticle>>(mParameters().LArGeantProducerLabel());

        // generate particle map
        std::cout << "Generating Particle Map..." << std::endl;
        mParticleMap.processEvent(mcParticles);

        if (mParameters().FillSingleNeutronCaptures()) 
        {
            std::cout << "Filling Single Neutron Captures..." << std::endl;
            auto mcEnergyDeposit = event.getValidHandle<std::vector<sim::SimEnergyDeposit>>(mParameters().IonAndScintProducerLabel());
            mSingleNeutronCaptures.processEvent(event_id, mParticleMap, mcParticles, mcEnergyDeposit);
        }
    }
    
    // end job
    void NeutronExtractor::endJob()
    {        
    }
}
DEFINE_ART_MODULE(neutron::NeutronExtractor)
