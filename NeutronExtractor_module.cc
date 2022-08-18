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
    }

    // analyze function
    void NeutronExtractor::analyze(art::Event const& event)
    {
    }
    
    // end job
    void NeutronExtractor::endJob()
    {        
    }
}
DEFINE_ART_MODULE(neutron::NeutronExtractor)
