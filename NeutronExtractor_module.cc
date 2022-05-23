/**
 * @file    NeutronExtractor_module.cc
 * @brief   A module for extracting truth/reco information about G4 particle trajectories
 *          and conducting some standard analysis tasks. 
 *          Generated at Mon Oct 11 11:21:12 2021 using cetskelgen
 * @ingroup NeutronExtractor
 * @author  Nicholas Carrara (nmcarrara@ucdavis.edu),
 *          Yashwanth Bezawada
**/
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
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"
#include "art_root_io/TFileService.h"
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

#include <TTree.h>
#include <TH1.h>
#include "TH1F.h"
#include "TGeoMaterial.h"
#include "TGeoElement.h"
#include <cmath>

#include "Configuration.h"
#include "DetectorGeometry.h"
#include "NeutronCapture.h"

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

        art::InputTag mLArGeantProducerLabel;
        art::InputTag mIonAndScintProducerLabel;
        art::InputTag mCluster3DProducerLabel;
        art::InputTag mSimChannelProducerLabel;
        art::InputTag mSimChannelInstanceProducerLabel;

        bool mFillNeutronCapture;

        /// ROOT output through art::TFileService
        /** We will save different TTrees to different TFiles specified 
         *  by the directories for each type.
         */ 
        art::ServiceHandle<art::TFileService> mTFileService;
        /// TTrees
        TTree *mMetaTree;

        // geometry information
        DetectorGeometry* mGeometry = DetectorGeometry::getInstance("NeutronExtractor");
        // MC neutron captures
        NeutronCapture mNeutronCapture;
    };

    // constructor
    NeutronExtractor::NeutronExtractor(const Parameters& config)
    : EDAnalyzer(config)
    , mParameters(config)
    {
        mLArGeantProducerLabel =    mParameters().LArGeantProducerLabel();
        mIonAndScintProducerLabel = mParameters().IonAndScintProducerLabel();
        mCluster3DProducerLabel =   mParameters().Cluster3DProducerLabel();
        mSimChannelProducerLabel =  mParameters().SimChannelProducerLabel();
        mSimChannelInstanceProducerLabel = mParameters().SimChannelInstanceProducerLabel();
        
        mMetaTree = mTFileService->make<TTree>("meta", "meta");
    }

    // begin job
    void NeutronExtractor::beginJob()
    {
        mGeometry->FillTTree();
    }

    // analyze function
    void NeutronExtractor::analyze(art::Event const& event)
    {
        /**
         * @details For each event, we will look through the various
         * available data products and send event info to the 
         * corresponding submodules that process them, starting with MCParticles
         * 
         */
        art::Handle<std::vector<simb::MCParticle>> particleHandle;
        if (!event.getByLabel(mLArGeantProducerLabel, particleHandle))
        {
            // if there are no particles for the event truth, then
            // we are in big trouble haha.  throw an exception
            throw cet::exception("NeutronExtractor")
                << " No simb::MCParticle objects in this event - "
                << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
        }
        // get the list of MC particles from Geant4
        std::cout << "Collecting MC Particles.." << std::endl;
        auto mcParticles = event.getValidHandle<std::vector<simb::MCParticle>>(mLArGeantProducerLabel);
        if (mFillNeutronCapture) 
        {
            auto const clockData(art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(event)); 
            auto mcSimChannels = 
                event.getValidHandle<std::vector<sim::SimChannel>>(
                    art::InputTag(mSimChannelProducerLabel.label(), mSimChannelInstanceProducerLabel.label())
                );
            auto recoSpacePoints = event.getValidHandle<std::vector<recob::SpacePoint>>(mCluster3DProducerLabel);
            art::FindManyP<recob::Hit> hitsFromSpsCluster3DAssn(recoSpacePoints, event, mCluster3DProducerLabel); 
            mNeutronCapture.processEvent(
                clockData,
                mcParticles, 
                mcSimChannels,
                recoSpacePoints,
                hitsFromSpsCluster3DAssn
            );
        }
    }
    
    // end job
    void NeutronExtractor::endJob()
    {
        // save configuration parameters
        mMetaTree->Branch("LArGeantProducerLabel",    &mLArGeantProducerLabel);
        mMetaTree->Branch("IonAndScintProducerLabel", &mIonAndScintProducerLabel);
        mMetaTree->Branch("Cluster3DProducerLabel",   &mCluster3DProducerLabel);
        mMetaTree->Branch("SimChannelProducerLabel",   &mSimChannelProducerLabel);
        mMetaTree->Branch("SimChannelInstanceProducerLabel", &mSimChannelInstanceProducerLabel);
 
        mMetaTree->Fill();
    }
}
DEFINE_ART_MODULE(neutron::NeutronExtractor)
