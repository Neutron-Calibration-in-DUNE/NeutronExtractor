/**
 * @file SingleNeutronCaptures.cxx
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2022-08-18
 */
#include "SingleNeutronCaptures.h"

namespace neutron
{
    SingleNeutronCaptures::SingleNeutronCaptures()
    {
        mSingleNeutronCaptureTTree = mTFileService->make<TTree>("single_neutron_capture", "single_neutron_capture");

        mSingleNeutronCaptureTTree->Branch("track_id", &mSingleNeutronCapture.track_id);
        mSingleNeutronCaptureTTree->Branch("parent_track_id", &mSingleNeutronCapture.parent_track_id);
        mSingleNeutronCaptureTTree->Branch("parent_pdg", &mSingleNeutronCapture.parent_pdg);

        mSingleNeutronCaptureTTree->Branch("mc_t", &mSingleNeutronCapture.mc_t);
        mSingleNeutronCaptureTTree->Branch("mc_x", &mSingleNeutronCapture.mc_x);
        mSingleNeutronCaptureTTree->Branch("mc_y", &mSingleNeutronCapture.mc_y);
        mSingleNeutronCaptureTTree->Branch("mc_z", &mSingleNeutronCapture.mc_z);
        mSingleNeutronCaptureTTree->Branch("mc_energy", &mSingleNeutronCapture.mc_energy);
        mSingleNeutronCaptureTTree->Branch("mc_process", &mSingleNeutronCapture.mc_process);
    }

    SingleNeutronCaptures::~SingleNeutronCaptures()
    {
    }

    void SingleNeutronCaptures::ResetSingleNeutronCapture(Int_t number_trajectory_points)
    {
        mSingleNeutronCapture.mc_t.resize(number_trajectory_points);
        mSingleNeutronCapture.mc_x.resize(number_trajectory_points);
        mSingleNeutronCapture.mc_y.resize(number_trajectory_points);
        mSingleNeutronCapture.mc_z.resize(number_trajectory_points);
        mSingleNeutronCapture.mc_energy.resize(number_trajectory_points);
        mSingleNeutronCapture.mc_process.resize(number_trajectory_points);

        mSingleNeutronCapture.gamma_indices.clear();
    }

    void SingleNeutronCaptures::processEvent(
        ParticleMap particleMap,
        const art::ValidHandle<std::vector<simb::MCParticle>>& mcParticles
    )
    {
        if (mcParticles.isValid())
        {
            for (auto particle : *mcParticles)
            {
                if (particle.PdgCode() == 2112)
                {
                    // set basic quantities
                    ResetSingleNeutronCapture(particle.NumberTrajectoryPoints());
                    mSingleNeutronCapture.track_id = particle.TrackId();
                    mSingleNeutronCapture.parent_track_id = particle.Mother();
                    mSingleNeutronCapture.parent_pdg = particleMap.GetParticleParentPDG(particle.TrackId());
                
                    mSingleNeutronCaptureTTree->Fill();
                }
            }
        }
    }
}