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
        mSingleNeutronCaptureTTree->Branch("mc_creation_process", &mSingleNeutronCapture.mc_creation_process);
        mSingleNeutronCaptureTTree->Branch("mc_ending_process", &mSingleNeutronCapture.mc_ending_process);

        mSingleNeutronCaptureTTree->Branch("mc_t", &mSingleNeutronCapture.mc_t);
        mSingleNeutronCaptureTTree->Branch("mc_x", &mSingleNeutronCapture.mc_x);
        mSingleNeutronCaptureTTree->Branch("mc_y", &mSingleNeutronCapture.mc_y);
        mSingleNeutronCaptureTTree->Branch("mc_z", &mSingleNeutronCapture.mc_z);
        mSingleNeutronCaptureTTree->Branch("mc_energy", &mSingleNeutronCapture.mc_energy);
        
        mSingleNeutronCaptureTTree->Branch("mc_daughter_track_ids", &mSingleNeutronCapture.mc_daughter_track_ids);
        mSingleNeutronCaptureTTree->Branch("mc_daughter_pdgs", &mSingleNeutronCapture.mc_daughter_pdgs);
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

        mSingleNeutronCapture.mc_daughter_track_ids.clear();
        mSingleNeutronCapture.mc_daughter_pdgs.clear();
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
                    mSingleNeutronCapture.mc_creation_process = particle.Process();
                    mSingleNeutronCapture.mc_ending_process = particle.EndProcess();
                    
                    // get all trajectory info
                    for (size_t ii = 0; ii < particle.NumberTrajectoryPoints(); ii++)
                    {
                        mSingleNeutronCapture.mc_t[ii] = particle.Vt(ii);
                        mSingleNeutronCapture.mc_x[ii] = particle.Vx(ii);
                        mSingleNeutronCapture.mc_y[ii] = particle.Vy(ii);
                        mSingleNeutronCapture.mc_z[ii] = particle.Vz(ii);
                        mSingleNeutronCapture.mc_energy[ii] = particle.E(ii);
                    }

                    for (Int_t jj = 0; jj < particle.NumberDaughters(); jj++)
                    {
                        mSingleNeutronCapture.mc_daughter_track_ids.emplace_back(particle.Daughter(jj));
                        mSingleNeutronCapture.mc_daughter_pdgs.emplace_back(particleMap.GetParticlePDG(particle.Daughter(jj)));
                    }

                    mSingleNeutronCaptureTTree->Fill();
                }
            }
        }
    }
}