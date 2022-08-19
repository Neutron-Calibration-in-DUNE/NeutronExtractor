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

        mSingleNeutronCaptureTTree->Branch("event", &mSingleNeutronCapture.event);
        mSingleNeutronCaptureTTree->Branch("track_id", &mSingleNeutronCapture.track_id);
        mSingleNeutronCaptureTTree->Branch("parent_track_id", &mSingleNeutronCapture.parent_track_id);
        mSingleNeutronCaptureTTree->Branch("parent_pdg", &mSingleNeutronCapture.parent_pdg);
        mSingleNeutronCaptureTTree->Branch("mc_creation_process", &mSingleNeutronCapture.mc_creation_process);
        mSingleNeutronCaptureTTree->Branch("mc_ending_process", &mSingleNeutronCapture.mc_ending_process);
        mSingleNeutronCaptureTTree->Branch("capture_tpc", &mSingleNeutronCapture.capture_tpc);
        mSingleNeutronCaptureTTree->Branch("capture_tpc_lar", &mSingleNeutronCapture.capture_tpc_lar);

        mSingleNeutronCaptureTTree->Branch("mc_t", &mSingleNeutronCapture.mc_t);
        mSingleNeutronCaptureTTree->Branch("mc_x", &mSingleNeutronCapture.mc_x);
        mSingleNeutronCaptureTTree->Branch("mc_y", &mSingleNeutronCapture.mc_y);
        mSingleNeutronCaptureTTree->Branch("mc_z", &mSingleNeutronCapture.mc_z);
        mSingleNeutronCaptureTTree->Branch("mc_energy", &mSingleNeutronCapture.mc_energy);
        mSingleNeutronCaptureTTree->Branch("mc_volume", &mSingleNeutronCapture.mc_volume);
        mSingleNeutronCaptureTTree->Branch("mc_material", &mSingleNeutronCapture.mc_material);
        
        mSingleNeutronCaptureTTree->Branch("mc_daughter_track_ids", &mSingleNeutronCapture.mc_daughter_track_ids);
        mSingleNeutronCaptureTTree->Branch("mc_daughter_pdgs", &mSingleNeutronCapture.mc_daughter_pdgs);

        mSingleNeutronCaptureTTree->Branch("edep_x", &mSingleNeutronCapture.edep_x);
        mSingleNeutronCaptureTTree->Branch("edep_y", &mSingleNeutronCapture.edep_y);
        mSingleNeutronCaptureTTree->Branch("edep_z", &mSingleNeutronCapture.edep_z);
        mSingleNeutronCaptureTTree->Branch("edep_energy", &mSingleNeutronCapture.edep_energy);
        mSingleNeutronCaptureTTree->Branch("edep_volume", &mSingleNeutronCapture.edep_volume);
        mSingleNeutronCaptureTTree->Branch("edep_material", &mSingleNeutronCapture.edep_material);
        mSingleNeutronCaptureTTree->Branch("edep_track_id", &mSingleNeutronCapture.edep_track_id);
        mSingleNeutronCaptureTTree->Branch("edep_pdg", &mSingleNeutronCapture.edep_pdg);
    }

    SingleNeutronCaptures::~SingleNeutronCaptures()
    {
    }

    void SingleNeutronCaptures::processEvent(
        Int_t event_id,
        ParticleMap particleMap,
        const art::ValidHandle<std::vector<simb::MCParticle>>& mcParticles,
        const art::ValidHandle<std::vector<sim::SimEnergyDeposit>>& mcEnergyDeposits
    )
    {
        if (mcParticles.isValid() && mcEnergyDeposits.isValid())
        {
            mSingleNeutronCaptureList.clear();
            for (auto particle : *mcParticles)
            {
                if (particle.PdgCode() == 2112)
                {
                    // set basic quantities
                    SingleNeutronCapture singleNeutronCapture(particle.NumberTrajectoryPoints());
                    singleNeutronCapture.event = event_id;
                    singleNeutronCapture.track_id = particle.TrackId();
                    singleNeutronCapture.parent_track_id = particle.Mother();
                    singleNeutronCapture.parent_pdg = particleMap.GetParticleParentPDG(particle.TrackId());
                    singleNeutronCapture.mc_creation_process = particle.Process();
                    singleNeutronCapture.mc_ending_process = particle.EndProcess();
                    
                    // get all trajectory info
                    for (size_t ii = 0; ii < particle.NumberTrajectoryPoints(); ii++)
                    {
                        singleNeutronCapture.mc_t[ii] = particle.T(ii);
                        singleNeutronCapture.mc_x[ii] = particle.Vx(ii);
                        singleNeutronCapture.mc_y[ii] = particle.Vy(ii);
                        singleNeutronCapture.mc_z[ii] = particle.Vz(ii);
                        singleNeutronCapture.mc_energy[ii] = particle.E(ii);

                        DetectorVolume volume = mGeometry->getVolume(
                            particle.Vx(ii), particle.Vy(ii), particle.Vz(ii)
                        );
                        singleNeutronCapture.mc_volume[ii] = volume.volume_name;
                        singleNeutronCapture.mc_material[ii] = volume.material_name;
                    }
                    DetectorVolume ending_volume = mGeometry->getVolume(
                        particle.EndX(), particle.EndY(), particle.EndZ()
                    );
                    if (ending_volume.volume_type == 2) {
                        singleNeutronCapture.capture_tpc = true;
                    }
                    else {
                        singleNeutronCapture.capture_tpc = false;
                    }
                    if (
                        ending_volume.volume_type == 2 && 
                        ending_volume.material_name == "LAr"
                    ) {
                        singleNeutronCapture.capture_tpc_lar = true;
                    }
                    else {
                        singleNeutronCapture.capture_tpc_lar = false;
                    }

                    for (Int_t jj = 0; jj < particle.NumberDaughters(); jj++)
                    {
                        singleNeutronCapture.mc_daughter_track_ids.emplace_back(particle.Daughter(jj));
                        singleNeutronCapture.mc_daughter_pdgs.emplace_back(particleMap.GetParticlePDG(particle.Daughter(jj)));
                    }

                    mSingleNeutronCaptureList.emplace_back(singleNeutronCapture);
                    mSingleNeutronCaptureMap[particle.TrackId()] = mSingleNeutronCaptureList.size() - 1;
                }
                else if (particle.PdgCode() == 22)
                {
                    
                }
                else if (particle.PdgCode() == 11)
                {

                }
            }
            for (auto energyDeposit : *mcEnergyDeposits)
            {
                if (particleMap.GetParticleAncestorPDG(energyDeposit.TrackID()) == 2112)
                {
                    Int_t neutronIndex = mSingleNeutronCaptureMap[
                        particleMap.GetParticleAncestorTrackID(energyDeposit.TrackID())
                    ];
                    mSingleNeutronCaptureList[neutronIndex].edep_x.emplace_back(energyDeposit.StartX());
                    mSingleNeutronCaptureList[neutronIndex].edep_y.emplace_back(energyDeposit.StartY());
                    mSingleNeutronCaptureList[neutronIndex].edep_z.emplace_back(energyDeposit.StartZ());
                    mSingleNeutronCaptureList[neutronIndex].edep_energy.emplace_back(energyDeposit.E());
                    DetectorVolume volume = mGeometry->getVolume(
                        energyDeposit.StartX(), energyDeposit.StartY(), energyDeposit.StartZ()
                    );
                    mSingleNeutronCaptureList[neutronIndex].edep_volume.emplace_back(volume.volume_name);
                    mSingleNeutronCaptureList[neutronIndex].edep_material.emplace_back(volume.material_name);
                    mSingleNeutronCaptureList[neutronIndex].edep_track_id.emplace_back(energyDeposit.TrackID());
                    mSingleNeutronCaptureList[neutronIndex].edep_pdg.emplace_back(
                        particleMap.GetParticlePDG(energyDeposit.TrackID())
                    );
                }
            }

            for (size_t ii = 0; ii < mSingleNeutronCaptureList.size(); ii++)
            {
                mSingleNeutronCapture = mSingleNeutronCaptureList[ii];
                mSingleNeutronCaptureTTree->Fill();
            }
        }
    }
}