/**
 * @file GammaTable.cxx
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2022-05-26
 */
#include "GammaTable.h"

namespace neutron
{
    GammaTable::GammaTable()
    {
    }

    GammaTable::~GammaTable()
    {}

    void GammaTable::processEvent(
        detinfo::DetectorClocksData const& clockData,
        const art::ValidHandle<std::vector<simb::MCParticle>>& mcParticles,
        const art::ValidHandle<std::vector<sim::SimEnergyDeposit>>& mcEnergyDeposits,
        const art::ValidHandle<std::vector<sim::SimChannel>>& mcChannels,
        const art::ValidHandle<std::vector<recob::SpacePoint>>& recoSpacePoints,
        const art::FindManyP<recob::Hit>& hitPandoraSPsAssn //to associate space points from pandora to hits
    )
    {
        if (mcParticles.isValid() and mcChannels.isValid() and recoSpacePoints.isValid())
        {
            mGammas.clear();
            /**
             * We first iterate through all particles and create a map of 
             * parent-daughter pairs for track ids.  This way we can search
             * recursively for ancestors of each particle.
             */
            std::map<Int_t, Int_t> parentDaughterMap;
            std::map<Int_t, Int_t> particlePDGMap;

            std::vector<int> neutron_captures;

            for (auto particle : *mcParticles)
            {
                parentDaughterMap[particle.TrackId()] = particle.Mother();
                particlePDGMap[particle.TrackId()] = particle.PdgCode();
                // check if the particle is a neutron
                if (particle.PdgCode() == 2112)
                {
                    if (particle.EndProcess() == "nCapture")
                    {
                        DetectorVolume ending_volume = fGeometry->getVolume(
                            particle.EndX(), particle.EndY(), particle.EndZ()
                        );
                        if (ending_volume.material_name == "LAr" and ending_volume.volume_type == 2) {
                            neutron_captures.emplace_back(particle.TrackId());
                        }
                    }
                }
                else if (particle.PdgCode() == 22)
                {
                    for(size_t i = 0; i < neutron_captures.size(); i++)
                    {
                        if (neutron_captures[i] == particle.Mother())
                        {
                            mGammas.emplace_back(
                                Gamma(
                                    particle.TrackId(), particle.Mother(), particle.E(), 
                                    particle.Vx(), particle.Vy(), particle.Vz(),
                                    particle.EndX(), particle.EndY(), particle.EndZ()
                                )
                            );
                            mGammaMap[particle.TrackId()] = mGammas.size()-1;
                        }
                    }
                }
                else if (particle.PdgCode() == 11)
                {
                    for(size_t i = 0; i < mGammas.size(); i++)
                    {
                        if (mGammas[i].track_id == particle.Mother())
                        {
                            mGammas[i].daughter_ids.emplace_back(particle.TrackId());
                            mGammas[i].daughter_energy.emplace_back(particle.E());
                            mGammas[i].daughter_start_x.emplace_back(particle.Vx());
                            mGammas[i].daughter_start_y.emplace_back(particle.Vy());
                            mGammas[i].daughter_start_z.emplace_back(particle.Vz());
                            mGammas[i].daughter_edep_energy.emplace_back(0);
                            mGammas[i].daughter_edep_x.emplace_back(-999);
                            mGammas[i].daughter_edep_y.emplace_back(-999);
                            mGammas[i].daughter_edep_z.emplace_back(-999);
                            mGammas[i].daughter_edep_num_electrons.emplace_back(0);
                            mGammas[i].daughter_edep_num_photons.emplace_back(0);
                            
                            mGammas[i].daughter_reco_sp_x.emplace_back(-999);
                            mGammas[i].daughter_reco_sp_y.emplace_back(-999);
                            mGammas[i].daughter_reco_sp_z.emplace_back(-999);
                            mGammas[i].daughter_reco_peak_time.emplace_back(0);
                            mGammas[i].daughter_reco_peak_time_sigma.emplace_back(0);
                            mGammas[i].daughter_reco_rms.emplace_back(0);
                            mGammas[i].daughter_reco_peak_amplitude.emplace_back(0);
                            mGammas[i].daughter_reco_peak_amplitude_sigma.emplace_back(0);
                            mGammas[i].daughter_reco_summed_adc.emplace_back(0);
                            mGammaMap[particle.TrackId()] = i;
                        }
                        for (size_t j = 0; j < mGammas[i].daughter_ids.size(); j++)
                        {
                            if (mGammas[i].daughter_ids[j] == particle.Mother())
                            {
                                mGammas[i].daughter_ids.emplace_back(particle.TrackId());
                                mGammas[i].daughter_energy.emplace_back(particle.E());
                                mGammas[i].daughter_start_x.emplace_back(particle.Vx());
                                mGammas[i].daughter_start_y.emplace_back(particle.Vy());
                                mGammas[i].daughter_start_z.emplace_back(particle.Vz());
                                mGammas[i].daughter_edep_energy.emplace_back(0);
                                mGammas[i].daughter_edep_x.emplace_back(-999);
                                mGammas[i].daughter_edep_y.emplace_back(-999);
                                mGammas[i].daughter_edep_z.emplace_back(-999);
                                mGammas[i].daughter_edep_num_electrons.emplace_back(0);
                                mGammas[i].daughter_edep_num_photons.emplace_back(0);
                                
                                mGammas[i].daughter_reco_sp_x.emplace_back(-999);
                                mGammas[i].daughter_reco_sp_y.emplace_back(-999);
                                mGammas[i].daughter_reco_sp_z.emplace_back(-999);
                                mGammas[i].daughter_reco_peak_time.emplace_back(0);
                                mGammas[i].daughter_reco_peak_time_sigma.emplace_back(0);
                                mGammas[i].daughter_reco_rms.emplace_back(0);
                                mGammas[i].daughter_reco_peak_amplitude.emplace_back(0);
                                mGammas[i].daughter_reco_peak_amplitude_sigma.emplace_back(0);
                                mGammas[i].daughter_reco_summed_adc.emplace_back(0);
                                mGammaMap[particle.TrackId()] = i;
                            }
                        }
                    }
                }
            }
            
            for (auto energyDeposit : *mcEnergyDeposits)
            {
                if (energyDeposit.PdgCode() == 11)
                {
                    auto particle_exists = std::find(
                        mGammaMap.begin(), 
                        mGammaMap.end(), 
                        energyDeposit.TrackID()
                    );
                    if (particle_exists != mGammaMap.end())
                    {
                        auto gamma_index = mGammaMap[energyDeposit.TrackID()];
                        for (size_t j = 0; j < mGammas[gamma_index].daughter_ids.size(); j++)
                        {
                            if (energyDeposit.TrackID() == mGammas[gamma_index].daughter_ids[j])
                            {
                                mGammas[gamma_index].daughter_edep_energy[j] = energyDeposit.Energy();
                                mGammas[gamma_index].daughter_edep_x[j] = energyDeposit.StartX();
                                mGammas[gamma_index].daughter_edep_y[j] = energyDeposit.StartY();
                                mGammas[gamma_index].daughter_edep_z[j] = energyDeposit.StartZ();
                                mGammas[gamma_index].daughter_edep_num_electrons[j] = energyDeposit.NumElectrons();
                                mGammas[gamma_index].daughter_edep_num_photons[j] = energyDeposit.NumPhotons();
                            }
                        }
                    }
                }
            }
            for (size_t i = 0; i < mGammas.size(); i++)
            {
                std::cout << "\ngamma: " << i << "\n";
                std::cout << "\tneutron_id: " << mGammas[i].neutron_id << "\n";
                std::cout << "\ttrack_id: " << mGammas[i].track_id << "\n";
                std::cout << "\tenergy: " << mGammas[i].energy << "\n";
                std::cout << "\tstart_x: " << mGammas[i].start_x << "\n";
                std::cout << "\tstart_y: " << mGammas[i].start_y << "\n";
                std::cout << "\tstart_z: " << mGammas[i].start_z << "\n";
                std::cout << "\tend_x: " << mGammas[i].end_x << "\n";
                std::cout << "\tend_y: " << mGammas[i].end_y << "\n";
                std::cout << "\tend_z: " << mGammas[i].end_z << "\n";
                for (size_t j = 0; j < mGammas[i].daughter_ids.size(); j++)
                {
                    std::cout << "\t\tdaughter_id: " << mGammas[i].daughter_ids[j] << "\n";
                    std::cout << "\t\tdaughter_energy: " << mGammas[i].daughter_energy[j] << "\n";
                    std::cout << "\t\tdaughter_start_x: " << mGammas[i].daughter_start_x[j] << "\n";
                    std::cout << "\t\tdaughter_start_y: " << mGammas[i].daughter_start_y[j] << "\n";
                    std::cout << "\t\tdaughter_start_z: " << mGammas[i].daughter_start_z[j] << "\n";
                    std::cout << "\t\tdaughter_edep_energy: " << mGammas[i].daughter_edep_energy[j] << "\n";
                    std::cout << "\t\tdaughter_edep_x: " << mGammas[i].daughter_edep_x[j] << "\n";
                    std::cout << "\t\tdaughter_edep_y: " << mGammas[i].daughter_edep_y[j] << "\n";
                    std::cout << "\t\tdaughter_edep_z: " << mGammas[i].daughter_edep_z[j] << "\n";
                    std::cout << "\t\tdaughter_edep_num_electrons: " << mGammas[i].daughter_edep_num_electrons[j] << "\n";
                    std::cout << "\t\tdaughter_edep_num_photons: " << mGammas[i].daughter_edep_num_photons[j] << "\n";
                }
                std::cout << "\tdaughter_total_energy: " << std::accumulate(mGammas[i].daughter_energy.begin(), mGammas[i].daughter_energy.end(), 0) << "\n";
            }
        }
    }
}