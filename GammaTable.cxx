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
        mGammaTree = mTFileService->make<TTree>("gammas", "gammas");
        mGammaTree->Branch("track_id", &mGamma.track_id);
        mGammaTree->Branch("neutron_id", &mGamma.neutron_id);
        mGammaTree->Branch("energy", &mGamma.energy);
        mGammaTree->Branch("start_x", &mGamma.start_x);
        mGammaTree->Branch("start_y", &mGamma.start_y);
        mGammaTree->Branch("start_z", &mGamma.start_z);
        mGammaTree->Branch("end_x", &mGamma.end_x);
        mGammaTree->Branch("end_y", &mGamma.end_y);
        mGammaTree->Branch("end_z", &mGamma.end_z);
        mGammaTree->Branch("gamma_x", &mGamma.gamma_x);
        mGammaTree->Branch("gamma_y", &mGamma.gamma_y);
        mGammaTree->Branch("gamma_z", &mGamma.gamma_z);
        mGammaTree->Branch("daughter_ids", &mGamma.daughter_ids);
        mGammaTree->Branch("daughter_level", &mGamma.daughter_level);
        mGammaTree->Branch("daughter_energy", &mGamma.daughter_energy);
        mGammaTree->Branch("daughter_start_x", &mGamma.daughter_start_x);
        mGammaTree->Branch("daughter_start_y", &mGamma.daughter_start_y);
        mGammaTree->Branch("daughter_start_z", &mGamma.daughter_start_z);
        mGammaTree->Branch("daughter_edep_energy", &mGamma.daughter_edep_energy);
        mGammaTree->Branch("daughter_edep_x", &mGamma.daughter_edep_x);
        mGammaTree->Branch("daughter_edep_y", &mGamma.daughter_edep_y);
        mGammaTree->Branch("daughter_edep_z", &mGamma.daughter_edep_z);
        mGammaTree->Branch("daughter_edep_num_electrons", &mGamma.daughter_edep_num_electrons);
        mGammaTree->Branch("daughter_edep_num_photons", &mGamma.daughter_edep_num_photons);

        mGammaTree->Branch("daughter_reco_sp_x", &mGamma.daughter_reco_sp_x);
        mGammaTree->Branch("daughter_reco_sp_y", &mGamma.daughter_reco_sp_y);
        mGammaTree->Branch("daughter_reco_sp_z", &mGamma.daughter_reco_sp_z);
        mGammaTree->Branch("daughter_reco_sp_x_sigma", &mGamma.daughter_reco_sp_x_sigma);
        mGammaTree->Branch("daughter_reco_sp_y_sigma", &mGamma.daughter_reco_sp_y_sigma);
        mGammaTree->Branch("daughter_reco_sp_z_sigma", &mGamma.daughter_reco_sp_z_sigma);
        mGammaTree->Branch("daughter_reco_sp_chisq", &mGamma.daughter_reco_sp_chisq);

        mGammaTree->Branch("daughter_reco_peak_time", &mGamma.daughter_reco_peak_time);
        mGammaTree->Branch("daughter_reco_peak_time_sigma", &mGamma.daughter_reco_peak_time_sigma);
        mGammaTree->Branch("daughter_reco_rms", &mGamma.daughter_reco_rms);
        mGammaTree->Branch("daughter_reco_peak_amplitude", &mGamma.daughter_reco_peak_amplitude);
        mGammaTree->Branch("daughter_reco_peak_amplitude_sigma", &mGamma.daughter_reco_peak_amplitude_sigma);
        mGammaTree->Branch("daughter_reco_summed_adc", &mGamma.daughter_reco_summed_adc);
        mGammaTree->Branch("daughter_reco_view", &mGamma.daughter_reco_view);

        mGammaStatisticsTree = mTFileService->make<TTree>("gamma_statistics", "gamma_statistics");
        mGammaStatisticsTree->Branch("total_num_gammas", &mGammaStatistics.total_num_gammas);
        mGammaStatisticsTree->Branch("energy", &mGammaStatistics.energy);
        mGammaStatisticsTree->Branch("num_gammas_mc", &mGammaStatistics.num_gammas_mc);
        mGammaStatisticsTree->Branch("num_gammas_reco", &mGammaStatistics.num_gammas_reco);
        mGammaStatisticsTree->Branch("num_mc_points", &mGammaStatistics.num_mc_points);
        mGammaStatisticsTree->Branch("num_reco_points", &mGammaStatistics.num_reco_points);
    }

    GammaTable::~GammaTable()
    {}

    void GammaTable::processEvent(
        detinfo::DetectorClocksData const& clockData,
        const art::ValidHandle<std::vector<simb::MCParticle>>& mcParticles,
        const art::ValidHandle<std::vector<sim::SimEnergyDeposit>>& mcEnergyDeposits,
        const art::ValidHandle<std::vector<sim::SimChannel>>& mcChannels,
        const art::ValidHandle<std::vector<recob::SpacePoint>>& recoSpacePoints,
        const art::FindManyP<recob::Hit>& hitPandoraSPsAssn 
    )
    {
        if (mcParticles.isValid() and mcChannels.isValid() and recoSpacePoints.isValid())
        {
            std::vector<Gamma> gammas;
            std::map<int, int> gamma_map;
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
                            gammas.emplace_back(
                                Gamma(
                                    particle.TrackId(), particle.Mother(), round(particle.E()*10e6)/10e6, 
                                    particle.Vx(), particle.Vy(), particle.Vz(),
                                    particle.EndX(), particle.EndY(), particle.EndZ()
                                )
                            );
                            for(size_t k = 0; k < particle.NumberTrajectoryPoints(); k++)
                            {
                                gammas.back().gamma_x.emplace_back(particle.Vx(k));
                                gammas.back().gamma_y.emplace_back(particle.Vy(k));
                                gammas.back().gamma_z.emplace_back(particle.Vz(k));
                            }
                            gamma_map[particle.TrackId()] = gammas.size()-1;
                        }
                    }
                }
                else if (particle.PdgCode() == 11)
                {
                    for(size_t i = 0; i < gammas.size(); i++)
                    {
                        if (gammas[i].track_id == particle.Mother())
                        {
                            gammas[i].daughter_ids.emplace_back(particle.TrackId());
                            gammas[i].daughter_level.emplace_back(0);
                            gammas[i].daughter_energy.emplace_back(round(particle.E()*10e6)/10e6);
                            gammas[i].daughter_start_x.emplace_back(particle.Vx());
                            gammas[i].daughter_start_y.emplace_back(particle.Vy());
                            gammas[i].daughter_start_z.emplace_back(particle.Vz());

                            gammas[i].daughter_edep_energy.emplace_back(std::vector<Double_t>());
                            gammas[i].daughter_edep_x.emplace_back(std::vector<Double_t>());
                            gammas[i].daughter_edep_y.emplace_back(std::vector<Double_t>());
                            gammas[i].daughter_edep_z.emplace_back(std::vector<Double_t>());
                            gammas[i].daughter_edep_num_electrons.emplace_back(std::vector<Int_t>());
                            gammas[i].daughter_edep_num_photons.emplace_back(std::vector<Int_t>());
                            
                            gammas[i].daughter_reco_sp_x.emplace_back(std::vector<Double_t>());
                            gammas[i].daughter_reco_sp_y.emplace_back(std::vector<Double_t>());
                            gammas[i].daughter_reco_sp_z.emplace_back(std::vector<Double_t>());
                            gammas[i].daughter_reco_sp_x_sigma.emplace_back(std::vector<Double_t>());
                            gammas[i].daughter_reco_sp_y_sigma.emplace_back(std::vector<Double_t>());
                            gammas[i].daughter_reco_sp_z_sigma.emplace_back(std::vector<Double_t>());
                            gammas[i].daughter_reco_sp_chisq.emplace_back(std::vector<Double_t>());

                            gammas[i].daughter_reco_peak_time.emplace_back(std::vector<Double_t>());
                            gammas[i].daughter_reco_peak_time_sigma.emplace_back(std::vector<Double_t>());
                            gammas[i].daughter_reco_rms.emplace_back(std::vector<Double_t>());
                            gammas[i].daughter_reco_peak_amplitude.emplace_back(std::vector<Double_t>());
                            gammas[i].daughter_reco_peak_amplitude_sigma.emplace_back(std::vector<Double_t>());
                            gammas[i].daughter_reco_summed_adc.emplace_back(std::vector<Double_t>());
                            gammas[i].daughter_reco_view.emplace_back(std::vector<Int_t>());
                            gamma_map[particle.TrackId()] = i;
                        }
                        for (size_t j = 0; j < gammas[i].daughter_ids.size(); j++)
                        {
                            if (gammas[i].daughter_ids[j] == particle.Mother())
                            {
                                gammas[i].daughter_ids.emplace_back(particle.TrackId());
                                gammas[i].daughter_level.emplace_back(gammas[i].daughter_level[j]+1);
                                gammas[i].daughter_energy.emplace_back(round(particle.E()*10e6)/10e6);
                                gammas[i].daughter_start_x.emplace_back(particle.Vx());
                                gammas[i].daughter_start_y.emplace_back(particle.Vy());
                                gammas[i].daughter_start_z.emplace_back(particle.Vz());

                                gammas[i].daughter_edep_energy.emplace_back(std::vector<Double_t>());
                                gammas[i].daughter_edep_x.emplace_back(std::vector<Double_t>());
                                gammas[i].daughter_edep_y.emplace_back(std::vector<Double_t>());
                                gammas[i].daughter_edep_z.emplace_back(std::vector<Double_t>());
                                gammas[i].daughter_edep_num_electrons.emplace_back(std::vector<Int_t>());
                                gammas[i].daughter_edep_num_photons.emplace_back(std::vector<Int_t>());
                                
                                gammas[i].daughter_reco_sp_x.emplace_back(std::vector<Double_t>());
                                gammas[i].daughter_reco_sp_y.emplace_back(std::vector<Double_t>());
                                gammas[i].daughter_reco_sp_z.emplace_back(std::vector<Double_t>());
                                gammas[i].daughter_reco_sp_x_sigma.emplace_back(std::vector<Double_t>());
                                gammas[i].daughter_reco_sp_y_sigma.emplace_back(std::vector<Double_t>());
                                gammas[i].daughter_reco_sp_z_sigma.emplace_back(std::vector<Double_t>());
                                gammas[i].daughter_reco_sp_chisq.emplace_back(std::vector<Double_t>());

                                gammas[i].daughter_reco_peak_time.emplace_back(std::vector<Double_t>());
                                gammas[i].daughter_reco_peak_time_sigma.emplace_back(std::vector<Double_t>());
                                gammas[i].daughter_reco_rms.emplace_back(std::vector<Double_t>());
                                gammas[i].daughter_reco_peak_amplitude.emplace_back(std::vector<Double_t>());
                                gammas[i].daughter_reco_peak_amplitude_sigma.emplace_back(std::vector<Double_t>());
                                gammas[i].daughter_reco_summed_adc.emplace_back(std::vector<Double_t>());
                                gammas[i].daughter_reco_view.emplace_back(std::vector<Int_t>());
                                gamma_map[particle.TrackId()] = i;
                            }
                        }
                    }
                }
            }
            
            for (auto energyDeposit : *mcEnergyDeposits)
            {
                if (energyDeposit.PdgCode() == 11)
                {
                    if (gamma_map.find(energyDeposit.TrackID()) != gamma_map.end())
                    {
                        auto gamma_index = gamma_map[energyDeposit.TrackID()];
                        for (size_t j = 0; j < gammas[gamma_index].daughter_ids.size(); j++)
                        {
                            if (energyDeposit.TrackID() == gammas[gamma_index].daughter_ids[j])
                            {
                                gammas[gamma_index].daughter_edep_energy[j].emplace_back(energyDeposit.Energy());
                                gammas[gamma_index].daughter_edep_x[j].emplace_back(energyDeposit.StartX());
                                gammas[gamma_index].daughter_edep_y[j].emplace_back(energyDeposit.StartY());
                                gammas[gamma_index].daughter_edep_z[j].emplace_back(energyDeposit.StartZ());
                                gammas[gamma_index].daughter_edep_num_electrons[j].emplace_back(energyDeposit.NumElectrons());
                                gammas[gamma_index].daughter_edep_num_photons[j].emplace_back(energyDeposit.NumPhotons());
                                gammas[gamma_index].num_edep_points += 1;
                            }
                        }
                    }
                }
            }
            std::vector<art::Ptr<recob::SpacePoint>> pointsList;
            art::fill_ptr_vector(pointsList, recoSpacePoints);            
            for (size_t i = 0; i < pointsList.size(); i++)
            {
                auto& spsHit = hitPandoraSPsAssn.at(i);
                bool space_point = false;
                Int_t track_id;
                for (auto hit : spsHit)
                {  
                    track_id = TruthMatchUtils::TrueParticleID(
                        clockData, hit, false
                    );
                    if (gamma_map.find(track_id) != gamma_map.end())
                    {
                        auto gamma_index = gamma_map[track_id];
                        for (size_t j = 0; j < gammas[gamma_index].daughter_ids.size(); j++)
                        {
                            if (track_id == gammas[gamma_index].daughter_ids[j])
                            {
                                space_point = true;
                                // collect results
                                auto xyz = pointsList[i]->XYZ();
                                auto xyz_sigma = pointsList[i]->ErrXYZ();
                                // check if point is in active volume
                                // Determine if edep is within the desired volume
                                gammas[gamma_index].daughter_reco_sp_x[j].emplace_back(xyz[0]);
                                gammas[gamma_index].daughter_reco_sp_y[j].emplace_back(xyz[1]);
                                gammas[gamma_index].daughter_reco_sp_z[j].emplace_back(xyz[2]);
                                gammas[gamma_index].daughter_reco_sp_x_sigma[j].emplace_back(xyz_sigma[0]);
                                gammas[gamma_index].daughter_reco_sp_y_sigma[j].emplace_back(xyz_sigma[4]);
                                gammas[gamma_index].daughter_reco_sp_z_sigma[j].emplace_back(xyz_sigma[8]);
                                gammas[gamma_index].daughter_reco_sp_chisq[j].emplace_back(pointsList[i]->Chisq());

                                gammas[gamma_index].daughter_reco_peak_time[j].emplace_back(hit->PeakTime());
                                gammas[gamma_index].daughter_reco_peak_time_sigma[j].emplace_back(hit->SigmaPeakTime());
                                gammas[gamma_index].daughter_reco_rms[j].emplace_back(hit->RMS());
                                gammas[gamma_index].daughter_reco_peak_amplitude[j].emplace_back(hit->PeakAmplitude());
                                gammas[gamma_index].daughter_reco_peak_amplitude_sigma[j].emplace_back(hit->SigmaPeakAmplitude());
                                gammas[gamma_index].daughter_reco_summed_adc[j].emplace_back(hit->SummedADC());
                                gammas[gamma_index].daughter_reco_view[j].emplace_back(hit->View());
                            }
                        }
                    }
                }
                if (space_point) {
                    gammas[gamma_map[track_id]].num_reco_points += 1;
                }
            }
            for (size_t i = 0; i < gammas.size(); i++)
            {
                mGamma = gammas[i];
                mGammaTree->Fill();
            }
            GammaStatistics gamma_statistics;
            for(size_t i = 0; i < gammas.size(); i++)
            {
                bool energy_exists = false;
                for (size_t j = 0; j < gamma_statistics.energy.size(); j++)
                {
                    if (gamma_statistics.energy[j] == gammas[i].energy)
                    {
                        if (gammas[i].num_edep_points > 0) {
                            gamma_statistics.num_gammas_mc[j] += 1;
                            gamma_statistics.num_mc_points[j].emplace_back(gammas[i].num_edep_points);
                        }
                        if (gammas[i].num_reco_points > 0) {
                            gamma_statistics.num_gammas_reco[j] += 1;
                            gamma_statistics.num_reco_points[j].emplace_back(gammas[i].num_reco_points);
                        }
                        energy_exists = true;
                    }
                }
                if (!energy_exists)
                {
                    gamma_statistics.energy.emplace_back(gammas[i].energy);
                    if (gammas[i].num_edep_points > 0) {
                        gamma_statistics.num_gammas_mc.emplace_back(1);
                        gamma_statistics.num_mc_points.emplace_back(std::vector<Int_t>({gammas[i].num_edep_points}));
                    }
                    else {
                        gamma_statistics.num_gammas_mc.emplace_back(0);
                        gamma_statistics.num_mc_points.emplace_back(std::vector<Int_t>({gammas[i].num_edep_points}));
                    }
                    if (gammas[i].num_reco_points > 0) {
                        gamma_statistics.num_gammas_reco.emplace_back(1);
                        gamma_statistics.num_reco_points.emplace_back(std::vector<Int_t>({gammas[i].num_reco_points}));
                    }
                    else {
                        gamma_statistics.num_gammas_reco.emplace_back(0);
                        gamma_statistics.num_reco_points.emplace_back(std::vector<Int_t>({gammas[i].num_reco_points}));
                    }
                }
            }
            gamma_statistics.total_num_gammas = gammas.size();
            mGammaStatistics = gamma_statistics;
            mGammaStatisticsTree->Fill();
        }
    }
}