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

        mGammaStatisticsTree = mTFileService->make<TTree>("gamma_statistics", "gamma_statistics");
        mGammaStatisticsTree->Branch("total_num_gammas", &mGammaStatistics.total_num_gammas);
        mGammaStatisticsTree->Branch("energy", &mGammaStatistics.energy);
        mGammaStatisticsTree->Branch("num_gammas_mc", &mGammaStatistics.num_gammas_mc);
        mGammaStatisticsTree->Branch("num_gammas_reco", &mGammaStatistics.num_gammas_reco);
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
                                    particle.TrackId(), particle.Mother(), round(particle.E()*10e6)/10e6, 
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
                            mGammas[i].daughter_level.emplace_back(0);
                            mGammas[i].daughter_energy.emplace_back(round(particle.E()*10e6)/10e6);
                            mGammas[i].daughter_start_x.emplace_back(particle.Vx());
                            mGammas[i].daughter_start_y.emplace_back(particle.Vy());
                            mGammas[i].daughter_start_z.emplace_back(particle.Vz());

                            mGammas[i].daughter_edep_energy.emplace_back(std::vector<Double_t>());
                            mGammas[i].daughter_edep_x.emplace_back(std::vector<Double_t>());
                            mGammas[i].daughter_edep_y.emplace_back(std::vector<Double_t>());
                            mGammas[i].daughter_edep_z.emplace_back(std::vector<Double_t>());
                            mGammas[i].daughter_edep_num_electrons.emplace_back(std::vector<Int_t>());
                            mGammas[i].daughter_edep_num_photons.emplace_back(std::vector<Int_t>());
                            
                            mGammas[i].daughter_reco_sp_x.emplace_back(std::vector<Double_t>());
                            mGammas[i].daughter_reco_sp_y.emplace_back(std::vector<Double_t>());
                            mGammas[i].daughter_reco_sp_z.emplace_back(std::vector<Double_t>());
                            mGammas[i].daughter_reco_sp_x_sigma.emplace_back(std::vector<Double_t>());
                            mGammas[i].daughter_reco_sp_y_sigma.emplace_back(std::vector<Double_t>());
                            mGammas[i].daughter_reco_sp_z_sigma.emplace_back(std::vector<Double_t>());
                            mGammas[i].daughter_reco_sp_chisq.emplace_back(std::vector<Double_t>());

                            mGammas[i].daughter_reco_peak_time.emplace_back(std::vector<Double_t>());
                            mGammas[i].daughter_reco_peak_time_sigma.emplace_back(std::vector<Double_t>());
                            mGammas[i].daughter_reco_rms.emplace_back(std::vector<Double_t>());
                            mGammas[i].daughter_reco_peak_amplitude.emplace_back(std::vector<Double_t>());
                            mGammas[i].daughter_reco_peak_amplitude_sigma.emplace_back(std::vector<Double_t>());
                            mGammas[i].daughter_reco_summed_adc.emplace_back(std::vector<Double_t>());
                            mGammaMap[particle.TrackId()] = i;
                        }
                        for (size_t j = 0; j < mGammas[i].daughter_ids.size(); j++)
                        {
                            if (mGammas[i].daughter_ids[j] == particle.Mother())
                            {
                                mGammas[i].daughter_ids.emplace_back(particle.TrackId());
                                mGammas[i].daughter_level.emplace_back(mGammas[i].daughter_level[j]+1);
                                mGammas[i].daughter_energy.emplace_back(round(particle.E()*10e6)/10e6);
                                mGammas[i].daughter_start_x.emplace_back(particle.Vx());
                                mGammas[i].daughter_start_y.emplace_back(particle.Vy());
                                mGammas[i].daughter_start_z.emplace_back(particle.Vz());

                                mGammas[i].daughter_edep_energy.emplace_back(std::vector<Double_t>());
                                mGammas[i].daughter_edep_x.emplace_back(std::vector<Double_t>());
                                mGammas[i].daughter_edep_y.emplace_back(std::vector<Double_t>());
                                mGammas[i].daughter_edep_z.emplace_back(std::vector<Double_t>());
                                mGammas[i].daughter_edep_num_electrons.emplace_back(std::vector<Int_t>());
                                mGammas[i].daughter_edep_num_photons.emplace_back(std::vector<Int_t>());
                                
                                mGammas[i].daughter_reco_sp_x.emplace_back(std::vector<Double_t>());
                                mGammas[i].daughter_reco_sp_y.emplace_back(std::vector<Double_t>());
                                mGammas[i].daughter_reco_sp_z.emplace_back(std::vector<Double_t>());
                                mGammas[i].daughter_reco_sp_x_sigma.emplace_back(std::vector<Double_t>());
                                mGammas[i].daughter_reco_sp_y_sigma.emplace_back(std::vector<Double_t>());
                                mGammas[i].daughter_reco_sp_z_sigma.emplace_back(std::vector<Double_t>());
                                mGammas[i].daughter_reco_sp_chisq.emplace_back(std::vector<Double_t>());

                                mGammas[i].daughter_reco_peak_time.emplace_back(std::vector<Double_t>());
                                mGammas[i].daughter_reco_peak_time_sigma.emplace_back(std::vector<Double_t>());
                                mGammas[i].daughter_reco_rms.emplace_back(std::vector<Double_t>());
                                mGammas[i].daughter_reco_peak_amplitude.emplace_back(std::vector<Double_t>());
                                mGammas[i].daughter_reco_peak_amplitude_sigma.emplace_back(std::vector<Double_t>());
                                mGammas[i].daughter_reco_summed_adc.emplace_back(std::vector<Double_t>());
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
                    if (mGammaMap.find(energyDeposit.TrackID()) != mGammaMap.end())
                    {
                        auto gamma_index = mGammaMap[energyDeposit.TrackID()];
                        for (size_t j = 0; j < mGammas[gamma_index].daughter_ids.size(); j++)
                        {
                            if (energyDeposit.TrackID() == mGammas[gamma_index].daughter_ids[j])
                            {
                                mGammas[gamma_index].daughter_edep_energy[j].emplace_back(energyDeposit.Energy());
                                mGammas[gamma_index].daughter_edep_x[j].emplace_back(energyDeposit.StartX());
                                mGammas[gamma_index].daughter_edep_y[j].emplace_back(energyDeposit.StartY());
                                mGammas[gamma_index].daughter_edep_z[j].emplace_back(energyDeposit.StartZ());
                                mGammas[gamma_index].daughter_edep_num_electrons[j].emplace_back(energyDeposit.NumElectrons());
                                mGammas[gamma_index].daughter_edep_num_photons[j].emplace_back(energyDeposit.NumPhotons());
                                mGammas[gamma_index].num_edep_points += 1;
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
                    if (mGammaMap.find(track_id) != mGammaMap.end())
                    {
                        auto gamma_index = mGammaMap[track_id];
                        for (size_t j = 0; j < mGammas[gamma_index].daughter_ids.size(); j++)
                        {
                            if (track_id == mGammas[gamma_index].daughter_ids[j])
                            {
                                space_point = true;
                                // collect results
                                auto xyz = pointsList[i]->XYZ();
                                auto xyz_sigma = pointsList[i]->ErrXYZ();
                                // check if point is in active volume
                                // Determine if edep is within the desired volume
                                mGammas[gamma_index].daughter_reco_sp_x[j].emplace_back(xyz[0]);
                                mGammas[gamma_index].daughter_reco_sp_y[j].emplace_back(xyz[1]);
                                mGammas[gamma_index].daughter_reco_sp_z[j].emplace_back(xyz[2]);
                                mGammas[gamma_index].daughter_reco_sp_x_sigma[j].emplace_back(xyz_sigma[0]);
                                mGammas[gamma_index].daughter_reco_sp_y_sigma[j].emplace_back(xyz_sigma[4]);
                                mGammas[gamma_index].daughter_reco_sp_z_sigma[j].emplace_back(xyz_sigma[8]);
                                mGammas[gamma_index].daughter_reco_sp_chisq[j].emplace_back(pointsList[i]->Chisq());

                                mGammas[gamma_index].daughter_reco_peak_time[j].emplace_back(hit->PeakTime());
                                mGammas[gamma_index].daughter_reco_peak_time_sigma[j].emplace_back(hit->SigmaPeakTime());
                                mGammas[gamma_index].daughter_reco_rms[j].emplace_back(hit->RMS());
                                mGammas[gamma_index].daughter_reco_peak_amplitude[j].emplace_back(hit->PeakAmplitude());
                                mGammas[gamma_index].daughter_reco_peak_amplitude_sigma[j].emplace_back(hit->SigmaPeakAmplitude());
                                mGammas[gamma_index].daughter_reco_summed_adc[j].emplace_back(hit->SummedADC());
                            }
                        }
                    }
                }
                if (space_point) {
                    mGammas[mGammaMap[track_id]].num_reco_points += 1;
                }
            }
            for (size_t i = 0; i < mGammas.size(); i++)
            {
                mGamma = mGammas[i];
                mGammaTree->Fill();
            }
            analyzeEvent();
        }
    }

    void GammaTable::analyzeEvent()
    {
        for(size_t i = 0; i < mGammas.size(); i++)
        {
            bool energy_exists = false;
            for (size_t j = 0; j < mGammaStatistics.energy.size(); j++)
            {
                if (mGammaStatistics.energy[j] == mGammas[i].energy)
                {
                    if (mGammas[i].num_edep_points > 0) {
                        mGammaStatistics.num_gammas_mc[j] += 1;
                    }
                    if (mGammas[i].num_reco_points > 0) {
                        mGammaStatistics.num_gammas_reco[j] += 1;
                    }
                    energy_exists = true;
                }
            }
            if (!energy_exists)
            {
                mGammaStatistics.energy.emplace_back(mGammas[i].energy);
                if (mGammas[i].num_edep_points > 0) {
                    mGammaStatistics.num_gammas_mc.emplace_back(1);
                }
                else {
                    mGammaStatistics.num_gammas_mc.emplace_back(0);
                }
                if (mGammas[i].num_reco_points > 0) {
                    mGammaStatistics.num_gammas_reco.emplace_back(1);
                }
                else {
                    mGammaStatistics.num_gammas_reco.emplace_back(0);
                }
            }
        }
        mGammaStatistics.total_num_gammas += mGammas.size();
    }

    void GammaTable::endJob()
    {
        mGammaStatisticsTree->Fill();
    }
}