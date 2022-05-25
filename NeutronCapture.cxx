/**
 * @file NeutronCapture.cxx
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2022-05-23
 */
#include "NeutronCapture.h"

namespace neutron
{
    NeutronCapture::NeutronCapture()
    {
        mNeutronCaptureRecoTree = fTFileService->make<TTree>("neutron_capture_reco", "neutron_capture_reco");
        mNeutronCaptureRecoTree->Branch("space_point_x", &mNeutronCaptureReco.SpacePointX);
        mNeutronCaptureRecoTree->Branch("space_point_y", &mNeutronCaptureReco.SpacePointX);
        mNeutronCaptureRecoTree->Branch("space_point_z", &mNeutronCaptureReco.SpacePointX);
        mNeutronCaptureRecoTree->Branch("space_point_chisq_x", &mNeutronCaptureReco.SpacePointChiSqX);
        mNeutronCaptureRecoTree->Branch("space_point_chisq_y", &mNeutronCaptureReco.SpacePointChiSqY);
        mNeutronCaptureRecoTree->Branch("space_point_chisq_z", &mNeutronCaptureReco.SpacePointChiSqZ);
        mNeutronCaptureRecoTree->Branch("space_point_chisq",   &mNeutronCaptureReco.SpacePointChiSq);

        mNeutronCaptureRecoTree->Branch("edep_track_id",    &mNeutronCaptureReco.EdepTrackID);
        mNeutronCaptureRecoTree->Branch("neutron_track_id", &mNeutronCaptureReco.NeutronTrackID);
        mNeutronCaptureRecoTree->Branch("gamma_track_id",   &mNeutronCaptureReco.GammaTrackID);
        mNeutronCaptureRecoTree->Branch("gamma_energy",     &mNeutronCaptureReco.GammaEnergy);
        mNeutronCaptureRecoTree->Branch("edep_energy",     &mNeutronCaptureReco.EdepEnergy);
        mNeutronCaptureRecoTree->Branch("edep_num_electrons",     &mNeutronCaptureReco.EdepNumElectrons);
        mNeutronCaptureRecoTree->Branch("edep_num_photons",     &mNeutronCaptureReco.EdepNumPhotons);

        mNeutronCaptureRecoTree->Branch("peak_time",        &mNeutronCaptureReco.PeakTime);
        mNeutronCaptureRecoTree->Branch("sigma_peak_time",  &mNeutronCaptureReco.SigmaPeakTime);
        mNeutronCaptureRecoTree->Branch("rms",              &mNeutronCaptureReco.RMS);
        mNeutronCaptureRecoTree->Branch("peak_amplitude",   &mNeutronCaptureReco.PeakAmplitude);
        mNeutronCaptureRecoTree->Branch("sigma_peak_amplitude", &mNeutronCaptureReco.SigmaPeakAmplitude);
        mNeutronCaptureRecoTree->Branch("summed_adc",       &mNeutronCaptureReco.SummedADC);
    }

    NeutronCapture::~NeutronCapture()
    {}

    void NeutronCapture::processEvent(
        detinfo::DetectorClocksData const& clockData,
        const art::ValidHandle<std::vector<simb::MCParticle>>& mcParticles,
        const art::ValidHandle<std::vector<sim::SimEnergyDeposit>>& mcEnergyDeposits,
        const art::ValidHandle<std::vector<sim::SimChannel>>& mcChannels,
        const art::ValidHandle<std::vector<recob::SpacePoint>>& recoSpacePoints,
        const art::FindManyP<recob::Hit>& hitPandoraSPsAssn //to associate space points from pandora to hits
    )
    {
        NeutronCaptureReco NeutronCaptureReco;
        if (mcParticles.isValid() and mcChannels.isValid() and recoSpacePoints.isValid())
        {
            /**
             * We first iterate through all particles and create a map of 
             * parent-daughter pairs for track ids.  This way we can search
             * recursively for ancestors of each particle.
             */
            std::map<Int_t, Int_t> parentDaughterMap;
            std::map<Int_t, Int_t> particlePDGMap;

            std::vector<int> neutron_captures;
            std::vector<std::vector<int>> gamma_ids;
            std::vector<std::vector<double>> gamma_energy;

            std::map<Int_t, Int_t> neutronMap;
            std::map<Int_t, Int_t> gammaMap;
            std::map<Int_t, bool> neutronCapture;

            for (auto particle : *mcParticles)
            {
                parentDaughterMap[particle.TrackId()] = particle.Mother();
                particlePDGMap[particle.TrackId()] = particle.PdgCode();
            }
            for (auto particle : *mcParticles)
            {
                // check if the particle is a neutron
                if (particle.PdgCode() == 2112)
                {
                    if (particle.EndProcess() == "nCapture")
                    {
                        neutron_captures.emplace_back(particle.TrackId());
                        gamma_ids.emplace_back(std::vector<int>());
                        gamma_energy.emplace_back(std::vector<double>());
                    }
                }
                // check if the particle is a gamma
                else if (particle.PdgCode() == 22)
                {
                    for(size_t i = 0; i < neutron_captures.size(); i++)
                    {
                        if (neutron_captures[i] == particle.Mother())
                        {
                            gamma_ids[i].emplace_back(particle.TrackId());
                            gamma_energy[i].emplace_back(particle.E());
                        }
                    }
                }
                // otherwise see if the particle is from a neutron capture
                else
                {
                    Int_t mother = particle.Mother();
                    Int_t track_id = particle.TrackId();
                    while (mother != 0)
                    {
                        for(size_t i = 0; i < neutron_captures.size(); i++)
                        {
                            if (neutron_captures[i] == mother)
                            {
                                neutronMap[particle.TrackId()] = i;
                                for (size_t j = 0; j < gamma_ids[i].size(); j++)
                                {
                                    if (gamma_ids[i][j] == track_id)
                                    {
                                        gammaMap[particle.TrackId()] = j;
                                        neutronCapture[particle.TrackId()] = true;
                                        break;
                                    }
                                }
                                break;
                            }
                        }
                        track_id = mother;
                        mother = parentDaughterMap[track_id];
                    }
                }
            }
            std::cout << "HERE" << std::endl;
            std::vector<art::Ptr<recob::SpacePoint>> pointsList;
            art::fill_ptr_vector(pointsList, recoSpacePoints);            
            for (size_t i = 0; i < pointsList.size(); i++)
            {
                auto& spsHit = hitPandoraSPsAssn.at(i);
                for (auto hit : spsHit)
                {  
                    Int_t track_id = TruthMatchUtils::TrueParticleID(
                        clockData, hit, false
                    );
                    // check that track_id is present in parentDaughterMap
                    if (parentDaughterMap.find(track_id) == parentDaughterMap.end())
                    {
                        continue;
                    }
                    if (particlePDGMap.find(track_id) == particlePDGMap.end())
                    {
                        continue;
                    }
                    if (neutronCapture[track_id])
                    {
                        // collect results
                        auto xyz = pointsList[i]->XYZ();
                        auto chisqxyz = pointsList[i]->ErrXYZ();
                        auto chisq = pointsList[i]->Chisq();
                        // check if point is in active volume
                        // Determine if edep is within the desired volume
                        DetectorVolume edep_volume = fGeometry->getVolume(
                            xyz[0], xyz[1], xyz[2]
                        );
                        
                        NeutronCaptureReco.SpacePointX.emplace_back(xyz[0]);
                        NeutronCaptureReco.SpacePointY.emplace_back(xyz[1]);
                        NeutronCaptureReco.SpacePointZ.emplace_back(xyz[2]);
                        NeutronCaptureReco.SpacePointChiSqX.emplace_back(chisqxyz[0]);
                        NeutronCaptureReco.SpacePointChiSqY.emplace_back(chisqxyz[4]);
                        NeutronCaptureReco.SpacePointChiSqZ.emplace_back(chisqxyz[8]);
                        NeutronCaptureReco.SpacePointChiSq.emplace_back(chisq);

                        NeutronCaptureReco.EdepTrackID.emplace_back(track_id);
                        NeutronCaptureReco.NeutronTrackID.emplace_back(neutron_captures[neutronMap[track_id]]);
                        NeutronCaptureReco.GammaTrackID.emplace_back(gamma_ids[neutronMap[track_id]][gammaMap[track_id]]);
                        NeutronCaptureReco.GammaEnergy.emplace_back(gamma_energy[neutronMap[track_id]][gammaMap[track_id]]);
                        for (auto energyDeposit : *mcEnergyDeposits)
                        {
                            if (energyDeposit.TrackID() == track_id)
                            {
                                NeutronCaptureReco.EdepEnergy.emplace_back(energyDeposit.Energy());
                                NeutronCaptureReco.EdepNumElectrons.emplace_back(energyDeposit.NumElectrons());
                                NeutronCaptureReco.EdepNumPhotons.emplace_back(energyDeposit.NumPhotons());
                            }   
                            break;
                        }
                        NeutronCaptureReco.PeakTime.emplace_back(hit->PeakTime());
                        NeutronCaptureReco.SigmaPeakTime.emplace_back(hit->SigmaPeakTime());
                        NeutronCaptureReco.RMS.emplace_back(hit->RMS());
                        NeutronCaptureReco.PeakAmplitude.emplace_back(hit->PeakAmplitude());
                        NeutronCaptureReco.SigmaPeakAmplitude.emplace_back(hit->SigmaPeakAmplitude());
                        NeutronCaptureReco.SummedADC.emplace_back(hit->SummedADC());

                        break;
                    }
                }
            }
        }
        std::cout << "HERE" << std::endl;
        mNeutronCaptureReco = NeutronCaptureReco;
        mNeutronCaptureRecoTree->Fill();
        std::cout << "HERE" << std::endl;
    }
}