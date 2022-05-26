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
                            std::cout << "Gamma: " << particle.TrackId() << ", " << particle.E() << std::endl;
                            mGammas.emplace_back(
                                Gamma(
                                    particle.TrackId(), particle.Mother(), particle.E(), 
                                    particle.Vx(), particle.Vy(), particle.Vz(),
                                    particle.EndX(), particle.EndY(), particle.EndZ()
                                )
                            );
                        }
                    }
                }
            }
        }
    }
}