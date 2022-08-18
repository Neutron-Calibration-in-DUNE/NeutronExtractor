/**
 * @file ParticleMap.cxx
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2022-08-17
 */
#include "ParticleMap.h"

namespace neutron
{
    ParticleMap::ParticleMap()
    {
        mParticleMapTTree = mTFileService->make<TTree>("particle_map", "particle_map");
        mParticleMapTTree->Branch("pdg", &mParticlePDG);
        mParticleMapTTree->Branch("parent_pdg", &mParticleParentPDG);
        mParticleMapTTree->Branch("ancestor_pdg", &mParticleAncestorPDG);
        mParticleMapTTree->Branch("parent_track_id", &mParticleParentTrackID);
        mParticleMapTTree->Branch("ancestor_track_id", &mParticleAncestorTrackID);
        mParticleMapTTree->Branch("level", &mParticleLevel);
    }
    ParticleMap::~ParticleMap()
    {
    }

    void ParticleMap::ResetMaps()
    {
        mParticlePDG.clear();
        mParticleParentPDG.clear();
        mParticleAncestorPDG.clear();
        mParticleParentTrackID.clear();
        mParticleAncestorTrackID.clear();
        mParticleLevel.clear();

        mParticlePDG[0] = 0;
        mParticleParentPDG[0] = 0;
        mParticleAncestorPDG[0] = 0;
        mParticleParentTrackID[0] = 0;
        mParticleAncestorTrackID[0] = 0;
        mParticleLevel[0] = 0;
    }

    void ParticleMap::processEvent(const art::ValidHandle<std::vector<simb::MCParticle>>& mcParticles)
    {
        ResetMaps();
        if (mcParticles.isValid())
        {
            Int_t level = 0;
            Int_t mother = 0;
            Int_t track_id = 0;

            for (auto particle : *mcParticles)
            {
                mParticlePDG[particle.TrackId()] = particle.PdgCode();
                mParticlePDG[particle.Mother()] = mParticlePDG[particle.Mother()];
                mParticleParentTrackID[particle.TrackId()] = particle.Mother();

                if (particle.Mother() == 0) 
                {
                    mParticleAncestorPDG[particle.TrackId()] = 0;
                    mParticleAncestorTrackID[particle.TrackId()] = 0;
                    mParticleLevel[particle.TrackId()] = 0;
                }
                else
                {
                    level = 0;
                    mother = particle.Mother();
                    while (mother != 0)
                    {
                        level += 1;
                        track_id = mother;
                        mother = mParticleParentTrackID[track_id];
                    }
                    mParticleAncestorPDG[particle.TrackId()] = mParticlePDG[track_id];
                    mParticleAncestorTrackID[particle.TrackId()] = track_id;
                    mParticleLevel[particle.TrackId()] = level;
                }
            }
        }
        mParticleMapTTree->Fill();
    }
}