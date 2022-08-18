/**
 * @file ParticleMap.h
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2022-08-17
 */
#pragma once
#include "Core.h"

namespace neutron
{
    class ParticleMap
    {
    public:
        ParticleMap();
        ~ParticleMap();

        void ResetMaps();
        void processEvent(const art::ValidHandle<std::vector<simb::MCParticle>>& mcParticles);

        Int_t GetParticlePDG(Int_t trackID)         { return mParticlePDG[trackID]; }
        Int_t GetParticleParentPDG(Int_t trackID)   { return mParticleParentPDG[trackID]; }
        Int_t GetParticleAncestorPDG(Int_t trackID) { return mParticleAncestorPDG[trackID]; }

        Int_t GetParticleParentTrackID(Int_t trackID)   { return mParticleParentTrackID[trackID]; }
        Int_t GetParticleAncestorTrackID(Int_t trackID) { return mParticleAncestorTrackID[trackID]; }

        Int_t GetParticleLevel(Int_t trackID) { return mParticleLevel[trackID]; }

    private:
        art::ServiceHandle<art::TFileService> mTFileService;
        TTree *mParticleMapTTree;

        std::map<Int_t, Int_t> mParticlePDG;
        std::map<Int_t, Int_t> mParticleParentPDG;
        std::map<Int_t, Int_t> mParticleAncestorPDG;
        
        std::map<Int_t, Int_t> mParticleParentTrackID;
        std::map<Int_t, Int_t> mParticleAncestorTrackID;

        std::map<Int_t, Int_t> mParticleLevel;
    };
}