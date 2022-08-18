/**
 * @file SingleNeutronCaptures.h
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2022-08-18
 */
#pragma once
#include "Core.h"
#include "ParticleMap.h"

namespace neutron
{
    struct SingleNeutronCapture
    {
        Int_t track_id;
        Int_t parent_track_id;
        Int_t parent_pdg;

        std::vector<Double_t> mc_t;
        std::vector<Double_t> mc_x;
        std::vector<Double_t> mc_y;
        std::vector<Double_t> mc_z;
        std::vector<Double_t> mc_energy;
        std::vector<std::string> mc_process;

        std::vector<Int_t> gamma_indices;
    };

    struct SingleNeutronCaptureGamma
    {

    };

    class SingleNeutronCaptures
    {
    public:
        SingleNeutronCaptures();
        ~SingleNeutronCaptures();

        void ResetSingleNeutronCapture(Int_t number_trajectory_points);

        void processEvent(
            ParticleMap particleMap,
            const art::ValidHandle<std::vector<simb::MCParticle>>& mcParticles
        );

    private:
        art::ServiceHandle<art::TFileService> mTFileService;
        TTree *mSingleNeutronCaptureTTree;
        SingleNeutronCapture mSingleNeutronCapture;
    };
}