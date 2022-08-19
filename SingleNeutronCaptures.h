/**
 * @file SingleNeutronCaptures.h
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2022-08-18
 */
#pragma once
#include "Core.h"
#include "DetectorGeometry.h"
#include "ParticleMap.h"

namespace neutron
{
    struct SingleNeutronCapture
    {
        Int_t event;
        Int_t track_id;
        Int_t parent_track_id;
        Int_t parent_pdg;
        std::string mc_creation_process;
        std::string mc_ending_process;
        bool capture_tpc;
        bool capture_tpc_lar;

        std::vector<Double_t> mc_t;
        std::vector<Double_t> mc_x;
        std::vector<Double_t> mc_y;
        std::vector<Double_t> mc_z;
        std::vector<Double_t> mc_energy;
        std::vector<std::string> mc_volume;
        std::vector<std::string> mc_material;

        std::vector<Int_t> mc_daughter_track_ids;
        std::vector<Int_t> mc_daughter_pdgs;

        std::vector<Int_t> mc_capture_gamma_track_ids;
        std::vector<Double_t> mc_capture_gamma_energies;

        std::vector<Double_t> mc_capture_gamma_t;
        std::vector<Double_t> mc_capture_gamma_x;
        std::vector<Double_t> mc_capture_gamma_y;
        std::vector<Double_t> mc_capture_gamma_z;
        std::vector<Double_t> mc_capture_gamma_energy;
        std::vector<std::string> mc_capture_gamma_volume;
        std::vector<std::string> mc_capture_gamma_material;
        std::vector<Int_t> mc_capture_gamma_track_id;

        std::vector<Double_t> edep_x;
        std::vector<Double_t> edep_y;
        std::vector<Double_t> edep_z;
        std::vector<Double_t> edep_energy;
        std::vector<std::string> edep_volume;
        std::vector<std::string> edep_material;
        std::vector<Int_t> edep_track_id;
        std::vector<Int_t> edep_pdg;
        std::vector<Int_t> edep_capture_gamma_track_ids;
        std::vector<Int_t> edep_capture_gamma_level;


        SingleNeutronCapture()
        {}
        SingleNeutronCapture(Int_t number_trajectory_points)
        {
            mc_t.resize(number_trajectory_points);
            mc_x.resize(number_trajectory_points);
            mc_y.resize(number_trajectory_points);
            mc_z.resize(number_trajectory_points);
            mc_energy.resize(number_trajectory_points);
            mc_volume.resize(number_trajectory_points);
            mc_material.resize(number_trajectory_points);
        }
    };

    class SingleNeutronCaptures
    {
    public:
        SingleNeutronCaptures();
        ~SingleNeutronCaptures();

        void processEvent(
            Int_t event_id,
            ParticleMap particleMap,
            const art::ValidHandle<std::vector<simb::MCParticle>>& mcParticles,
            const art::ValidHandle<std::vector<sim::SimEnergyDeposit>>& mcEnergyDeposits
        );

    private:
        art::ServiceHandle<art::TFileService> mTFileService;
        TTree *mSingleNeutronCaptureTTree;

        DetectorGeometry* mGeometry = DetectorGeometry::getInstance("SingleNeutronCaptures");

        SingleNeutronCapture mSingleNeutronCapture;
        std::vector<SingleNeutronCapture> mSingleNeutronCaptureList;
        std::map<Int_t, Int_t> mSingleNeutronCaptureMap;
    };
}