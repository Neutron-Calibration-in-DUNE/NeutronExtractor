/**
 * @file GammaTable.h
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2022-05-26
 */
#pragma once
#include <string>
#include <vector>
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/Simulation/SimPhotons.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RawData/BeamInfo.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/EndPoint2D.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larsim/Utils/TruthMatchUtils.h"

#include "DetectorGeometry.h"

namespace neutron
{
    struct Gamma
    {
        Int_t track_id;
        Int_t neutron_id;
        Double_t energy;
        Double_t start_x;
        Double_t start_y;
        Double_t start_z;
        Double_t end_x;
        Double_t end_y;
        Double_t end_z;
        std::vector<Int_t> daughter_ids = {};
        std::vector<Int_t> daughter_level = {};
        std::vector<Double_t> daughter_energy = {};
        std::vector<Double_t> daughter_start_x = {};
        std::vector<Double_t> daughter_start_y = {};
        std::vector<Double_t> daughter_start_z = {};
        std::vector<std::vector<Double_t>> daughter_edep_energy = {};
        std::vector<std::vector<Double_t>> daughter_edep_x = {};
        std::vector<std::vector<Double_t>> daughter_edep_y = {};
        std::vector<std::vector<Double_t>> daughter_edep_z = {};
        std::vector<std::vector<Int_t>> daughter_edep_num_electrons = {};
        std::vector<std::vector<Int_t>> daughter_edep_num_photons = {};
        Int_t num_edep_points = 0;
        
        std::vector<std::vector<Double_t>> daughter_reco_sp_x = {};
        std::vector<std::vector<Double_t>> daughter_reco_sp_y = {};
        std::vector<std::vector<Double_t>> daughter_reco_sp_z = {};
        std::vector<std::vector<Double_t>> daughter_reco_sp_x_sigma = {};
        std::vector<std::vector<Double_t>> daughter_reco_sp_y_sigma = {};
        std::vector<std::vector<Double_t>> daughter_reco_sp_z_sigma = {};
        std::vector<std::vector<Double_t>> daughter_reco_sp_chisq = {};
        Int_t num_reco_points = 0;

        std::vector<std::vector<Double_t>> daughter_reco_peak_time = {};
        std::vector<std::vector<Double_t>> daughter_reco_peak_time_sigma = {};
        std::vector<std::vector<Double_t>> daughter_reco_rms = {};
        std::vector<std::vector<Double_t>> daughter_reco_peak_amplitude = {};
        std::vector<std::vector<Double_t>> daughter_reco_peak_amplitude_sigma = {};
        std::vector<std::vector<Double_t>> daughter_reco_summed_adc = {};

        Gamma(
            Int_t track_id, Int_t neutron_id, Double_t energy,
            Double_t start_x, Double_t start_y, Double_t start_z,
            Double_t end_x, Double_t end_y, Double_t end_z
        )
        : track_id(track_id), neutron_id(neutron_id), energy(energy)
        , start_x(start_x), start_y(start_y), start_z(start_z)
        , end_x(end_x), end_y(end_y), end_z(end_z)
        {}
    };

    struct GammaStatistics
    {
        Int_t total_num_gammas = 0;
        std::vector<Double_t> energy = {};
        std::vector<Int_t> num_gammas_mc = {};
        std::vector<Int_t> num_gammas_reco = {};
    };

    class GammaTable
    {
    public:
        GammaTable();
        ~GammaTable();

        void processEvent(
            detinfo::DetectorClocksData const& clockData,
            const art::ValidHandle<std::vector<simb::MCParticle>>& mcParticles,
            const art::ValidHandle<std::vector<sim::SimEnergyDeposit>>& mcEnergyDeposits,
            const art::ValidHandle<std::vector<sim::SimChannel>>& mcChannels,
            const art::ValidHandle<std::vector<recob::SpacePoint>>& recoSpacePoints,
            const art::FindManyP<recob::Hit>& hitSpacePointAssn
        );

        void analyzeEvent();
        void endJob();

    private:
        /// ROOT output through art::TFileService
        /** We will save different TTrees to different TFiles specified 
         *  by the directories for each type.
         */ 
        art::ServiceHandle<art::TFileService> mTFileService;
        TTree *mGammaTree;
        TTree *mGammaStatisticsTree;
        // geometry information
        DetectorGeometry* fGeometry = DetectorGeometry::getInstance("GammaTable");

        Gamma mGamma;
        std::vector<Gamma> mGammas;
        std::map<int, int> mGammaMap;

        GammaStatistics mGammaStatistics;
    };
}