/**
 * @file    Neutron.h
 * @brief   A class for holding neutron information collected from LArSoft
 * @ingroup Neutron
 * @author  Nicholas Carrara (nmcarrara@ucdavis.edu),
**/
#pragma once

// art includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Utilities/Exception.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// special utility includes
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"
#include "art_root_io/TFileService.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

// necessary ROOT libraries
#include <TTree.h>
#include <TH1.h>
#include "TH1F.h"
#include "TGeoMaterial.h"
#include "TGeoElement.h"

// std includes
#include <string>
#include <vector>
#include <memory>

namespace neutron 
{

    struct CapturedNeutron
    {
        // meta information
        Int_t event;
        Int_t particle_id;
        // low level trajectory information
        std::vector<Double_t> t;
        std::vector<Double_t> x;
        std::vector<Double_t> y;
        std::vector<Double_t> z;
        std::vector<Double_t> E;
        std::vector<Double_t> px;
        std::vector<Double_t> py;
        std::vector<Double_t> pz;
        // // trajectory volume information
        // DetectorVolume ending_volume;
        // // inelastic scatter information
        // int number_of_inelastic = 0;
        // std::vector<int> inelastic_ids = {};
        // std::vector<int> inelastic_timesteps = {};
        // std::vector<double> inelastic_x = {};
        // std::vector<double> inelastic_y = {};
        // std::vector<double> inelastic_z = {};
        // // gamma information
        // int number_of_gammas = 0;
        // std::vector<int> gamma_ids = {};
        // std::vector<double> gamma_momentum_x = {};
        // std::vector<double> gamma_momentum_y = {};
        // std::vector<double> gamma_momentum_z = {};
        // std::vector<double> gamma_energy = {};
        // // electron information
        // std::vector<int> number_of_electrons = {};
        // std::vector<int> electron_ids = {};
        // std::vector<int> gamma_parent = {};
        // std::vector<double> electron_x = {};
        // std::vector<double> electron_y = {};
        // std::vector<double> electron_z = {};
        // std::vector<double> electron_energy = {};

        CapturedNeutron(){}
        ~CapturedNeutron(){}

        // constructor with simb::MCParticle
        CapturedNeutron(int eventId, simb::MCParticle particle)
        : event(eventId)
        , particle_id(particle.TrackId())
        {
            // collect trajectory information
            for(size_t i = 0; i < particle.NumberTrajectoryPoints(); i++)
            {
                t.emplace_back(particle.T(i));
                x.emplace_back(particle.Vx(i));
                y.emplace_back(particle.Vy(i));
                z.emplace_back(particle.Vz(i));
                E.emplace_back(particle.E(i));
                px.emplace_back(particle.Px(i));
                py.emplace_back(particle.Py(i));
                pz.emplace_back(particle.Pz(i));
            }
        }
    };

}
