/**
 * @file    Voxelizer.h
 * @brief   A class for voxelizing sim/reco.
 * @ingroup Voxelizer
 * @author  Nicholas Carrara (nmcarrara@ucdavis.edu),
**/
#pragma once

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

#include "DetectorGeometry.h"

namespace neutron
{
    // struct for containing voxel information
    struct Voxels
    {
        // voxelization info
        Double_t x_min; Double_t x_max;
        Double_t y_min; Double_t y_max;
        Double_t z_min; Double_t z_max;
        Double_t voxel_size;
        Int_t num_voxels_x;
        Int_t num_voxels_y;
        Int_t num_voxels_z;
        // coordinates
        std::vector<Int_t> x_id;
        std::vector<Int_t> y_id;
        std::vector<Int_t> z_id;

        // features
        std::vector<Double_t> values;

        // labels
        std::vector<Int_t> labels;

        Voxels(Double_t x_min, Double_t x_max, 
            Double_t y_min, Double_t y_max,
            Double_t z_min, Double_t z_max,
            Double_t voxel_size,
            Int_t num_voxels_x, Int_t num_voxels_y, Int_t num_voxels_z,
            std::vector<Int_t> x_id, std::vector<Int_t> y_id, std::vector<Int_t> z_id)
        : x_min(x_min), x_max(x_max), y_min(y_min), y_max(y_max), z_min(z_min), z_max(z_max),
        voxel_size(voxel_size), num_voxels_x(num_voxels_x), num_voxels_y(num_voxels_y),
        num_voxels_z(num_voxels_z), x_id(x_id), y_id(y_id), z_id(z_id)
        {}

        Int_t findVoxel(Int_t x, Int_t y, Int_t z)
        {
            for (Int_t i = 0; i < x_id.size(); i++)
            {
                if (x_id[i] == x and y_id[i] == y and z_id[i] == z)
                {
                    return i;
                }
            }
            return -1;
        }

        void consolidate(const bool discretizeFeatures)
        {
            std::vector<Int_t> x;
            std::vector<Int_t> y;
            std::vector<Int_t> z;
            std::vector<Double_t> val;
            for (Int_t i = 0; i < x_id.size(); i++)
            {
                for (Int_t j = 0; j < x.size(); j++)
                {
                    if (x_id[i] == x[j] and y_id[i] == y[j] and z_id[i] == z[j])
                    {
                        if (!discretizeFeatures)
                        {
                            val[j] += values[i];
                        }
                        break;
                    }
                }
                x.emplace_back(x_id[i]);
                y.emplace_back(y_id[i]);
                z.emplace_back(z_id[i]);
                if (discretizeFeatures)
                {
                    val.emplace_back(1);
                }
                else
                {
                    val.emplace_back(values[i]);
                }
            }
            x_id = x;
            y_id = y;
            z_id = z;
            values = val;
        }
    };

    ///////////////////////////////////////////////////////////////////////////////////////
    // class for constructing voxels
    // take (x,y,z) positions of energy depositions from MC truth, or from reconstructed
    // hits, and voxelize them according to some 'voxel_size' (4.7mm) - comparable to the
    // detector resolution, i.e. the spacing between the wires.
    // coordinates = [i,j,k], features = [edep_energy] or [1], 
    // labels = [mixed==2, cosmic==1, not cosmic==0] or [cosmic==1, not cosmic==0]
    ///////////////////////////////////////////////////////////////////////////////////////
    class Voxelizer
    {
    public:
        Voxelizer();
        Voxelizer(BoundingBox boundingBox);
        ~Voxelizer();

        void setBoundingBox(BoundingBox boundingBox) { fBoundingBox = boundingBox; }

        # function that returns voxels
        Voxels generateVoxels(Double_t voxelSize,
            std::vector<Double_t> x_values,
            std::vector<Double_t> y_values,
            std::vector<Double_t> z_values
        );

        Voxels generateLabeledNeutronCosmicVoxels(
            const Double_t voxelSize,
            const std::vector<Double_t> &neutron_x,
            const std::vector<Double_t> &neutron_y,
            const std::vector<Double_t> &neutron_z,
            const std::vector<Double_t> &neutron_edep_energy,
            const std::vector<Double_t> &muon_x,
            const std::vector<Double_t> &muon_y,
            const std::vector<Double_t> &muon_z,
            const std::vector<Double_t> &muon_edep_energy,
            const bool discretizeFeatures,
            const bool useMixedLabels
        );
        
    private:
        BoundingBox fBoundingBox;
    };

}