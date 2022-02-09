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

#include <iterator>
#include <map>

#include "DetectorGeometry.h"

namespace neutron
{
    // Voxel Struct (individual voxel); Used in generateSubvolVoxMap()
    struct voxStruct {
        Int_t x_id;
        Int_t y_id;
        Int_t z_id;
    };

    // struct for containing voxel information (all voxels in an event)
    struct Voxels
    {
        Int_t event_id;
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

        // edep ids
        std::vector<std::vector<Int_t>> neutron_edep_ids;
        std::vector<std::vector<Int_t>> muon_edep_ids;

        Voxels(Int_t event) 
        : event_id(event) 
        {}

        Voxels(Int_t event,
            Double_t x_min, Double_t x_max, 
            Double_t y_min, Double_t y_max,
            Double_t z_min, Double_t z_max,
            Double_t voxel_size,
            Int_t num_voxels_x, Int_t num_voxels_y, Int_t num_voxels_z,
            std::vector<Int_t> x_id, std::vector<Int_t> y_id, std::vector<Int_t> z_id
        )
        : event_id(event)
        , x_min(x_min)
        , x_max(x_max)
        , y_min(y_min)
        , y_max(y_max)
        , z_min(z_min)
        , z_max(z_max)
        , voxel_size(voxel_size)
        , num_voxels_x(num_voxels_x)
        , num_voxels_y(num_voxels_y)
        , num_voxels_z(num_voxels_z)
        , x_id(x_id)
        , y_id(y_id)
        , z_id(z_id)
        {}

        Voxels(Int_t event,
            Double_t x_min, Double_t x_max, 
            Double_t y_min, Double_t y_max,
            Double_t z_min, Double_t z_max,
            Double_t voxel_size,
            Int_t num_voxels_x, Int_t num_voxels_y, Int_t num_voxels_z,
            std::vector<Int_t> x_id, std::vector<Int_t> y_id, std::vector<Int_t> z_id,
            std::vector<std::vector<Int_t>> edep_ids
        )
        : event_id(event)
        , x_min(x_min)
        , x_max(x_max)
        , y_min(y_min)
        , y_max(y_max)
        , z_min(z_min)
        , z_max(z_max)
        , voxel_size(voxel_size)
        , num_voxels_x(num_voxels_x)
        , num_voxels_y(num_voxels_y)
        , num_voxels_z(num_voxels_z)
        , x_id(x_id)
        , y_id(y_id)
        , z_id(z_id)
        , neutron_edep_ids(edep_ids)
        {}

        // // Generate sub volume number
        // // We are dividing the space into sub volumes to better search for voxels
        // Int_t factorial(Int_t n)
        // {
        //     if (n == 0 || n == 1) return 1;
        //     else return n * factorial(n-1);
        // }
        // Int_t combination(Int_t n, Int_t r)
        // {
        //     return factorial(n) / (factorial(r) * factorial(n-r));
        // }
        // Int_t generateSubvolNum(Int_t xID, Int_t yID, Int_t zID)
        // {
        //     Int_t x = (int) (xID/100) + 1;
        //     Int_t y = (int) (yID/100) + 1;
        //     Int_t z = (int) (zID/100) + 1;
        //     return combination(x, 1) + combination(x + y + 1, 2) + combination(x + y + z + 2, 3);
        // }

        // // Generating the map between sub volumes and voxels
        // std::map<Int_t, std::vector<voxStruct>> subvolVoxMap;
        // // Need to do this before using findVoxel function
        // void generateSubvolVoxMap()
        // {
        //     subvolVoxMap.clear();

        //     // Filling the map
        //     for (size_t i = 0; i < x_id.size(); i++)
        //     {   
        //         Int_t subvol = generateSubvolNum(x_id[i], y_id[i], z_id[i]);
                
        //         voxStruct vox;
        //         vox.x_id = x_id[i];
        //         vox.y_id = y_id[i];
        //         vox.z_id = z_id[i];

        //         std::map<Int_t, std::vector<voxStruct>>::iterator subvolItr = subvolVoxMap.find(subvol);

        //         if(subvolItr != subvolVoxMap.end())
        //         {
        //             subvolVoxMap[subvol].insert(vox);
        //         }
        //         else
        //         {
        //             subvolVoxMap.insert( make_pair(subvol, std::vector<voxStruct>()) );
        //             subvolVoxMap[subvol].insert(vox);
        //         }
        //     }
        // }

        // //Do generateSubvolVoxMap() before using this
        // Int_t findVoxel(Int_t x, Int_t y, Int_t z)
        // {
        //     Int_t subvol = generateSubvolNum(x, y, z);

        //     std::map<Int_t, std::vector<voxStruct>>::iterator subvolItr = subvolVoxMap.find(subvol);

        //     if (subvolItr != subvolVoxMap.end())
        //     {
        //         for (int i = 0; i < (int) subvolItr->second.size(); i++)
        //         {
        //             if (subvolItr->second[i].x_id == x && subvolItr->second[i].y_id == y && subvolItr->second[i].z_id == z)
        //             {
        //                 return i;
        //             }
        //         }
        //     }

        //     return -1;
        // }

        Int_t findVoxel(Int_t x, Int_t y, Int_t z)
        {
            for (size_t i = 0; i < x_id.size(); i++)
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
            std::vector<std::vector<Int_t>> neutron_ids;
            std::vector<std::vector<Int_t>> muon_ids;
            for (size_t i = 0; i < x_id.size(); i++)
            {
                bool duplicate = false;
                for (size_t j = 0; j < x.size(); j++)
                {
                    if (x_id[i] == x[j] and y_id[i] == y[j] and z_id[i] == z[j])
                    {
                        if (!discretizeFeatures)
                        {
                            val[j] += values[i];
                        }
                        duplicate = true;
                        neutron_ids[j].insert(neutron_ids[j].end(), neutron_edep_ids[i].begin(), neutron_edep_ids[i].end());
                        muon_ids[j].insert(muon_ids[j].end(), muon_edep_ids[i].begin(), muon_edep_ids[i].end());
                    }
                }
                if (duplicate == false)
                {
                    x.emplace_back(x_id[i]);
                    y.emplace_back(y_id[i]);
                    z.emplace_back(z_id[i]);
                    neutron_ids.emplace_back(std::vector<Int_t>({neutron_edep_ids[i]}));
                    muon_ids.emplace_back(std::vector<Int_t>({muon_edep_ids[i]}));
                    if (discretizeFeatures)
                    {
                        val.emplace_back(1);
                    }
                    else
                    {
                        val.emplace_back(values[i]);
                    }
                }   
            }
            x_id = x;
            y_id = y;
            z_id = z;
            values = val;
            neutron_edep_ids = neutron_ids;
            muon_edep_ids = muon_ids;
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

        // function that returns voxels
        Voxels generateVoxels(
            Int_t event,
            Double_t voxelSize,
            std::vector<Double_t> x_values,
            std::vector<Double_t> y_values,
            std::vector<Double_t> z_values,
            std::vector<Int_t>    edep_ids,
            Int_t                 label
        );

        Voxels generateLabeledNeutronCosmicVoxels(
            Int_t event,
            const Double_t voxelSize,
            const std::vector<Double_t> &neutron_x,
            const std::vector<Double_t> &neutron_y,
            const std::vector<Double_t> &neutron_z,
            const std::vector<Double_t> &neutron_edep_energy,
            const std::vector<Int_t>    &neutron_edep_ids,
            const std::vector<Double_t> &muon_x,
            const std::vector<Double_t> &muon_y,
            const std::vector<Double_t> &muon_z,
            const std::vector<Double_t> &muon_edep_energy,
            const std::vector<Int_t>    &muon_edep_ids,
            const bool discretizeFeatures,
            const bool useMixedLabels
        );
        
    private:
        BoundingBox fBoundingBox;
    };

}