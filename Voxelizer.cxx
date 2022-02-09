#include "Voxelizer.h"

namespace neutron 
{
    Voxelizer::Voxelizer()
    {}

    Voxelizer::Voxelizer(BoundingBox boundingBox)
    : fBoundingBox(boundingBox)
    {

    }

    Voxelizer::~Voxelizer()
    {}

    Voxels Voxelizer::generateVoxels(
        Int_t event,
        Double_t voxelSize,
        std::vector<Double_t> x_values,
        std::vector<Double_t> y_values,
        std::vector<Double_t> z_values,
        std::vector<Int_t>    edep_ids,
        Int_t                 label
    )
    {
        // set up variables
        Double_t xMin = fBoundingBox.x_min;
        Double_t yMin = fBoundingBox.y_min;
        Double_t zMin = fBoundingBox.z_min;
        Double_t xRange = (fBoundingBox.x_max - fBoundingBox.x_min);
        Double_t yRange = (fBoundingBox.y_max - fBoundingBox.y_min);
        Double_t zRange = (fBoundingBox.z_max - fBoundingBox.z_min);
        Int_t numXVoxels = int(xRange / voxelSize);
        Int_t numYVoxels = int(yRange / voxelSize);
        Int_t numZVoxels = int(zRange / voxelSize);
        // to voxelize, first normalize: x' = (x - x_min)/xRange,
        // then multiply by number of voxels: x'' = x' * numXVoxels,
        // then cast to an int, x_voxel = int(x'')
        std::vector<Int_t> x_voxels(x_values.size());
        std::vector<Int_t> y_voxels(y_values.size());
        std::vector<Int_t> z_voxels(z_values.size());
        
        std::vector<std::vector<Int_t>> voxel_edep_ids(x_values.size());
        // iterate through values
        for (size_t i = 0; i < x_values.size(); i++)
        {
            x_voxels[i] = int((x_values[i] - xMin)/voxelSize);
            y_voxels[i] = int((y_values[i] - yMin)/voxelSize);
            z_voxels[i] = int((z_values[i] - zMin)/voxelSize);
            int temp = std::static_cast<int>(i);
            voxel_edep_ids[i] = std::vector<Int_t>({temp});
        }
        Voxels voxels(
            event,
            xMin, fBoundingBox.x_max, 
            yMin, fBoundingBox.y_max, 
            zMin, fBoundingBox.z_max, 
            voxelSize,
            numXVoxels, numYVoxels, numZVoxels,
            x_voxels, y_voxels, z_voxels);
        if (label == 0) {
            voxels.neutron_edep_ids = voxel_edep_ids;
        }
        else {
            voxels.muon_edep_ids = voxel_edep_ids;
        }
        return voxels;
    }

    Voxels Voxelizer::generateLabeledNeutronCosmicVoxels(
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
    )
    {
        Voxels neutronVoxels = generateVoxels(event, voxelSize,
            neutron_x, neutron_y, neutron_z, neutron_edep_ids, 0);
        Voxels muonVoxels = generateVoxels(event, voxelSize,
            muon_x, muon_y, muon_z, muon_edep_ids, 1);
        
        neutronVoxels.values = neutron_edep_energy;
        muonVoxels.values = muon_edep_energy;

        // consolidate edep energy values
        neutronVoxels.consolidate(discretizeFeatures);
        muonVoxels.consolidate(discretizeFeatures);
        
        for (size_t i = 0; i < neutronVoxels.x_id.size(); i++)
        {
            neutronVoxels.labels.emplace_back(0);
        }

        // construct training set
        //neutronVoxels.generateSubvolVoxMap(); // Do this before using findVoxel()

        for (size_t i = 0; i < muonVoxels.x_id.size(); i++)
        {
            // search neutron list for muon voxel
            Int_t index = neutronVoxels.findVoxel(
                muonVoxels.x_id[i],
                muonVoxels.y_id[i],
                muonVoxels.z_id[i]
            );
            // if muon voxel is in the neutron voxel list, then...
            if (index != -1)
            {
                // decide what to do with the energy
                if (!discretizeFeatures)
                {
                    neutronVoxels.values[index] += muonVoxels.values[i];
                }
                // decide what to do with the labels
                if (useMixedLabels)
                {
                    neutronVoxels.labels[index] = 2;
                }
                else
                {
                    if (muon_edep_energy[i] >= neutron_edep_energy[index])
                    {
                        neutronVoxels.labels[index] = 1;
                    }
                }
                // add the edep ids
                neutronVoxels.muon_edep_ids[index] = muonVoxels.muon_edep_ids[i];
            }
            else
            {
                neutronVoxels.x_id.emplace_back(muonVoxels.x_id[i]);
                neutronVoxels.y_id.emplace_back(muonVoxels.y_id[i]);
                neutronVoxels.z_id.emplace_back(muonVoxels.z_id[i]);
                neutronVoxels.values.emplace_back(muonVoxels.values[i]);
                neutronVoxels.labels.emplace_back(1);
                neutronVoxels.muon_edep_ids.emplace_back(muonVoxels.muon_edep_ids[i];
            }
        }
        return neutronVoxels;
    }
}