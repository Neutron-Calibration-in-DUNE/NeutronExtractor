/**
 * @file DetectorGeometry.cxx
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2022-05-21
 */
#include "DetectorGeometry.h"

namespace neutron 
{
    DetectorGeometry* DetectorGeometry::sInstance{nullptr};
    std::mutex DetectorGeometry::sMutex;

    DetectorGeometry *DetectorGeometry::getInstance(const std::string& name)
    {
        std::lock_guard<std::mutex> lock(sMutex);
        if (sInstance == nullptr)
        {
            sInstance = new DetectorGeometry(name);
        }
        return sInstance;
    }

    std::string DetectorGeometry::GetTPCName(const size_t i) 
    {
        if (i < sTPCNames.size()) { return sTPCNames[i]; }
        else { return sTPCNames[0]; }
    }
    BoundingBox DetectorGeometry::GetTPCBox(const size_t i) 
    {
        if (i < sTPCBoxes.size()) { return sTPCBoxes[i]; }
        else { return sTPCBoxes[0]; }
    }
    BoundingBox DetectorGeometry::GetActiveTPCBox(const size_t i) 
    {
        if (i < sActiveTPCBoxes.size()) { return sActiveTPCBoxes[i]; }
        else { return sActiveTPCBoxes[0]; }
    }
    double DetectorGeometry::GetTPCMass(const size_t i) 
    {
        if (i < sTPCMasses.size()) { return sTPCMasses[i]; }
        else { return sTPCMasses[0]; }
    }
    double DetectorGeometry::GetTPCDriftDistance(const size_t i) 
    {
        if (i < sTPCDriftDistances.size()) { return sTPCDriftDistances[i]; }
        else { return sTPCDriftDistances[0]; }
    }
    
    DetectorGeometry::DetectorGeometry(const std::string name)
    : sName(name)
    {
        sGeometryCore = lar::providerFrom<geo::Geometry>();
        auto const clock_data = 
            art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob();
        sTriggerOffset = trigger_offset(clock_data);

        sWorldName = sGeometryCore->GetWorldVolumeName();
        sWorldBox.setBox(sGeometryCore->WorldBox());

        // create name-volumetype map for world
        sMaterialPOI.SetCoordinates(sWorldBox.x_min,sWorldBox.y_min,sWorldBox.z_min);
        std::string volumeName = sGeometryCore->VolumeName(sMaterialPOI);
        sVolumeTypeMap[volumeName] = VolumeType::World;

        // collect detector info
        sDetectorName = sGeometryCore->DetectorName();
        sDetectorBox.setBox(
            -sGeometryCore->DetHalfWidth(), sGeometryCore->DetHalfWidth(),
            -sGeometryCore->DetHalfHeight(), sGeometryCore->DetHalfHeight(),
            0, sGeometryCore->DetLength()
        );

        // collect cryostat info
        // for now, assuming analysis is done over a single cryostat
        geo::CryostatGeo const& Cryo = sGeometryCore->Cryostat();
        sCryostatName = std::string(Cryo.ID());
        sCryostatBox.setBox(Cryo.Boundaries());

        // create name-volumetype map for cryostat
        sMaterialPOI.SetCoordinates(sCryostatBox.x_min,sCryostatBox.y_min,sCryostatBox.z_min);
        volumeName = sGeometryCore->VolumeName(sMaterialPOI);
        sVolumeTypeMap[volumeName] = VolumeType::Cryostat;

        // iterate over all TPCs
        sNumberOfTPCs  = sGeometryCore->TotalNTPC();
        for (geo::TPCGeo const& TPC : sGeometryCore->IterateTPCs())
        {
            sTPCNames.emplace_back(TPC.ID());
            sTPCBoxes.emplace_back(BoundingBox(TPC.BoundingBox()));
            sActiveTPCBoxes.emplace_back(BoundingBox(TPC.ActiveBoundingBox()));
            sTPCMasses.emplace_back(TPC.ActiveMass());
            sTPCDriftDistances.emplace_back(TPC.DriftDistance());
            // create name-volumetype map for this tpc
            sVolumeTypeMap[sGeometryCore->VolumeName(TPC.GetCenter())] = VolumeType::TPC;
        }

        // find the total TPC and total Active TPC volumes
        findTotalTPCBoxes();
        sTotalTPCMass = sGeometryCore->TotalMass();    
    }

    // get volume information for a point
    DetectorVolume DetectorGeometry::getVolume(std::vector<double> position)
    {
        return getVolume(position[0],position[1],position[2]);
    }
    // get volume information for a point
    DetectorVolume DetectorGeometry::getVolume(double x, double y, double z)
    {
        sMaterialPOI.SetCoordinates(x,y,z);
        std::string volumeName = sGeometryCore->VolumeName(sMaterialPOI);
        VolumeType volumeType = sVolumeTypeMap[volumeName];
        sMaterial = sGeometryService->Material(sMaterialPOI);
        double material = sMaterial->GetZ();
        std::string materialName = sMaterial->GetName();
        return DetectorVolume(volumeType, volumeName, materialName, material);
    }

    void DetectorGeometry::findTotalTPCBoxes()
    {
        double x_min = 0; double x_max = 0;
        double y_min = 0; double y_max = 0;
        double z_min = 0; double z_max = 0;
        for (size_t i = 0; i < sTPCBoxes.size(); i++) {
            if (sTPCBoxes[i].x_min < x_min) x_min = sTPCBoxes[i].x_min;
            if (sTPCBoxes[i].x_max > x_max) x_max = sTPCBoxes[i].x_max;
            if (sTPCBoxes[i].y_min < y_min) y_min = sTPCBoxes[i].y_min;
            if (sTPCBoxes[i].y_max > y_max) y_max = sTPCBoxes[i].y_max;
            if (sTPCBoxes[i].z_min < z_min) z_min = sTPCBoxes[i].z_min;
            if (sTPCBoxes[i].z_max > z_max) z_max = sTPCBoxes[i].z_max;
        }
        sTotalTPCBox.setBox(x_min, x_max, y_min, y_max, z_min, z_max);
        x_min = 0; x_max = 0;
        y_min = 0; y_max = 0;
        z_min = 0; z_max = 0;
        for (size_t i = 0; i < sActiveTPCBoxes.size(); i++) {
            if (sActiveTPCBoxes[i].x_min < x_min) x_min = sActiveTPCBoxes[i].x_min;
            if (sActiveTPCBoxes[i].x_max > x_max) x_max = sActiveTPCBoxes[i].x_max;
            if (sActiveTPCBoxes[i].y_min < y_min) y_min = sActiveTPCBoxes[i].y_min;
            if (sActiveTPCBoxes[i].y_max > y_max) y_max = sActiveTPCBoxes[i].y_max;
            if (sActiveTPCBoxes[i].z_min < z_min) z_min = sActiveTPCBoxes[i].z_min;
            if (sActiveTPCBoxes[i].z_max > z_max) z_max = sActiveTPCBoxes[i].z_max;
        }
        sTotalActiveTPCBox.setBox(x_min, x_max, y_min, y_max, z_min, z_max);
    }
    void DetectorGeometry::FillTTree()
    {
        // add geometry info
        sGeometryTTree = sTFileService->make<TTree>("detector_geometry", "detector_geometry");
        sGeometryTTree->Branch("world_name", &sWorldName);
        sGeometryTTree->Branch("world_box_ranges", &(sWorldBox), "x_min/D:x_max/D:y_min/D:y_max/D:z_min/D:z_max/D");
        sGeometryTTree->Branch("detector_name", &sDetectorName);
        sGeometryTTree->Branch("detector_box_ranges", &(sDetectorBox), "x_min/D:x_max/D:y_min/D:y_max/D:z_min/D:z_max/D");
        sGeometryTTree->Branch("cryostat_name", &sCryostatName);
        sGeometryTTree->Branch("cryostat_box_ranges", &(sCryostatBox), "x_min/D:x_max/D:y_min/D:y_max/D:z_min/D:z_max/D");
        sGeometryTTree->Branch("number_of_tpcs", &sNumberOfTPCs);
        sGeometryTTree->Branch("tpc_names", &sTPCNames);
        for (int i = 0; i < sNumberOfTPCs; i++) {
            sGeometryTTree->Branch(std::string("tpc_"+std::to_string(i)+"_name").c_str(), &(sTPCNames[i]));
            sGeometryTTree->Branch(std::string("tpc_"+std::to_string(i)+"_box_ranges").c_str(), &(sTPCBoxes[i]), "x_min/D:x_max/D:y_min/D:y_max/D:z_min/D:z_max/D");
            sGeometryTTree->Branch(std::string("tpc_"+std::to_string(i)+"_mass").c_str(), &(sTPCMasses[i]));
            sGeometryTTree->Branch(std::string("tpc_"+std::to_string(i)+"_drift_distance").c_str(), &(sTPCDriftDistances[i]));
        }
        sGeometryTTree->Branch("tpc_masses", &sTPCMasses);
        sGeometryTTree->Branch("tpc_drift_distances", &sTPCDriftDistances);
        sGeometryTTree->Branch("total_tpc_box_ranges", &(sTotalTPCBox), "x_min/D:x_max/D:y_min/D:y_max/D:z_min/D:z_max/D");
        sGeometryTTree->Branch("total_active_tpc_box_ranges", &(sTotalActiveTPCBox), "x_min/D:x_max/D:y_min/D:y_max/D:z_min/D:z_max/D");
        sGeometryTTree->Branch("total_tpc_mass", &sTotalTPCMass);
        sGeometryTTree->Fill();
    }
}