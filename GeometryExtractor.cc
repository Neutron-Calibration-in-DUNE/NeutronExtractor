/**
 * @file    GeometryExtractor.cc
 * @brief   A module for extracting geometry information from art root files.
 * @ingroup GeometryExtractor
 * @author  Nicholas Carrara (nmcarrara@ucdavis.edu),
**/

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

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"

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

namespace geometry
{
    ///////////////////////////////////////////////////////////////////////////////////////
    // struct for detector volume information
    ///////////////////////////////////////////////////////////////////////////////////////
    struct DetectorVolume
    {
        VolumeType volume_type;
        std::string volume_name;
        std::string material_name;
        double material;
        DetectorVolume() {}
        DetectorVolume(VolumeType volumeType, std::string volumeName, 
            std::string materialName, double material)
        : volume_type(volumeType)
        , volume_name(volumeName)
        , material_name(materialName)
        , material(material)
        {}
    };
    ///////////////////////////////////////////////////////////////////////////////////////
    // struct for bounding boxes
    ///////////////////////////////////////////////////////////////////////////////////////
    struct BoundingBox
    {
        double x_min = 0; double x_max = 0;
        double y_min = 0; double y_max = 0;
        double z_min = 0; double z_max = 0;

        double width()  { return x_max - x_min; }
        double height() { return y_max - y_min; }
        double length() { return z_max - z_min; }

        void setBox(geo::BoxBoundedGeo const& Box) {
            x_min = Box.MinX(); x_max = Box.MaxX();
            y_min = Box.MinY(); y_max = Box.MaxY();
            z_min = Box.MinZ(); z_max = Box.MaxZ();
        }
        void setBox(double x_min, double x_max,
                    double y_min, double y_max,
                    double z_min, double z_max)
        {
            x_min = x_min; x_max = x_max;
            y_min = y_min; y_max = y_max;
            z_min = z_min; z_max = z_max;
        }
        BoundingBox() {}
        BoundingBox(double xs[2], double ys[2], double zs[2])
        {
            x_min = xs[0]; x_max = xs[1];
            y_min = ys[0]; y_max = ys[1];
            z_min = zs[0]; z_max = zs[1];
        }
        BoundingBox(double vals[6])
        {
            x_min = vals[0]; x_max = vals[1];
            y_min = vals[2]; y_max = vals[3];
            z_min = vals[4]; z_max = vals[4];
        }
        BoundingBox(double xmin, double xmax,
                    double ymin, double ymax,
                    double zmin, double zmax)
        {
            x_min = xmin; x_max = xmax;
            y_min = ymin; y_max = ymax;
            z_min = zmin; z_max = zmax;
        }
        BoundingBox(geo::BoxBoundedGeo const& Box) {
            x_min = Box.MinX(); x_max = Box.MaxX();
            y_min = Box.MinY(); y_max = Box.MaxY();
            z_min = Box.MinZ(); z_max = Box.MaxZ();
        }
    };
    ///////////////////////////////////////////////////////////////////////////////////////
    // class for storting geometry information
    ///////////////////////////////////////////////////////////////////////////////////////
    class DetectorGeometry
    {
    public:
        DetectorGeometry();
        ~DetectorGeometry(){}

        // getters
        std::string GetWorldName()          { return fWorldName; }
        BoundingBox GetWorldBox()           { return fWorldBox; }
        std::string GetDetectorName()       { return fDetectorName; }
        BoundingBox GetDetectorBox()        { return fDetectorBox; }
        std::string GetCryostatName()       { return fCryostatName; }
        BoundingBox GetCryostatBox()        { return fCryostatBox; }
        int GetNumberOfTPCs()               { return fNumberOfTPCs; }
        std::vector<std::string> GetTPCNames()  { return fTPCNames; }
        std::string GetTPCName(const size_t i) {
            if (i < fTPCNames.size()) { return fTPCNames[i]; }
            else { return fTPCNames[0]; }
        }
        BoundingBox GetTPCBox(const size_t i) {
            if (i < fTPCBoxes.size()) { return fTPCBoxes[i]; }
            else { return fTPCBoxes[0]; }
        }
        BoundingBox GetActiveTPCBox(const size_t i) {
            if (i < fActiveTPCBoxes.size()) { return fActiveTPCBoxes[i]; }
            else { return fActiveTPCBoxes[0]; }
        }
        std::vector<double> GetTPCMasses() { return fTPCMasses; }
        double GetTPCMass(const size_t i) {
            if (i < fTPCMasses.size()) { return fTPCMasses[i]; }
            else { return fTPCMasses[0]; }
        }
        std::vector<double> GetTPCDriftDistances() { return fTPCDriftDistances; }
        double GetTPCDriftDistance(const size_t i) {
            if (i < fTPCDriftDistances.size()) { return fTPCDriftDistances[i]; }
            else { return fTPCDriftDistances[0]; }
        }
        BoundingBox GetTotalTPCBox()        { return fTotalTPCBox; }
        BoundingBox GetTotalActiveTPCBox()  { return fTotalActiveTPCBox; }
        double GetTotalTPCMass()            { return fTotalTPCMass; }
        // get volume information for a point
        DetectorVolume getVolume(std::vector<double> position);
        // function for finding total tpc volumes
        void findTotalTPCBoxes();
        // fill the geometry ttree
        void FillTTree();
        
    private:
        ////////////////////////////////////////////////
        // Information which is automatically stored
        // meta variables
        art::ServiceHandle<geo::Geometry> fGeometryService;
        geo::GeometryCore const* fGeometryCore;
        // ROOT 
        art::ServiceHandle<art::TFileService> fTFileService;
        TTree *fGeometryTree;
        size_t fTriggerOffset;
        // map from volume names to volume type
        std::map<std::string,VolumeType> fVolumeTypeMap;
        // world volume
        std::string fWorldName;
        BoundingBox fWorldBox;
        // detector volume
        std::string fDetectorName;
        BoundingBox fDetectorBox;
        // cryostat volume
        std::string fCryostatName;
        BoundingBox fCryostatBox;
        // tpc volumes
        int fNumberOfTPCs;
        std::vector<std::string> fTPCNames;
        std::vector<BoundingBox> fTPCBoxes;
        std::vector<BoundingBox> fActiveTPCBoxes;
        std::vector<double> fTPCMasses;
        std::vector<double> fTPCDriftDistances;
        // full tpc volume
        BoundingBox fTotalTPCBox;
        BoundingBox fTotalActiveTPCBox;
        double fTotalTPCMass;
        ////////////////////////////////////////////////
        // detector material variables
        ////////////////////////////////////////////////
        // we will need to ask Geant4 about material 
        // properties for the detector volume
        // at each point of interest.  This requires holding 
        // this information in a
        // TGeoMaterial object, which is part of ROOT.
        const TGeoMaterial *fMaterial;
        geo::Point_t fMaterialPOI;
    };
    ///////////////////////////////////////////////////////////////////////////////////////
    DetectorGeometry::DetectorGeometry()
    {
        // set up the geometry interface
        fGeometryCore = lar::providerFrom<geo::Geometry>();
        // initialize TTrees
        fGeometryTree = fTFileService->make<TTree>("Geometry", "Geometry");
        // get detector clock data
        auto const clock_data = 
            art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob();
        fTriggerOffset = trigger_offset(clock_data);

        // collect world info
        fWorldName = fGeometryCore->GetWorldVolumeName();
        fWorldBox.setBox(fGeometryCore->WorldBox());
        // create name-volumetype map for world
        fMaterialPOI.SetCoordinates(fWorldBox.x_min,fWorldBox.y_min,fWorldBox.z_min);
        std::string volumeName = fGeometryCore->VolumeName(fMaterialPOI);
        fVolumeTypeMap[volumeName] = VolumeType::World;
        // collect detector info
        fDetectorName = fGeometryCore->DetectorName();
        fDetectorBox.setBox(-fGeometryCore->DetHalfWidth(), fGeometryCore->DetHalfWidth(),
                            -fGeometryCore->DetHalfHeight(), fGeometryCore->DetHalfHeight(),
                            0, fGeometryCore->DetLength());
        // collect cryostat info
        // for now, assuming analysis is done over a single cryostat
        geo::CryostatGeo const& Cryo = fGeometryCore->Cryostat();
        fCryostatName = std::string(Cryo.ID());
        fCryostatBox.setBox(Cryo.Boundaries());
        // create name-volumetype map for cryostat
        fMaterialPOI.SetCoordinates(fCryostatBox.x_min,fCryostatBox.y_min,fCryostatBox.z_min);
        volumeName = fGeometryCore->VolumeName(fMaterialPOI);
        fVolumeTypeMap[volumeName] = VolumeType::Cryostat;
        // iterate over all TPCs
        fNumberOfTPCs  = fGeometryCore->TotalNTPC();
        for (geo::TPCGeo const& TPC : fGeometryCore->IterateTPCs())
        {
            fTPCNames.emplace_back(TPC.ID());
            fTPCBoxes.emplace_back(BoundingBox(TPC.BoundingBox()));
            fActiveTPCBoxes.emplace_back(BoundingBox(TPC.ActiveBoundingBox()));
            fTPCMasses.emplace_back(TPC.ActiveMass());
            fTPCDriftDistances.emplace_back(TPC.DriftDistance());
            // create name-volumetype map for this tpc
            fVolumeTypeMap[fGeometryCore->VolumeName(TPC.GetCenter())] = VolumeType::TPC;
        }
        // find the total TPC and total Active TPC volumes
        findTotalTPCBoxes();
        fTotalTPCMass = fGeometryCore->TotalMass();    
    }
    // get volume information for a point
    DetectorVolume DetectorGeometry::getVolume(std::vector<double> position)
    {

        fMaterialPOI.SetCoordinates(position[0],position[1],position[2]);
        // get the volume information
        //std::cout << "volname: " << fGeometryCore->VolumeName(fMaterialPOI) << std::endl;
        std::string volumeName = fGeometryCore->VolumeName(fMaterialPOI);
        VolumeType volumeType = fVolumeTypeMap[volumeName];
        // get the current material information
        fMaterial = fGeometryService->Material(fMaterialPOI);
        double material = fMaterial->GetZ();
        std::string materialName = fMaterial->GetName();
        //std::cout << "mat name: " << fMaterial->GetName() << std::endl;
        // return the constructed volume 
        return DetectorVolume(volumeType, volumeName, materialName, material);
    }
    // get total tpc volume information
    void DetectorGeometry::findTotalTPCBoxes()
    {
        double x_min = 0; double x_max = 0;
        double y_min = 0; double y_max = 0;
        double z_min = 0; double z_max = 0;
        for (size_t i = 0; i < fTPCBoxes.size(); i++) {
            if (fTPCBoxes[i].x_min < x_min) x_min = fTPCBoxes[i].x_min;
            if (fTPCBoxes[i].x_max > x_max) x_max = fTPCBoxes[i].x_max;
            if (fTPCBoxes[i].y_min < y_min) y_min = fTPCBoxes[i].y_min;
            if (fTPCBoxes[i].y_max > y_max) y_max = fTPCBoxes[i].y_max;
            if (fTPCBoxes[i].z_min < z_min) z_min = fTPCBoxes[i].z_min;
            if (fTPCBoxes[i].z_max > z_max) z_max = fTPCBoxes[i].z_max;
        }
        fTotalTPCBox.setBox(x_min, x_max, y_min, y_max, z_min, z_max);
        x_min = 0; x_max = 0;
        y_min = 0; y_max = 0;
        z_min = 0; z_max = 0;
        for (size_t i = 0; i < fActiveTPCBoxes.size(); i++) {
            if (fActiveTPCBoxes[i].x_min < x_min) x_min = fActiveTPCBoxes[i].x_min;
            if (fActiveTPCBoxes[i].x_max > x_max) x_max = fActiveTPCBoxes[i].x_max;
            if (fActiveTPCBoxes[i].y_min < y_min) y_min = fActiveTPCBoxes[i].y_min;
            if (fActiveTPCBoxes[i].y_max > y_max) y_max = fActiveTPCBoxes[i].y_max;
            if (fActiveTPCBoxes[i].z_min < z_min) z_min = fActiveTPCBoxes[i].z_min;
            if (fActiveTPCBoxes[i].z_max > z_max) z_max = fActiveTPCBoxes[i].z_max;
        }
        fTotalActiveTPCBox.setBox(x_min, x_max, y_min, y_max, z_min, z_max);
    }
    void DetectorGeometry::FillTTree()
    {
        // add geometry info
        fGeometryTree->Branch("world_name", &fWorldName);
        fGeometryTree->Branch("world_box_ranges", &(fWorldBox), "x_min/D:x_max/D:y_min/D:y_max/D:z_min/D:z_max/D");
        fGeometryTree->Branch("detector_name", &fDetectorName);
        fGeometryTree->Branch("detector_box_ranges", &(fDetectorBox), "x_min/D:x_max/D:y_min/D:y_max/D:z_min/D:z_max/D");
        fGeometryTree->Branch("cryostat_name", &fCryostatName);
        fGeometryTree->Branch("cryostat_box_ranges", &(fCryostatBox), "x_min/D:x_max/D:y_min/D:y_max/D:z_min/D:z_max/D");
        fGeometryTree->Branch("number_of_tpcs", &fNumberOfTPCs);
        fGeometryTree->Branch("tpc_names", &fTPCNames);
        for (int i = 0; i < fNumberOfTPCs; i++) {
            fGeometryTree->Branch(std::string("tpc_"+std::to_string(i)+"_name").c_str(), &(fTPCNames[i]));
            fGeometryTree->Branch(std::string("tpc_"+std::to_string(i)+"_box_ranges").c_str(), &(fTPCBoxes[i]), "x_min/D:x_max/D:y_min/D:y_max/D:z_min/D:z_max/D");
            fGeometryTree->Branch(std::string("tpc_"+std::to_string(i)+"_mass").c_str(), &(fTPCMasses[i]));
            fGeometryTree->Branch(std::string("tpc_"+std::to_string(i)+"_drift_distance").c_str(), &(fTPCDriftDistances[i]));
        }
        fGeometryTree->Branch("tpc_masses", &fTPCMasses);
        fGeometryTree->Branch("tpc_drift_distances", &fTPCDriftDistances);
        fGeometryTree->Branch("total_tpc_box_ranges", &(fTotalTPCBox), "x_min/D:x_max/D:y_min/D:y_max/D:z_min/D:z_max/D");
        fGeometryTree->Branch("total_active_tpc_box_ranges", &(fTotalActiveTPCBox), "x_min/D:x_max/D:y_min/D:y_max/D:z_min/D:z_max/D");
        fGeometryTree->Branch("total_tpc_mass", &fTotalTPCMass);
        fGeometryTree->Fill();
    }
}