/**
 * @file    GeometryExtractor.h
 * @brief   A class for extracting geometry information from art root files.
 * @ingroup GeometryExtractor
 * @author  Nicholas Carrara (nmcarrara@ucdavis.edu),
**/
#pragma once
#include "Core.h"
namespace neutron 
{
    enum MaterialList {

    };
    enum VolumeType {
        World,
        Cryostat,
        TPC,
    };

    struct DetectorVolume
    {
        VolumeType volume_type;
        std::string volume_name;
        std::string material_name;
        double material;

        DetectorVolume() {}
        DetectorVolume(
            VolumeType volumeType, std::string volumeName, 
            std::string materialName, double material
        )
        : volume_type(volumeType)
        , volume_name(volumeName)
        , material_name(materialName)
        , material(material)
        {}
    };
    // struct for bounding boxes
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
        void setBox(double xmin, double xmax,
                    double ymin, double ymax,
                    double zmin, double zmax)
        {
            x_min = xmin; x_max = xmax;
            y_min = ymin; y_max = ymax;
            z_min = zmin; z_max = zmax;
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
    // class for storting geometry information
    class DetectorGeometry
    {
    private:
        static DetectorGeometry * sInstance;
        static std::mutex sMutex;

    protected:
        DetectorGeometry(const std::string name);
        ~DetectorGeometry() {}
        std::string sName;

    public:
        // this singleton cannot be cloned
        DetectorGeometry(DetectorGeometry &other) = delete;
        // singleton should also not be assignable
        void operator=(const DetectorGeometry &) = delete;

        // static method that controls access to the singleton
        // instance
        static DetectorGeometry *getInstance(const std::string& name);

        std::string Name() const {
            return sName;
        }

        std::string GetWorldName()      { return sWorldName; }
        BoundingBox GetWorldBox()       { return sWorldBox; }
        std::string GetDetectorName()   { return sDetectorName; }
        BoundingBox GetDetectorBox()    { return sDetectorBox; }
        std::string GetCryostatName()   { return sCryostatName; }
        BoundingBox GetCryostatBox()    { return sCryostatBox; }
        int GetNumberOfTPCs()           { return sNumberOfTPCs; }

        std::vector<std::string> GetTPCNames()      { return sTPCNames; }
        std::vector<double> GetTPCMasses()          { return sTPCMasses; }
        std::vector<double> GetTPCDriftDistances()  { return sTPCDriftDistances; }
        BoundingBox GetTotalTPCBox()                { return sTotalTPCBox; }
        BoundingBox GetTotalActiveTPCBox()          { return sTotalActiveTPCBox; }
        double GetTotalTPCMass()                    { return sTotalTPCMass; }

        std::string GetTPCName(const size_t i);      
        BoundingBox GetTPCBox(const size_t i);       
        BoundingBox GetActiveTPCBox(const size_t i);
        double GetTPCMass(const size_t i);     
        double GetTPCDriftDistance(const size_t i);
        
        // get volume information for a point
        DetectorVolume getVolume(std::vector<double> position);
        DetectorVolume getVolume(double x, double y, double z);

        // function for finding total tpc volumes
        void findTotalTPCBoxes();
        // fill the geometry ttree
        void FillTTree();
        
    private:
        art::ServiceHandle<geo::Geometry> sGeometryService;
        geo::GeometryCore const* sGeometryCore;
        art::ServiceHandle<art::TFileService> sTFileService;

        TTree *sGeometryTTree;
        size_t sTriggerOffset;
        
        // map from volume names to volume type
        std::map<std::string,VolumeType> sVolumeTypeMap;

        // world volume
        std::string sWorldName;
        BoundingBox sWorldBox;

        // detector volume
        std::string sDetectorName;
        BoundingBox sDetectorBox;

        // cryostat volume
        std::string sCryostatName;
        BoundingBox sCryostatBox;

        // tpc volumes
        int sNumberOfTPCs;
        std::vector<std::string> sTPCNames;
        std::vector<BoundingBox> sTPCBoxes;
        std::vector<BoundingBox> sActiveTPCBoxes;
        std::vector<double> sTPCMasses;
        std::vector<double> sTPCDriftDistances;

        // full tpc volume
        BoundingBox sTotalTPCBox;
        BoundingBox sTotalActiveTPCBox;
        double sTotalTPCMass;

        ////////////////////////////////////////////////
        // detector material variables
        ////////////////////////////////////////////////
        // we will need to ask Geant4 about material 
        // properties for the detector volume
        // at each point of interest.  This requires holding 
        // this information in a
        // TGeoMaterial object, which is part of ROOT.
        const TGeoMaterial *sMaterial;
        geo::Point_t sMaterialPOI;
    };
}