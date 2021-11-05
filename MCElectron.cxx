#include "MCElectron.h"

namespace neutron
{
    MCElectron::MCElectron()
    {
        // initialize TTrees
        fMCElectronTree = fTFileService->make<TTree>("MCElectron", "MCElectron");
        // initialize number of Electrons
        fNumberOfElectrons = 0;
    }
    MCElectron::~MCElectron()
    {}

    void MCElectron::addElectron(Int_t eventId, simb::MCParticle particle)
    {
        // add new key to the event/particle map
        fElectronMap[std::make_pair(eventId,particle.TrackId())] = fNumberOfElectrons;
        fElectronMapKeys.emplace_back(std::vector<Int_t>({eventId,particle.TrackId()}));
        fNumberOfElectrons++;
        // add MCParticle values
        fEventId.emplace_back(eventId);
        fElectronTrackId.emplace_back(particle.TrackId());
        fElectronParentId.emplace_back(particle.Mother());
        fElectronStatusCode.emplace_back(particle.StatusCode());
        fElectronNumberOfDaughters.emplace_back(particle.NumberDaughters());
        std::vector<Int_t> daughters = {};
        for (Int_t i = 0; i < particle.NumberDaughters(); i++)
        {
            daughters.emplace_back(particle.Daughter(i));
        }
        fElectronDaughters.emplace_back(daughters);
        fElectronNumberOfTrajectoryPoints.emplace_back(particle.NumberTrajectoryPoints());
        // store trajectory information
        std::vector<Double_t> T; std::vector<Double_t> X; std::vector<Double_t> Y; 
        std::vector<Double_t> Z; std::vector<Double_t> E; std::vector<Double_t> Px; 
        std::vector<Double_t> Py; std::vector<Double_t> Pz;
        Double_t totalDistance = 0;
        std::vector<Double_t> distances = {0};
        bool enteredActiveVolume = false;
        bool exitActiveVolume = false;
        Int_t enteredActiveVolumeTime = -1;
        Int_t exitActiveVolumeTime = -1;
        for (size_t i = 0; i < particle.NumberTrajectoryPoints(); i++)
        {
            T.emplace_back(particle.T(i));
            X.emplace_back(particle.Vx(i));
            Y.emplace_back(particle.Vy(i));
            Z.emplace_back(particle.Vz(i));
            E.emplace_back(particle.E(i));
            Px.emplace_back(particle.Px(i));
            Py.emplace_back(particle.Py(i));
            Pz.emplace_back(particle.Pz(i));
            if (i != 0) {
                distances.emplace_back(EuclideanDistance(X[i],Y[i],Z[i],X[i-1],Y[i-1],Z[i-1]));
                totalDistance += distances[i];
            }
            DetectorVolume currentVolume = fGeometry->getVolume(
                particle.Vx(i), particle.Vy(i), particle.Vz(i)
            );
            if (currentVolume.volume_type == 2) {
                if (enteredActiveVolume == false) {
                    enteredActiveVolume = true;
                    enteredActiveVolumeTime = i;
                }
            }
            else {
                if (enteredActiveVolume == true && exitActiveVolume == false) {
                    exitActiveVolume = true;
                    exitActiveVolumeTime = i;
                }
            }
        }
        fElectronEnteredActiveVolume.emplace_back(enteredActiveVolume);
        fElectronEnteredActiveVolumeTime.emplace_back(enteredActiveVolumeTime);
        fElectronExitActiveVolume.emplace_back(exitActiveVolume);
        fElectronExitActiveVolumeTime.emplace_back(exitActiveVolumeTime);
        fElectronT.emplace_back(T);
        fElectronX.emplace_back(X);
        fElectronY.emplace_back(Y);
        fElectronZ.emplace_back(Z);
        fElectronE.emplace_back(E);
        fElectronPx.emplace_back(Px);
        fElectronPy.emplace_back(Py);
        fElectronPz.emplace_back(Pz);
        fElectronTotalDistance.emplace_back(totalDistance);
        fElectronDistances.emplace_back(distances);
        fElectronTotalDisplacement.emplace_back(EuclideanDistance(X[0],Y[0],Z[0],particle.EndX(),
            particle.EndY(),particle.EndZ())
        );
        // store other information
        fElectronProcess.emplace_back(particle.Process());
        fElectronEndProcess.emplace_back(particle.EndProcess());        
        // get volume information
        DetectorVolume beginVolume = fGeometry->getVolume(particle.Vx(0), particle.Vy(0), particle.Vz(0));
        DetectorVolume endVolume = fGeometry->getVolume(particle.EndX(), particle.EndY(), particle.EndZ());
        fElectronVolumeTypeBegin.emplace_back(beginVolume.volume_type);
        fElectronVolumeTypeEnd.emplace_back(endVolume.volume_type);
        fElectronVolumeNameBegin.emplace_back(beginVolume.volume_name);
        fElectronVolumeNameEnd.emplace_back(endVolume.volume_name);
        fElectronMaterialNameBegin.emplace_back(beginVolume.material_name);
        fElectronMaterialNameEnd.emplace_back(endVolume.material_name);
        fElectronMaterialBegin.emplace_back(beginVolume.material);
        fElectronMaterialEnd.emplace_back(endVolume.material);
    }

    void MCElectron::FillTTree()
    {
        // fill each Electron
        // temporary container for individual Electron info
        typedef struct {
            std::vector<Int_t> map_keys;
            Int_t event_id;
            Int_t track_id;
            Int_t parent_id;
            Int_t status_code;
            Int_t number_of_daughters;
            std::vector<Int_t> daughters;
            Int_t number_of_trajectory_points;
            std::vector<Double_t> t;
            std::vector<Double_t> x;
            std::vector<Double_t> y;
            std::vector<Double_t> z;
            std::vector<Double_t> E;
            std::vector<Double_t> px;
            std::vector<Double_t> py;
            std::vector<Double_t> pz;
            std::string process;
            std::string end_process;
            Int_t volume_type_begin;
            Int_t volume_type_end;
            std::string volume_name_begin;
            std::string volume_name_end;
            std::string material_name_begin;
            std::string material_name_end;
            Double_t material_begin;
            Double_t material_end;
            std::vector<Double_t> distances;
            Double_t total_distance;
            Double_t total_displacement;
            Bool_t entered_active_volume;
            Bool_t exit_active_volume;
            Int_t entered_active_volume_time;
            Int_t exit_active_volume_time;
        } Electron;
        Electron mc_Electron;
        // set up the branch structure
        fMCElectronTree->Branch("map_keys", &mc_Electron.map_keys);
        fMCElectronTree->Branch("event_id", &mc_Electron.event_id);
        fMCElectronTree->Branch("track_id", &mc_Electron.track_id);
        fMCElectronTree->Branch("parent_id", &mc_Electron.parent_id);
        fMCElectronTree->Branch("status_code", &mc_Electron.status_code);
        fMCElectronTree->Branch("number_of_daughters", &mc_Electron.number_of_daughters);
        fMCElectronTree->Branch("daughters", &mc_Electron.daughters);
        fMCElectronTree->Branch("number_of_trajectory_points", &mc_Electron.number_of_trajectory_points);
        fMCElectronTree->Branch("t", &mc_Electron.t);
        fMCElectronTree->Branch("x", &mc_Electron.x);
        fMCElectronTree->Branch("y", &mc_Electron.y);
        fMCElectronTree->Branch("z", &mc_Electron.z);
        fMCElectronTree->Branch("E", &mc_Electron.E);
        fMCElectronTree->Branch("px", &mc_Electron.px);
        fMCElectronTree->Branch("py", &mc_Electron.py);
        fMCElectronTree->Branch("pz", &mc_Electron.pz);
        fMCElectronTree->Branch("process", &mc_Electron.process);
        fMCElectronTree->Branch("end_process", &mc_Electron.end_process);
        fMCElectronTree->Branch("volume_type_begin", &mc_Electron.volume_type_begin);
        fMCElectronTree->Branch("volume_type_end", &mc_Electron.volume_type_end);
        fMCElectronTree->Branch("volume_name_begin", &mc_Electron.volume_name_begin);
        fMCElectronTree->Branch("volume_name_end", &mc_Electron.volume_name_end);
        fMCElectronTree->Branch("material_name_begin", &mc_Electron.material_name_begin);
        fMCElectronTree->Branch("material_name_end", &mc_Electron.material_name_end);
        fMCElectronTree->Branch("material_begin", &mc_Electron.material_begin);
        fMCElectronTree->Branch("material_end", &mc_Electron.material_end);
        fMCElectronTree->Branch("distances", &mc_Electron.distances);
        fMCElectronTree->Branch("total_distance", &mc_Electron.total_distance);
        fMCElectronTree->Branch("total_displacement", &mc_Electron.total_displacement);
        fMCElectronTree->Branch("entered_active_volume", &mc_Electron.entered_active_volume);
        fMCElectronTree->Branch("entered_active_volume_time", &mc_Electron.entered_active_volume_time);
        fMCElectronTree->Branch("exit_active_volume", &mc_Electron.exit_active_volume);
        fMCElectronTree->Branch("exit_active_volume_time", &mc_Electron.exit_active_volume_time);
        // iterate over all Electrons
        for (size_t i = 0; i < fNumberOfElectrons; i++)
        {
            mc_Electron.map_keys = fElectronMapKeys[i];
            mc_Electron.event_id = fEventId[i];
            mc_Electron.track_id = fElectronTrackId[i];
            mc_Electron.parent_id = fElectronParentId[i];
            mc_Electron.status_code = fElectronStatusCode[i];
            mc_Electron.number_of_daughters = fElectronNumberOfDaughters[i];
            mc_Electron.daughters = fElectronDaughters[i];
            mc_Electron.number_of_trajectory_points = fElectronNumberOfTrajectoryPoints[i];
            mc_Electron.t = fElectronT[i];
            mc_Electron.x = fElectronX[i];
            mc_Electron.y = fElectronY[i];
            mc_Electron.z = fElectronZ[i];
            mc_Electron.E = fElectronE[i];
            mc_Electron.px = fElectronPx[i];
            mc_Electron.py = fElectronPy[i];
            mc_Electron.pz = fElectronPz[i];
            mc_Electron.process = fElectronProcess[i];
            mc_Electron.end_process = fElectronEndProcess[i];
            mc_Electron.volume_type_begin = fElectronVolumeTypeBegin[i];
            mc_Electron.volume_type_end = fElectronVolumeTypeEnd[i];
            mc_Electron.volume_name_begin = fElectronVolumeNameBegin[i];
            mc_Electron.volume_name_end = fElectronVolumeNameEnd[i];
            mc_Electron.material_name_begin = fElectronMaterialNameBegin[i];
            mc_Electron.material_name_end = fElectronMaterialNameEnd[i];
            mc_Electron.material_begin = fElectronMaterialBegin[i];
            mc_Electron.material_end = fElectronMaterialEnd[i];
            mc_Electron.distances = fElectronDistances[i];
            mc_Electron.total_distance = fElectronTotalDistance[i];
            mc_Electron.total_displacement = fElectronTotalDisplacement[i];
            mc_Electron.entered_active_volume = fElectronEnteredActiveVolume[i];
            mc_Electron.entered_active_volume_time = fElectronEnteredActiveVolumeTime[i];
            mc_Electron.exit_active_volume = fElectronExitActiveVolume[i];
            mc_Electron.exit_active_volume_time = fElectronExitActiveVolumeTime[i];
            fMCElectronTree->Fill();
        }
    }
}