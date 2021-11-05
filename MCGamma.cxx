#include "MCGamma.h"

namespace neutron
{
    MCGamma::MCGamma()
    {
        // initialize TTrees
        fMCGammaTree = fTFileService->make<TTree>("MCGamma", "MCGamma");
        // initialize number of Gammas
        fNumberOfGammas = 0;
    }
    MCGamma::~MCGamma()
    {}

    void MCGamma::addGamma(Int_t eventId, simb::MCParticle particle)
    {
        // add new key to the event/particle map
        fGammaMap[std::make_pair(eventId,particle.TrackId())] = fNumberOfGammas;
        fGammaMapKeys.emplace_back(std::vector<Int_t>({eventId,particle.TrackId()}));
        fNumberOfGammas++;
        // add new gamma entry
        fNumberOfCaptureGammas.emplace_back(0);
        fCaptureGammaInitialEnergy.emplace_back(std::vector<Double_t>());
        fCaptureGammaInitialX.emplace_back(std::vector<Double_t>());
        fCaptureGammaFinalX.emplace_back(std::vector<Double_t>());
        fCaptureGammaInitialY.emplace_back(std::vector<Double_t>());
        fCaptureGammaFinalY.emplace_back(std::vector<Double_t>());
        fCaptureGammaInitialZ.emplace_back(std::vector<Double_t>());
        fCaptureGammaFinalZ.emplace_back(std::vector<Double_t>());
        fCaptureGammaInitialPx.emplace_back(std::vector<Double_t>());
        fCaptureGammaInitialPy.emplace_back(std::vector<Double_t>());
        fCaptureGammaInitialPz.emplace_back(std::vector<Double_t>());
        // add MCParticle values
        fEventId.emplace_back(eventId);
        fGammaTrackId.emplace_back(particle.TrackId());
        fGammaParentId.emplace_back(particle.Mother());
        fGammaStatusCode.emplace_back(particle.StatusCode());
        fGammaNumberOfDaughters.emplace_back(particle.NumberDaughters());
        std::vector<Int_t> daughters = {};
        for (Int_t i = 0; i < particle.NumberDaughters(); i++)
        {
            daughters.emplace_back(particle.Daughter(i));
        }
        fGammaDaughters.emplace_back(daughters);
        fGammaNumberOfTrajectoryPoints.emplace_back(particle.NumberTrajectoryPoints());
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
        fGammaEnteredActiveVolume.emplace_back(enteredActiveVolume);
        fGammaEnteredActiveVolumeTime.emplace_back(enteredActiveVolumeTime);
        fGammaExitActiveVolume.emplace_back(exitActiveVolume);
        fGammaExitActiveVolumeTime.emplace_back(exitActiveVolumeTime);
        fGammaT.emplace_back(T);
        fGammaX.emplace_back(X);
        fGammaY.emplace_back(Y);
        fGammaZ.emplace_back(Z);
        fGammaE.emplace_back(E);
        fGammaPx.emplace_back(Px);
        fGammaPy.emplace_back(Py);
        fGammaPz.emplace_back(Pz);
        fGammaTotalDistance.emplace_back(totalDistance);
        fGammaDistances.emplace_back(distances);
        fGammaTotalDisplacement.emplace_back(EuclideanDistance(X[0],Y[0],Z[0],particle.EndX(),
            particle.EndY(),particle.EndZ())
        );
        // store other information
        fGammaProcess.emplace_back(particle.Process());
        fGammaEndProcess.emplace_back(particle.EndProcess());        
        // get volume information
        DetectorVolume beginVolume = fGeometry->getVolume(particle.Vx(0), particle.Vy(0), particle.Vz(0));
        DetectorVolume endVolume = fGeometry->getVolume(particle.EndX(), particle.EndY(), particle.EndZ());
        fGammaVolumeTypeBegin.emplace_back(beginVolume.volume_type);
        fGammaVolumeTypeEnd.emplace_back(endVolume.volume_type);
        fGammaVolumeNameBegin.emplace_back(beginVolume.volume_name);
        fGammaVolumeNameEnd.emplace_back(endVolume.volume_name);
        fGammaMaterialNameBegin.emplace_back(beginVolume.material_name);
        fGammaMaterialNameEnd.emplace_back(endVolume.material_name);
        fGammaMaterialBegin.emplace_back(beginVolume.material);
        fGammaMaterialEnd.emplace_back(endVolume.material);
    }

    void MCGamma::FillTTree()
    {
        // fill each Gamma
        // temporary container for individual Gamma info
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
        } Gamma;
        Gamma mc_Gamma;
        // set up the branch structure
        fMCGammaTree->Branch("map_keys", &mc_Gamma.map_keys);
        fMCGammaTree->Branch("event_id", &mc_Gamma.event_id);
        fMCGammaTree->Branch("track_id", &mc_Gamma.track_id);
        fMCGammaTree->Branch("parent_id", &mc_Gamma.parent_id);
        fMCGammaTree->Branch("status_code", &mc_Gamma.status_code);
        fMCGammaTree->Branch("number_of_daughters", &mc_Gamma.number_of_daughters);
        fMCGammaTree->Branch("daughters", &mc_Gamma.daughters);
        fMCGammaTree->Branch("number_of_trajectory_points", &mc_Gamma.number_of_trajectory_points);
        fMCGammaTree->Branch("t", &mc_Gamma.t);
        fMCGammaTree->Branch("x", &mc_Gamma.x);
        fMCGammaTree->Branch("y", &mc_Gamma.y);
        fMCGammaTree->Branch("z", &mc_Gamma.z);
        fMCGammaTree->Branch("E", &mc_Gamma.E);
        fMCGammaTree->Branch("px", &mc_Gamma.px);
        fMCGammaTree->Branch("py", &mc_Gamma.py);
        fMCGammaTree->Branch("pz", &mc_Gamma.pz);
        fMCGammaTree->Branch("process", &mc_Gamma.process);
        fMCGammaTree->Branch("end_process", &mc_Gamma.end_process);
        fMCGammaTree->Branch("volume_type_begin", &mc_Gamma.volume_type_begin);
        fMCGammaTree->Branch("volume_type_end", &mc_Gamma.volume_type_end);
        fMCGammaTree->Branch("volume_name_begin", &mc_Gamma.volume_name_begin);
        fMCGammaTree->Branch("volume_name_end", &mc_Gamma.volume_name_end);
        fMCGammaTree->Branch("material_name_begin", &mc_Gamma.material_name_begin);
        fMCGammaTree->Branch("material_name_end", &mc_Gamma.material_name_end);
        fMCGammaTree->Branch("material_begin", &mc_Gamma.material_begin);
        fMCGammaTree->Branch("material_end", &mc_Gamma.material_end);
        fMCGammaTree->Branch("distances", &mc_Gamma.distances);
        fMCGammaTree->Branch("total_distance", &mc_Gamma.total_distance);
        fMCGammaTree->Branch("total_displacement", &mc_Gamma.total_displacement);
        fMCGammaTree->Branch("entered_active_volume", &mc_Gamma.entered_active_volume);
        fMCGammaTree->Branch("entered_active_volume_time", &mc_Gamma.entered_active_volume_time);
        fMCGammaTree->Branch("exit_active_volume", &mc_Gamma.exit_active_volume);
        fMCGammaTree->Branch("exit_active_volume_time", &mc_Gamma.exit_active_volume_time);
        // iterate over all Gammas
        for (size_t i = 0; i < fNumberOfGammas; i++)
        {
            mc_Gamma.map_keys = fGammaMapKeys[i];
            mc_Gamma.event_id = fEventId[i];
            mc_Gamma.track_id = fGammaTrackId[i];
            mc_Gamma.parent_id = fGammaParentId[i];
            mc_Gamma.status_code = fGammaStatusCode[i];
            mc_Gamma.number_of_daughters = fGammaNumberOfDaughters[i];
            mc_Gamma.daughters = fGammaDaughters[i];
            mc_Gamma.number_of_trajectory_points = fGammaNumberOfTrajectoryPoints[i];
            mc_Gamma.t = fGammaT[i];
            mc_Gamma.x = fGammaX[i];
            mc_Gamma.y = fGammaY[i];
            mc_Gamma.z = fGammaZ[i];
            mc_Gamma.E = fGammaE[i];
            mc_Gamma.px = fGammaPx[i];
            mc_Gamma.py = fGammaPy[i];
            mc_Gamma.pz = fGammaPz[i];
            mc_Gamma.process = fGammaProcess[i];
            mc_Gamma.end_process = fGammaEndProcess[i];
            mc_Gamma.volume_type_begin = fGammaVolumeTypeBegin[i];
            mc_Gamma.volume_type_end = fGammaVolumeTypeEnd[i];
            mc_Gamma.volume_name_begin = fGammaVolumeNameBegin[i];
            mc_Gamma.volume_name_end = fGammaVolumeNameEnd[i];
            mc_Gamma.material_name_begin = fGammaMaterialNameBegin[i];
            mc_Gamma.material_name_end = fGammaMaterialNameEnd[i];
            mc_Gamma.material_begin = fGammaMaterialBegin[i];
            mc_Gamma.material_end = fGammaMaterialEnd[i];
            mc_Gamma.distances = fGammaDistances[i];
            mc_Gamma.total_distance = fGammaTotalDistance[i];
            mc_Gamma.total_displacement = fGammaTotalDisplacement[i];
            mc_Gamma.entered_active_volume = fGammaEnteredActiveVolume[i];
            mc_Gamma.entered_active_volume_time = fGammaEnteredActiveVolumeTime[i];
            mc_Gamma.exit_active_volume = fGammaExitActiveVolume[i];
            mc_Gamma.exit_active_volume_time = fGammaExitActiveVolumeTime[i];
            fMCGammaTree->Fill();
        }
    }
}