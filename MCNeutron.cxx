#include "MCNeutron.h"

namespace neutron
{
    MCNeutron::MCNeutron()
    {
        // initialize TTrees
        fMCNeutronTree = fTFileService->make<TTree>("MCNeutron", "MCNeutron");
        // initialize number of neutrons
        fNumberOfNeutrons = 0;
    }
    MCNeutron::~MCNeutron()
    {}

    void MCNeutron::addNeutron(Int_t eventId, simb::MCParticle particle)
    {
        // add new key to the event/particle map
        fNeutronMap[std::make_pair(eventId,particle.TrackId())] = fNumberOfNeutrons;
        fNeutronMapKeys.emplace_back(std::vector<Int_t>({eventId,particle.TrackId()}));
        fNumberOfNeutrons++;
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
        fNeutronTrackId.emplace_back(particle.TrackId());
        fNeutronParentId.emplace_back(particle.Mother());
        fNeutronStatusCode.emplace_back(particle.StatusCode());
        fNeutronNumberOfDaughters.emplace_back(particle.NumberDaughters());
        std::vector<Int_t> daughters = {};
        for (Int_t i = 0; i < particle.NumberDaughters(); i++)
        {
            daughters.emplace_back(particle.Daughter(i));
        }
        fNeutronDaughters.emplace_back(daughters);
        fNeutronNumberOfTrajectoryPoints.emplace_back(particle.NumberTrajectoryPoints());
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
        fNeutronEnteredActiveVolume.emplace_back(enteredActiveVolume);
        fNeutronEnteredActiveVolumeTime.emplace_back(enteredActiveVolumeTime);
        fNeutronExitActiveVolume.emplace_back(exitActiveVolume);
        fNeutronExitActiveVolumeTime.emplace_back(exitActiveVolumeTime);
        fNeutronT.emplace_back(T);
        fNeutronX.emplace_back(X);
        fNeutronY.emplace_back(Y);
        fNeutronZ.emplace_back(Z);
        fNeutronE.emplace_back(E);
        fNeutronPx.emplace_back(Px);
        fNeutronPy.emplace_back(Py);
        fNeutronPz.emplace_back(Pz);
        fNeutronTotalDistance.emplace_back(totalDistance);
        fNeutronDistances.emplace_back(distances);
        fNeutronTotalDisplacement.emplace_back(EuclideanDistance(X[0],Y[0],Z[0],particle.EndX(),
            particle.EndY(),particle.EndZ())
        );
        // store other information
        fNeutronProcess.emplace_back(particle.Process());
        fNeutronEndProcess.emplace_back(particle.EndProcess());
        // see if neutron is a primary or not and then store its lineage
        std::vector<Int_t> lineage = {particle.TrackId()};
        if (particle.EndProcess() == "nCapture" && particle.Mother() != 0)
        {
            fNeutronInelastic.emplace_back(lineage);
            Int_t index = getNeutronIndex(eventId, particle.Mother());
            lineage.emplace_back(fNeutronTrackId[index]);
            fNeutronInelastic[index] = lineage;
            Int_t mother = fNeutronParentId[index];
            while (mother != 0)
            {
                index = getNeutronIndex(eventId, particle.Mother());
                lineage.emplace_back(fNeutronTrackId[index]);
                fNeutronInelastic[index] = lineage;
                mother = fNeutronParentId[index];
            }
        }
        else
        {
            fNeutronInelastic.emplace_back(lineage);
        }
        // get volume information
        DetectorVolume beginVolume = fGeometry->getVolume(particle.Vx(0), particle.Vy(0), particle.Vz(0));
        DetectorVolume endVolume = fGeometry->getVolume(particle.EndX(), particle.EndY(), particle.EndZ());
        fNeutronVolumeTypeBegin.emplace_back(beginVolume.volume_type);
        fNeutronVolumeTypeEnd.emplace_back(endVolume.volume_type);
        fNeutronVolumeNameBegin.emplace_back(beginVolume.volume_name);
        fNeutronVolumeNameEnd.emplace_back(endVolume.volume_name);
        fNeutronMaterialNameBegin.emplace_back(beginVolume.material_name);
        fNeutronMaterialNameEnd.emplace_back(endVolume.material_name);
        fNeutronMaterialBegin.emplace_back(beginVolume.material);
        fNeutronMaterialEnd.emplace_back(endVolume.material);
    }

    void MCNeutron::addGamma(Int_t eventId, simb::MCParticle gamma)
    {
        // get neutron index
        Int_t neutronIndex = fNeutronMap[std::make_pair(eventId,gamma.Mother())];
        // update number of capture gamms
        fNumberOfCaptureGammas[neutronIndex]++;
        // add gamma info
        fCaptureGammaInitialEnergy[neutronIndex].emplace_back(gamma.E(0));
        fCaptureGammaInitialX[neutronIndex].emplace_back(gamma.Vx(0));
        fCaptureGammaFinalX[neutronIndex].emplace_back(gamma.EndX());
        fCaptureGammaInitialY[neutronIndex].emplace_back(gamma.Vy(0));
        fCaptureGammaFinalY[neutronIndex].emplace_back(gamma.EndY());
        fCaptureGammaInitialZ[neutronIndex].emplace_back(gamma.Vz(0));
        fCaptureGammaFinalZ[neutronIndex].emplace_back(gamma.EndZ());
        fCaptureGammaInitialPx[neutronIndex].emplace_back(gamma.Px(0));
        fCaptureGammaInitialPy[neutronIndex].emplace_back(gamma.Py(0));
        fCaptureGammaInitialPz[neutronIndex].emplace_back(gamma.Pz(0));
    }

    void MCNeutron::FillTTree()
    {
        // fill each neutron
        // temporary container for individual neutron info
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
            std::vector<Int_t> inelastic;
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
            Int_t number_of_capture_gammas;
            std::vector<Double_t> capture_gamma_initial_energy;
            std::vector<Double_t> capture_gamma_initial_x;
            std::vector<Double_t> capture_gamma_final_x;
            std::vector<Double_t> capture_gamma_initial_y;
            std::vector<Double_t> capture_gamma_final_y;
            std::vector<Double_t> capture_gamma_initial_z;
            std::vector<Double_t> capture_gamma_final_z;
            std::vector<Double_t> capture_gamma_initial_px;
            std::vector<Double_t> capture_gamma_initial_py;
            std::vector<Double_t> capture_gamma_initial_pz;
        } NEUTRON;
        NEUTRON mc_neutron;
        // set up the branch structure
        fMCNeutronTree->Branch("map_keys", &mc_neutron.map_keys);
        fMCNeutronTree->Branch("event_id", &mc_neutron.event_id);
        fMCNeutronTree->Branch("track_id", &mc_neutron.track_id);
        fMCNeutronTree->Branch("parent_id", &mc_neutron.parent_id);
        fMCNeutronTree->Branch("status_code", &mc_neutron.status_code);
        fMCNeutronTree->Branch("number_of_daughters", &mc_neutron.number_of_daughters);
        fMCNeutronTree->Branch("daughters", &mc_neutron.daughters);
        fMCNeutronTree->Branch("number_of_trajectory_points", &mc_neutron.number_of_trajectory_points);
        fMCNeutronTree->Branch("t", &mc_neutron.t);
        fMCNeutronTree->Branch("x", &mc_neutron.x);
        fMCNeutronTree->Branch("y", &mc_neutron.y);
        fMCNeutronTree->Branch("z", &mc_neutron.z);
        fMCNeutronTree->Branch("E", &mc_neutron.E);
        fMCNeutronTree->Branch("px", &mc_neutron.px);
        fMCNeutronTree->Branch("py", &mc_neutron.py);
        fMCNeutronTree->Branch("pz", &mc_neutron.pz);
        fMCNeutronTree->Branch("process", &mc_neutron.process);
        fMCNeutronTree->Branch("end_process", &mc_neutron.end_process);
        fMCNeutronTree->Branch("inelastic", &mc_neutron.inelastic);
        fMCNeutronTree->Branch("volume_type_begin", &mc_neutron.volume_type_begin);
        fMCNeutronTree->Branch("volume_type_end", &mc_neutron.volume_type_end);
        fMCNeutronTree->Branch("volume_name_begin", &mc_neutron.volume_name_begin);
        fMCNeutronTree->Branch("volume_name_end", &mc_neutron.volume_name_end);
        fMCNeutronTree->Branch("material_name_begin", &mc_neutron.material_name_begin);
        fMCNeutronTree->Branch("material_name_end", &mc_neutron.material_name_end);
        fMCNeutronTree->Branch("material_begin", &mc_neutron.material_begin);
        fMCNeutronTree->Branch("material_end", &mc_neutron.material_end);
        fMCNeutronTree->Branch("distances", &mc_neutron.distances);
        fMCNeutronTree->Branch("total_distance", &mc_neutron.total_distance);
        fMCNeutronTree->Branch("total_displacement", &mc_neutron.total_displacement);
        fMCNeutronTree->Branch("entered_active_volume", &mc_neutron.entered_active_volume);
        fMCNeutronTree->Branch("entered_active_volume_time", &mc_neutron.entered_active_volume_time);
        fMCNeutronTree->Branch("exit_active_volume", &mc_neutron.exit_active_volume);
        fMCNeutronTree->Branch("exit_active_volume_time", &mc_neutron.exit_active_volume_time);
        fMCNeutronTree->Branch("number_of_capture_gammas", &mc_neutron.number_of_capture_gammas);
        fMCNeutronTree->Branch("capture_gamma_initial_energy", &mc_neutron.capture_gamma_initial_energy);
        fMCNeutronTree->Branch("capture_gamma_initial_x", &mc_neutron.capture_gamma_initial_x);
        fMCNeutronTree->Branch("capture_gamma_final_x", &mc_neutron.capture_gamma_final_x);
        fMCNeutronTree->Branch("capture_gamma_initial_y", &mc_neutron.capture_gamma_initial_y);
        fMCNeutronTree->Branch("capture_gamma_final_y", &mc_neutron.capture_gamma_final_y);
        fMCNeutronTree->Branch("capture_gamma_initial_z", &mc_neutron.capture_gamma_initial_z);
        fMCNeutronTree->Branch("capture_gamma_final_z", &mc_neutron.capture_gamma_final_z);
        fMCNeutronTree->Branch("capture_gamma_initial_px", &mc_neutron.capture_gamma_initial_px);
        fMCNeutronTree->Branch("capture_gamma_initial_py", &mc_neutron.capture_gamma_initial_py);
        fMCNeutronTree->Branch("capture_gamma_initial_pz", &mc_neutron.capture_gamma_initial_pz);
        // iterate over all neutrons
        for (size_t i = 0; i < fNumberOfNeutrons; i++)
        {
            mc_neutron.map_keys = fNeutronMapKeys[i];
            mc_neutron.event_id = fEventId[i];
            mc_neutron.track_id = fNeutronTrackId[i];
            mc_neutron.parent_id = fNeutronParentId[i];
            mc_neutron.status_code = fNeutronStatusCode[i];
            mc_neutron.number_of_daughters = fNeutronNumberOfDaughters[i];
            mc_neutron.daughters = fNeutronDaughters[i];
            mc_neutron.number_of_trajectory_points = fNeutronNumberOfTrajectoryPoints[i];
            mc_neutron.t = fNeutronT[i];
            mc_neutron.x = fNeutronX[i];
            mc_neutron.y = fNeutronY[i];
            mc_neutron.z = fNeutronZ[i];
            mc_neutron.E = fNeutronE[i];
            mc_neutron.px = fNeutronPx[i];
            mc_neutron.py = fNeutronPy[i];
            mc_neutron.pz = fNeutronPz[i];
            mc_neutron.process = fNeutronProcess[i];
            mc_neutron.end_process = fNeutronEndProcess[i];
            mc_neutron.inelastic = fNeutronInelastic[i];
            mc_neutron.volume_type_begin = fNeutronVolumeTypeBegin[i];
            mc_neutron.volume_type_end = fNeutronVolumeTypeEnd[i];
            mc_neutron.volume_name_begin = fNeutronVolumeNameBegin[i];
            mc_neutron.volume_name_end = fNeutronVolumeNameEnd[i];
            mc_neutron.material_name_begin = fNeutronMaterialNameBegin[i];
            mc_neutron.material_name_end = fNeutronMaterialNameEnd[i];
            mc_neutron.material_begin = fNeutronMaterialBegin[i];
            mc_neutron.material_end = fNeutronMaterialEnd[i];
            mc_neutron.distances = fNeutronDistances[i];
            mc_neutron.total_distance = fNeutronTotalDistance[i];
            mc_neutron.total_displacement = fNeutronTotalDisplacement[i];
            mc_neutron.entered_active_volume = fNeutronEnteredActiveVolume[i];
            mc_neutron.entered_active_volume_time = fNeutronEnteredActiveVolumeTime[i];
            mc_neutron.exit_active_volume = fNeutronExitActiveVolume[i];
            mc_neutron.exit_active_volume_time = fNeutronExitActiveVolumeTime[i];
            mc_neutron.number_of_capture_gammas = fNumberOfCaptureGammas[i];
            mc_neutron.capture_gamma_initial_energy = fCaptureGammaInitialEnergy[i];
            mc_neutron.capture_gamma_initial_x = fCaptureGammaInitialX[i];
            mc_neutron.capture_gamma_final_x = fCaptureGammaFinalX[i];
            mc_neutron.capture_gamma_initial_y = fCaptureGammaInitialY[i];
            mc_neutron.capture_gamma_final_y = fCaptureGammaFinalY[i];
            mc_neutron.capture_gamma_initial_z = fCaptureGammaInitialZ[i];
            mc_neutron.capture_gamma_final_z = fCaptureGammaFinalZ[i];
            mc_neutron.capture_gamma_initial_px = fCaptureGammaInitialPx[i];
            mc_neutron.capture_gamma_initial_py = fCaptureGammaInitialPy[i];
            mc_neutron.capture_gamma_initial_pz = fCaptureGammaInitialPz[i];
            fMCNeutronTree->Fill();
        }
    }
}