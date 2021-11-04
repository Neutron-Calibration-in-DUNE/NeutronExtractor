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
        fNeutronMap[std::make_pair(eventId,particle.TrackId())] = fNumberOfNeutrons;
        fNumberOfNeutrons++;
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
        std::vector<Double_t> T; std::vector<Double_t> X; std::vector<Double_t> Y; std::vector<Double_t> Z;
        std::vector<Double_t> E; std::vector<Double_t> Px; std::vector<Double_t> Py; std::vector<Double_t> Pz;
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
        }
        fNeutronT.emplace_back(T);
        fNeutronX.emplace_back(X);
        fNeutronY.emplace_back(Y);
        fNeutronZ.emplace_back(Z);
        fNeutronE.emplace_back(E);
        fNeutronPx.emplace_back(Px);
        fNeutronPy.emplace_back(Py);
        fNeutronPz.emplace_back(Pz);
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
            fNeutronInelastic.emplace_back(std::vector<Int_t>());
        }

        // if (particle.Mother() != 0)
        // {
        //     lineage.emplace_back(particle.Mother());
        //     Int_t index = getNeutronIndex(eventId, particle.Mother());
        //     Int_t mother = fNeutronParentId[index];
        //     while (mother != 0)
        //     {
        //         lineage.emplace_back(fNeutronTrackId[index]);
        //         index = getNeutronIndex(eventId, particle.Mother());
        //         mother = fNeutronParentId[index];
        //     }
        // }
        // fNeutronInelastic.emplace_back(lineage);
    }

    void MCNeutron::FillTTree()
    {
        // temporary container for individual neutron info
        typedef struct {
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
        } NEUTRON;
        NEUTRON mc_neutron;
        // set up the branch structure
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
        // iterate over all neutrons
        for (size_t i = 0; i < fNumberOfNeutrons; i++)
        {
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
            fMCNeutronTree->Fill();
        }
    }
}