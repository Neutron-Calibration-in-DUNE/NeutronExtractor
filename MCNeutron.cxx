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
        fNumberOfNeutrons++;
        fEventId.emplace_back(eventId);
        fNeutronTrackId.emplace_back(particle.TrackId());
        fNeutronParentId.emplace_back(particle.Mother());
        fNeutronNumberOfTrajectoryPoints.emplace_back(particle.NumberTrajectoryPoints());
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
    }

    void MCNeutron::FillTTree()
    {
        // temporary container for individual neutron info
        typedef struct {
            Int_t event_id;
            Int_t track_id;
            Int_t parent_id;
            Int_t number_of_trajectory_points;
            std::vector<Double_t> t;
            std::vector<Double_t> x;
            std::vector<Double_t> y;
            std::vector<Double_t> z;
            std::vector<Double_t> E;
            std::vector<Double_t> px;
            std::vector<Double_t> py;
            std::vector<Double_t> pz;
        } NEUTRON;
        NEUTRON mc_neutron;
        // set up the branch structure
        fMCNeutronTree->Branch("event_id", &mc_neutron.event_id);
        fMCNeutronTree->Branch("track_id", &mc_neutron.track_id);
        fMCNeutronTree->Branch("parent_id", &mc_neutron.parent_id);
        fMCNeutronTree->Branch("number_of_trajectory_points", &mc_neutron.number_of_trajectory_points);
        fMCNeutronTree->Branch("t", &mc_neutron.t);
        fMCNeutronTree->Branch("x", &mc_neutron.x);
        fMCNeutronTree->Branch("y", &mc_neutron.y);
        fMCNeutronTree->Branch("z", &mc_neutron.z);
        fMCNeutronTree->Branch("E", &mc_neutron.E);
        fMCNeutronTree->Branch("px", &mc_neutron.px);
        fMCNeutronTree->Branch("py", &mc_neutron.py);
        fMCNeutronTree->Branch("pz", &mc_neutron.pz);
        // iterate over all neutrons
        for (size_t i = 0; i < fNumberOfNeutrons; i++)
        {
            mc_neutron.event_id = fEventId[i];
            mc_neutron.track_id = fNeutronTrackId[i];
            mc_neutron.parent_id = fNeutronParentId[i];
            mc_neutron.number_of_trajectory_points = fNeutronNumberOfTrajectoryPoints[i];
            mc_neutron.t = fNeutronT[i];
            mc_neutron.x = fNeutronX[i];
            mc_neutron.y = fNeutronY[i];
            mc_neutron.z = fNeutronZ[i];
            mc_neutron.E = fNeutronE[i];
            mc_neutron.px = fNeutronPx[i];
            mc_neutron.py = fNeutronPy[i];
            mc_neutron.pz = fNeutronPz[i];
            fMCNeutronTree->Fill();
        }
    }
}