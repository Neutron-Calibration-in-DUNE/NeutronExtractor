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
        // iterate over all neutrons
        for (size_t i = 0; i < fNumberOfNeutrons; i++)
        {
            mc_neutron.event_id = fEventId[i];
            mc_neutron.track_id = fNeutronTrackId[i];
            mc_neutron.number_of_trajectory_points = fNeutronNumberOfTrajectoryPoints[i];
            fMCNeutronTree->Fill();
        }
    }
}