/**
 * @file    Configuration.h
 * @author  Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief   A struct for holding LArSoft configuration parameters
 *          for the NeutronExtractor module.
 * @version 0.1
 * @date 2022-02-15
 */
#pragma once

// art framework includes
#include "art/Framework/Core/EDAnalyzer.h"
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

namespace neutron
{
    /**
     * @brief A collection of fhicl parameters for the NeutronExtractor
     * module.  Each of these must be specified, which have default values
     * in the accompanying fhicl file.
     */
    struct Configuration
    {
        fhicl::Atom<art::InputTag> LArGeantProducerLabel
        {
            fhicl::Name("LArGeantProducerLabel"),
            fhicl::Comment("Tag of the input data product for the largeant side of the simulation.")
        };
        fhicl::Atom<art::InputTag> IonAndScintProducerLabel
        {
            fhicl::Name("IonAndScintProducerLabel"),
            fhicl::Comment("Tag of the input data product for the IonAndScint side of the simulation.")
        };
        fhicl::Atom<art::InputTag> Cluster3DProducerLabel
        {
            fhicl::Name("Cluster3DProducerLabel"),
            fhicl::Comment("Tag of the input data product for the Cluster3D side of the simulation.")
        };

        fhicl::Atom<bool> FillNeutronCapture
        {
            fhicl::Name("FillNeutronCapture"),
            fhicl::Comment("Whether to fill the neutron captures reco product.")
        };
    };

    using Parameters = art::EDAnalyzer::Table<Configuration>;
}