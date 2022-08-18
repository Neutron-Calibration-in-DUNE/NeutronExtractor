/**
 * @file SingleNeutronCaptures.h
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2022-08-18
 */
#pragma once
#include "Core.h"

namespace neutron
{
    class SingleNeutronCaptures
    {
    public:
        SingleNeutronCaptures();
        ~SingleNeutronCaptures();

        void processEvent();
        
    private:
    };
}