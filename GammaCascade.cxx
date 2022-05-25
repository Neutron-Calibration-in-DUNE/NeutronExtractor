/**
 * @file GammaCascade.cxx
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @author Luca Pagani [lpagani@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2022-05-25
 */
#include "GammaCascade.h"

namespace neutron
{
    Gamma::Gamma(){}
    Gamma::~Gamma(){}
    Gamma::Gamma(double energy, double rel_intensity, 
            double tot_intensity, double conv_coefficient)
    : energy(energy), rel_intensity(rel_intensity)
    , tot_intensity(tot_intensity), conv_coefficient(conv_coefficient) 
    {}

    void Level::addGamma(double energy, double rel_intensity, 
            double tot_intensity, double conv_coefficient)
    {
        gammas.emplace_back(Gamma(energy, rel_intensity, tot_intensity, conv_coefficient));
    }
    Level::Level(){}
    Level::~Level(){}
    Level::Level(double energy, std::string spin_parity,
            std::string half_life)
    : energy(energy), spin_parity(spin_parity)
    , half_life(half_life)
    {}
    void Level::constructProbabilities()
    {
        double levelSum = 0;
        for (size_t i = 0; i < gammas.size(); i++) {
            levelSum += gammas[i].rel_intensity;
        }
        for (size_t i = 0; i < gammas.size(); i++) {
            probabilities.emplace_back(gammas[i].rel_intensity/levelSum);
            if (i == 0) {
                cumulative_probabilities.emplace_back(probabilities[i]);
            }
            else {
                cumulative_probabilities.emplace_back(
                    probabilities[i] + cumulative_probabilities[i-1]
                );
            }
        }
    }

    Event::Event(){}
    Event::~Event(){}
    Event::Event(std::vector<double> energies)
    : energies(std::move(energies))
    {}

    ENDFParser::ENDFParser(){}
    ENDFParser::~ENDFParser(){}
    ENDFParser::ENDFParser(std::string inputFile)
    : fInputFile(inputFile)
    {
        parse();
        constructProbabilities();
        constructNextLevels();
        print();
    }
    void ENDFParser::parse()
    {
        char reader[80];
        char levelEnergy[10];
        char gammaEnergy[10];
        char spinParity[18];
        char halfLife[10];
        char relativePhotonIntensity[8];
        char conversionCoefficient[7];
        char totalTransitionIntensity[10];

        fFile = fopen(fInputFile.c_str(), "r");

        if (!fFile){
            std::cout << "ERROR: Input file " << fInputFile << " not found!" << std::endl;
            exit(0);
        }

        for (int i = 0; i < 8; i++){
            relativePhotonIntensity[i] = 0;
        }
        for (int i = 0; i < 7; i++){
            conversionCoefficient[i] = 0;
        }
        for (int i = 0; i < 10; i++)
        {
            levelEnergy[i] = 0;
            gammaEnergy[i] = 0;
            halfLife[i] = 0;
            totalTransitionIntensity[i] = 0;
        }
        for (int i = 0; i < 18; i++)
        {
            spinParity[i] = 0;
        }
        while (fgets(reader, 80, fFile)!=NULL)
        {
            if (gammaEnergy[0] !=0)
            {
                double lev_energy = strtod(levelEnergy, NULL);

                size_t l = 0;
                bool found_match = false;
                for (; l < fLevels.size(); ++l) 
                {
                    const Level& lev = fLevels.at(l);
                    if ( lev_energy == lev.energy ) 
                    {
                        found_match = true;
                        break;
                    }
                }

                if (found_match) {
                    // Add a new Gamma to the existing level
                    fLevels.at(l).addGamma(strtod(gammaEnergy, NULL),
                        strtod(relativePhotonIntensity, NULL),
                        strtod(totalTransitionIntensity, NULL),
                        strtod(conversionCoefficient, NULL)
                    );
                }
                else 
                {
                    // Create a new Level and add the first Gamma to it
                    fLevels.emplace_back(Level(strtod(levelEnergy, NULL),
                        std::string(spinParity, 18),
                        std::string(halfLife, 10)) );

                    fLevels.back().addGamma(strtod(gammaEnergy, NULL),
                        strtod(relativePhotonIntensity, NULL),
                        strtod(totalTransitionIntensity, NULL),
                        strtod(conversionCoefficient, NULL)
                    );
                }

                for (int i = 0; i < 10; i++) {
                    gammaEnergy[i] = 0;
                }
            }
            else
            {
                if ((reader[5] == ' ')&&(reader[6] == ' ')&&(reader[7] == 'L'))
                {
                    strncpy(levelEnergy, &reader[9], 9);
                    strncpy(spinParity, &reader[21], 17);
                    strncpy(halfLife, &reader[39], 9);
                }
                if ((reader[5] == ' ')&&(reader[6] == ' ')&&(reader[7] == 'G'))
                {
                    strncpy(gammaEnergy, &reader[9], 9);
                    strncpy(relativePhotonIntensity, &reader[21], 7);
                    strncpy(totalTransitionIntensity, &reader[64], 9);

                    if(reader[55] == ' ')
                    {
                        for(int i = 1; i < 7; i++)
                        {
                            conversionCoefficient[0] = '0';
                            conversionCoefficient[i] = 0;
                        }
                    }
                    else
                    {
                        strncpy(conversionCoefficient, &reader[55], 6);
                    }
                }
            }
        }
        fclose (fFile);
    }
    const Level ENDFParser::operator[](size_t index)
    {
        if (index < fLevels.size()) {
            return fLevels[index];
        }
        else {
            return fLevels[0];
        }
    }
    const Level ENDFParser::max() {
        return fLevels[fLevels.size() - 1];
    }
    const Level ENDFParser::min() {
        return fLevels[0];
    }
    void ENDFParser::constructProbabilities()
    {
        for (size_t i = 0; i < fLevels.size(); i++)
        {
            fLevels[i].constructProbabilities();
        }
    }
    void ENDFParser::constructNextLevels()
    {
        for (size_t l_idx = 0; l_idx < fLevels.size(); l_idx++) {
            for (size_t g_idx = 0; g_idx < fLevels[l_idx].gammas.size(); g_idx++) {
                size_t match_idx = 0; 
                double match_diff = std::numeric_limits<double>::max();
                for (size_t t_idx = 0; t_idx < fLevels.size(); t_idx++) {
                    double temp_diff = fabs(fLevels[l_idx].energy - 
                        fLevels[l_idx].gammas[g_idx].energy - 
                        fLevels[t_idx].energy);
                    if (temp_diff < match_diff) 
                    {
                        match_idx = t_idx;
                        match_diff = temp_diff;
                    }
                }
                fLevels[l_idx].next_level.emplace_back(match_idx);
            }
        }
    }
    void ENDFParser::print()
    {
        for (size_t level = 0; level < fLevels.size(); level++)
        {
            std::cout << "Level " << level << " - " << fLevels[level].energy << " MeV\n";
            for (size_t gamma = 0; gamma < fLevels[level].gammas.size(); gamma++)
            {
                std::cout << "  - " << gamma << ") ";
                std::cout << fLevels[level].gammas[gamma].energy << " --> level ";
                std::cout << fLevels[level].next_level[gamma] << " - p(l-";
                std::cout << fLevels[level].next_level[gamma] << "|g-" << gamma << ") = ";
                std::cout << fLevels[level].probabilities[gamma] << "\n"; 
            }
        }
    }

    Cascade::Cascade()
    : fGenerator(fDevice())
    {}
    Cascade::~Cascade(){}
    Cascade::Cascade(std::string inputFile,
        size_t numberOfEvents,
        std::string outputFile)
    : fParser(ENDFParser(inputFile))
    , fNumberOfEvents(numberOfEvents)
    , fOutputFile(outputFile)
    , fGenerator(fDevice())
    {
        generate(fNumberOfEvents);
    }
    
    void Cascade::generate(size_t numberOfEvents)
    {
        for (auto i = 0; i < numberOfEvents; i++) {
            fEvents.emplace_back(generate());
        }
    }

    Event Cascade::generate()
    {
        std::vector<double> energies;
        double startingEnergy = fParser.max().energy;
        size_t mc_level_index = fParser.numLevels() - 1;

        while (startingEnergy >= fParser.min().energy)
        {
            double mc_random = fUniform(fGenerator);
            for (size_t mc_gamma_index = 0; 
                 mc_gamma_index < fParser[mc_level_index].gammas.size(); 
                 mc_gamma_index++
            )
            {
                if (mc_random <= fParser[mc_level_index].cumulative_probabilities[mc_gamma_index])
                {
                    energies.emplace_back(fParser[mc_level_index].gammas[mc_gamma_index].energy);
                    mc_level_index = fParser[mc_level_index].next_level[mc_gamma_index];
                    startingEnergy -= energies.back();
                    break;
                }
            }
        }
        return Event(energies);
    }

    void Cascade::save()
    {
        std::ofstream output;
        output.open(fOutputFile);
        for (size_t i = 0; i < fEvents.size(); i++)
        {
            for (size_t j = 0; j < fEvents[i].energies.size(); j++)
            {
                output << fEvents[i].energies[j];
                if (j < fEvents[i].energies.size()-1) {
                    output << ",";
                }
            }
            output << "\n";
        }
        output.close();
    }

    void Cascade::save_levels(std::string output_file)
    {
        std::ofstream output;
        output.open(output_file);
        for (size_t i = 0; i < fParser.numLevels(); i++)
        {
            for (size_t j = 0; j < fParser[i].gammas.size(); j++)
            {
                output << i << "," << fParser[i].energy << "," << j << "," << fParser[i].gammas[j].energy;
                output << "," << fParser[i].gammas[j].rel_intensity << ",";
                output << fParser[i].probabilities[j] << "," << fParser[i].cumulative_probabilities[j];
                output << "," << fParser[i].next_level[j] << "\n";
            }
        }
        output.close();
    }
}