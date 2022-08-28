
#include "prokaryote_body_gene.h"

//Class to hold info of elipse

MainBodyPartGene::MainBodyPartGene(shared_ptr<VectorFunction> function, double distance, int number_points,
    vector<vector<double>> elipses_values, int resolution)
{
    this->function = function;
    this->distance = distance;
    this->number_points = number_points;
    this->elipses_values = elipses_values;
    this->resolution = resolution;
}

EndBodyPartGene::EndBodyPartGene(shared_ptr<VectorFunction> function, double distance, int number_points,
    vector<double> elipse_value, int resolution, double multiplier, string position)
{
    this->function = function;
    this->distance = distance;
    this->number_points = number_points;
    this->elipse_value = elipse_value;
    this->resolution = resolution;
    this->multiplier = multiplier;
    this->position = position;
}

unique_ptr<MainBodyPartGene> ConvertEndToMain(unique_ptr<EndBodyPartGene> gene) 
{
    vector<vector<double>> elipses_values;
    double adder = (1 - gene->multiplier) / (gene->number_points);
    
    if (gene->position == "begin") {
        double run_multiplier = 1;
        for (int j = 0; j < gene->number_points; j++) { run_multiplier *= (1 - adder * (j + 1));}

        elipses_values.push_back({ gene->elipse_value[0] * 0.0001,
                gene->elipse_value[1], gene->elipse_value[2] * 0.0001 });
        
        for (int i = 0; i < gene->number_points; i++)
        {
            elipses_values.push_back({ gene->elipse_value[0] * run_multiplier,
                gene->elipse_value[1], gene->elipse_value[2] * run_multiplier });
            run_multiplier = run_multiplier / (gene->multiplier + (adder * (i)));
        }
        elipses_values.push_back({ gene->elipse_value[0] * 1,
                gene->elipse_value[1], gene->elipse_value[2] * 1 });
    }
    else {
        double run_multiplier = (1 - adder);
        elipses_values.push_back({ gene->elipse_value[0] * 1,
                gene->elipse_value[1], gene->elipse_value[2] * 1 });
        for (int i = 0; i < gene->number_points; i++)
        {
            elipses_values.push_back({ gene->elipse_value[0] * run_multiplier,
                gene->elipse_value[1], gene->elipse_value[2] * run_multiplier });
            run_multiplier *= (1 - (i+2)*adder);
        }
        elipses_values.push_back({ gene->elipse_value[0] * 0.0001,
                gene->elipse_value[1], gene->elipse_value[2] * 0.0001 });
    }
    
    return make_unique<MainBodyPartGene>(MainBodyPartGene(gene->function, gene->distance, gene->number_points + 2, elipses_values, gene->resolution));
}