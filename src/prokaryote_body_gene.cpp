
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
    vector<double> elipse_value, int resolution, string position)
{
    this->function = function;
    this->distance = distance;
    this->number_points = number_points;
    this->elipse_value = elipse_value;
    this->resolution = resolution;
    this->position = position;
}

unique_ptr<MainBodyPartGene> ConvertEndToMain(unique_ptr<EndBodyPartGene> gene) 
{
    vector<vector<double>> elipses_values;

    double adder_a = gene->elipse_value[0] / gene->number_points;
    double adder_b = gene->elipse_value[2] / gene->number_points;

    if (gene->position == "begin") {
        double run_distance_a = gene->elipse_value[0]-adder_a;
        double run_distance_b = gene->elipse_value[2]-adder_b;

        elipses_values.push_back({ gene->elipse_value[0] * 0.001,
                gene->elipse_value[1], gene->elipse_value[2] * 0.001 });
        
        for (int i = 0; i < gene->number_points; i++)
        {
            elipses_values.push_back({ sqrt(pow(gene->elipse_value[0], 2) - pow((run_distance_a),2)),
                gene->elipse_value[1], sqrt(pow(gene->elipse_value[2], 2) - pow((run_distance_b),2)) });
            run_distance_a -= adder_a;
            run_distance_b -= adder_b;
        }
    }
    else {
        double run_distance_a = 0;
        double run_distance_b = 0;

        for (int i = 0; i < gene->number_points; i++)
        {
            elipses_values.push_back({ sqrt(pow(gene->elipse_value[0], 2) - pow((run_distance_a),2)),
                gene->elipse_value[1], sqrt(pow(gene->elipse_value[2], 2) - pow((run_distance_b),2)) });
            run_distance_a += adder_a;
            run_distance_b += adder_b;
        }
        elipses_values.push_back({ gene->elipse_value[0] * 0.0001,
                gene->elipse_value[1], gene->elipse_value[2] * 0.0001 });
    }
    
    return make_unique<MainBodyPartGene>(MainBodyPartGene(gene->function, gene->distance, gene->number_points + 1, elipses_values, gene->resolution));
}