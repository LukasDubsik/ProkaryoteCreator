#ifndef PROKARYOTE_BODY_GENE_H
#define PROKARYOTE_BODY_GENE_H

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <math.h>

#include "function.h"

using namespace std;

class MainBodyPartGene
{
public:
	//Function defining outline of body points
	shared_ptr<VectorFunction> function;
	//Start point of function and distance of points
	double distance;
	//Other important values
	int number_points;
	vector<vector<double>> elipses_values;
	int resolution;

	MainBodyPartGene(shared_ptr<VectorFunction> function, double distance, int number_points,
	vector<vector<double>> elipses_values, int resolution);

};

class EndBodyPartGene
{
public:
	//Function defining outline of body points
	shared_ptr<VectorFunction> function;
	//Start point of function and distance of points
	double distance;
	//Other important values
	int number_points;
	vector<double> elipse_value;
	int resolution;

	double multiplier;
	//Begin or end
	string position;

	EndBodyPartGene(shared_ptr<VectorFunction> function, double distance, int number_points,
		vector<double> elipses_values, int resolution, double multiplier, string position);

};

unique_ptr<MainBodyPartGene> ConvertEndToMain(unique_ptr<EndBodyPartGene> gene);

#endif

