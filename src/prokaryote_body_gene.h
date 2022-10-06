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
	//distance between individual points on the function
	double distance;
	//Number of elipses to be present
	int number_points;
	//Vector of elipses values, which are 1. length of a axis 2. angle of a axis from baseline 3. length of b axis
	vector<vector<double>> elipses_values;
	//Number of points per elipse, higher the number the longer generation takes, but bacterioa look "smoother"
	int resolution;

	MainBodyPartGene(shared_ptr<VectorFunction> function, double distance, int number_points,
	vector<vector<double>> elipses_values, int resolution);
	MainBodyPartGene(shared_ptr<VectorFunction> function, double distance, int number_points,
		vector<double> elipses_values, int resolution);

};

class EndBodyPartGene
{
public:
	//Function defining outline of body points
	shared_ptr<VectorFunction> function;
	//Distance between individual points on the function
	double distance;
	//Number of elipses to be present
	int number_points;
	//Elipse values, which are 1. length of a axis 2. angle of a axis from baseline 3. length of b axis
	//This is then duplicated by the number of points to create number of similar, but differing in axis length, elipses
	vector<double> elipse_value;
	//Number of points per elipse, higher the number the longer generation takes, but bacterioa look "smoother"
	int resolution;
	//Begin or end position
	string position;

	EndBodyPartGene(shared_ptr<VectorFunction> function, double distance, int number_points,
		vector<double> elipses_values, int resolution, string position);

};

unique_ptr<MainBodyPartGene> ConvertEndToMain(unique_ptr<EndBodyPartGene> gene);

#endif

