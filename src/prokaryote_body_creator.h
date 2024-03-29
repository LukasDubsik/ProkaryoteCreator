#ifndef PROKARYOTE_BODY_CREATOR_H
#define PROKARYOTE_BODY_CREATOR_H

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <math.h>
#include <mutex>

#include "prokaryote_body_gene.h"

#include "function.h"

using namespace std;

//Class containing points and lines of elipse
class Elipse
{
public:
    vector<double> a;
    vector<double> b;
    vector<double> c;
    double a_angle;
    int n_angles;

    vector<vector<double>> points;
    vector<vector<int>> lines;

    Elipse(vector<vector<double>> poi, vector<vector<int>> lin, vector<double> a, vector<double> b, vector<double> c, double a_angle, int n_angles);
};
//These classes may look the same, but are used for different purposes, so for differencing and future application
class ProkaryotePartContainer 
{
public:
    shared_ptr<VectorFunction> function;
    vector<unique_ptr<Elipse>> elipses;
    double origin_point;
    vector<double> origin;
    double end_point;
    vector<double> end;

    vector<double> derivative_end = { 0, 0, 0 };
    double end_angle = 0;

    ProkaryotePartContainer(shared_ptr<VectorFunction> function, vector<unique_ptr<Elipse>> elipses, double origin_point, vector<double> origin, double end_point, vector<double> end);
    ProkaryotePartContainer();
};

class ProkaryoteBodyContainer 
{
public:
    vector<unique_ptr<Elipse>> elipses;
    double origin_point;
    double end_point;
    vector<double> origin;
    vector<double> end;

    ProkaryoteBodyContainer(ProkaryotePartContainer part);
};

//Get points upon which create the elipses
double FindEndPoint(shared_ptr<VectorFunction> function, double length, double a, int n_iterations = 20, double threshold = 1e-3);
vector<double> GetSeparators(shared_ptr<VectorFunction> function, double length, double starter_point, int n_points);

ProkaryotePartContainer AssembleContainerBody(unique_ptr<MainBodyPartGene> gene); //Cap is created in the same way, just precise gene

ProkaryotePartContainer ConnectParts(ProkaryotePartContainer prokaryote_part_1, ProkaryotePartContainer prokaryote_part_2, int counter, int run);
ProkaryoteBodyContainer AssembleParts(vector<ProkaryotePartContainer> prokaryote_parts);

#endif
