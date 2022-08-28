
#include <future>
#include <thread>

#include "prokaryote_body_creator.h"

#define _USE_MATH_DEFINES

# define M_PI           3.14159265358979323846  /* pi */

std::mutex mtx;

//Class to hold info of elipse
Elipse::Elipse(vector<vector<double>> poi, vector<vector<int>> lin, vector<double> a, vector<double> b, vector<double> c, double a_angle, int n_angles)
{
    points = poi;
    lines = lin;
    this->a = a;
    this->b = b;
    this->c = c;
    this->a_angle = a_angle;
    this->n_angles = n_angles;
}


ProkaryotePartContainer::ProkaryotePartContainer(shared_ptr<VectorFunction> function, vector<unique_ptr<Elipse>> elipses, double origin_point, vector<double> origin, double end_point, vector<double> end)
{
    this->function = function;
    this->elipses = move(elipses);
    this->origin_point = origin_point;
    this->origin = origin;
    this->end_point = end_point;
    this->end = end;
}

ProkaryotePartContainer::ProkaryotePartContainer() {};

ProkaryoteBodyContainer::ProkaryoteBodyContainer(ProkaryotePartContainer part)
{
    this->elipses = move(part.elipses);
    this->origin_point = part.origin_point;
    this->end_point = part.end_point;
    this->origin = part.origin;
    this->end = part.end;
}

double FindEndPoint(shared_ptr<VectorFunction> function, double length, double a, int n_iterations, double threshold)
{
    // Using the Secant method to estimate root of distance function
    double gamma = 0.0001;
    std::lock_guard<std::mutex> lck(mtx);
    //Checking if given arguments are valid
    if (function->function_size != 3) { throw invalid_argument("Function intended to be for 3D space, need to have functions present."); }

    double n0 = a;
    double n1 = n0 + gamma;
    double n2 = 0;
    
    for (int i = 0; i <= n_iterations; i++)
    {
        n2 = ((n0) * (function->SizeFunctionValue( a, n1) - length) - (n1) * (function->SizeFunctionValue( a, n0)
            - length)) / (((function->SizeFunctionValue( a, n1) - length) - (function->SizeFunctionValue( a, n0) - length)));
        if ((abs(function->SizeFunctionValue( a, n2) - length)) < threshold) { return n2; }
        else { n0 = n1;  n1 = n2; }
    }
    
    return n2;
}


vector<double> VectorOnPlane(vector<double> d, vector<double> a, double length1, double length2, double angle)
{
    double f0 =  cos(-angle) * length1 * length2; double f1 = f0 * d[0];
    double f2 = a[0] * d[2] - d[0] * a[2];
    double f3 = -a[0] * d[1] + a[1] * d[0]; double f4 = -d[1] * f1; double f5 = - d[1] * f2 - f3 * d[2];
    double f6 = f3 * d[0];
    double f7 = pow(f3, 2) * pow(f5, 2) + pow(f6, 2) * pow(f2, 2) + pow(f6, 2) * pow(f3, 2);
    double f8 = 2 * pow(f3, 2) * f4 * f5 + 2 * pow(f6, 2) * f1 * f2;
    double f9 = pow(f3, 2) * pow(f4, 2) + pow(f6, 2) * pow(f1, 2) - pow(f6, 2) * pow(f3, 2) * pow(length1, 2);
    double k3 = (-f8 + sqrt(pow(f8, 2) - 4 * f7 * f9)) / (2 * f7); double k2 = (f1 + k3 * f2) / f3; double k1 = (f4 + k3 * f5) / f6;

    return { k1, k2, k3 };
}

//Function to get the separated points on function
vector<double> GetSeparators(shared_ptr<VectorFunction> function, double length, double starter_point, int n_points)
{
    //Take the starter points and runner points
    vector<double> return_points;
    double current_point = starter_point;
    double residual_point;
    //For requested number of points
    for (int i = 0; i <= n_points; i++)
    {
        //Gets the endpoint using the size function
        residual_point = FindEndPoint(function, length, current_point);
        return_points.emplace_back(residual_point);
        //Update the current point
        current_point = residual_point;
    }
    return return_points;
}

Elipse ElipsePoints(shared_ptr<VectorFunction> function, double variable_value, double a_length, double a_angle, double b_length, vector<double> c, int number_angles)
{
    //Define the basal vector for a elipse part creation
    map<string, double> values = { {function->variable_name, variable_value} };
    double theta = 0;
    double limit = 1e-10;
    
    vector<double> d = { function->functions[0].PartialDerivative(values, function->variable_name)+limit, function->functions[1].PartialDerivative(values, function->variable_name)+limit,
    function->functions[2].PartialDerivative(values, function->variable_name) +limit };
    
    double base1 = 1; double base2 = (-(base1*d[0] + d[2]*c[2])) / d[1]; double base3 = c[2];
    double base_length = sqrt(pow(base1, 2) + pow(base2, 2) + pow(base3, 2));

    double f1 = d[1] * base3 - d[2] * base2; double f2 = base1 * d[2] - d[0] * base3; double f3 = d[0] * base2 - d[1] * base1;
    
    vector<double> a = VectorOnPlane(d, { base1, base2, base3 }, a_length, base_length, a_angle);
    vector<double> b = VectorOnPlane(d, a, b_length, a_length, M_PI/2);
    vector<double> run_vector = a;
    
    vector<vector<double>> points_begin;
    vector<vector<double>> points_end;
    vector<vector<double>> points;
    vector<vector<int>> lines;

    double x_line; double y_line; double z_line; double t2; double t1;
    double k1; double k2; double k3; double k4; double k5; double k6;
    double x1; double y1; double z1;
    double x2; double y2; double z2;

    for ( int i = 0; i < number_angles; i++)
    {
        x_line = run_vector[0]; y_line = run_vector[1]; z_line = run_vector[2];
        k1 = y_line * a[0] - x_line * a[1]; k2 = y_line * b[0] - x_line * b[1];
        k3 = z_line * a[0] - x_line * a[2]; k4 = z_line * b[0] - x_line * b[2];
        t2 = atan((y_line * a[0] - x_line * a[1]) / (x_line * b[1] - y_line * b[0]));
        t1 = (a[0] * cos(t2) + b[0] * sin(t2)) / x_line;

        x1 = c[0] + t1 * x_line; y1 = c[1] + t1 * y_line; z1 = c[2] + t1 * z_line;
        x2 = c[0] - t1 * x_line; y2 = c[1] - t1 * y_line; z2 = c[2] - t1 * z_line;

        if (f1 * (x1 - c[0]) + f2 * (y1 - c[1]) + f3 * (z1 - c[2]) > 0) { points_begin.push_back({ x1, y1, z1 }); }
        else (points_end.push_back({ x1, y1, z1 }));
        if (f1 * (x2 - c[0]) + f2 * (y2 - c[1]) + f3 * (z2 - c[2]) > 0) { points_begin.push_back({ x2, y2, z2 }); }
        else (points_end.push_back({ x2, y2, z2 }));

        theta += ((1.0 / (number_angles)) * M_PI);
        
        run_vector = VectorOnPlane(d, a, 1, a_length, theta);
    }

    points.reserve(points_begin.size() + points_end.size()); // preallocate memory
    points.insert(points.end(), points_begin.begin(), points_begin.end());
    points.insert(points.end(), points_end.begin(), points_end.end());

    for (int i = 0; i <= points.size() - 2; i++) { lines.push_back({i, i + 1}); }
    lines.push_back({ (int)points.size() - 1, 0 });
    //lines.push_back({ 2*number_angles - 1, 0 });

    Elipse elipse = Elipse(points, lines, a, b, c, a_angle, number_angles);
    return elipse;
}

ProkaryotePartContainer AssembleContainerBody(unique_ptr<MainBodyPartGene> gene) 
{
    //Initilize container
    ProkaryotePartContainer container;
    container.origin_point = 0;
    container.origin = gene->function->get(0);
    
    double current_point = 0;
    vector<std::future<Elipse>> futures;
    
    for (int i = 0; i < gene->number_points; i++) 
    {
        current_point = FindEndPoint(gene->function, gene->distance, current_point);
        futures.emplace_back(std::async(std::launch::deferred, ElipsePoints, gene->function, current_point, gene->elipses_values[i][0],
            gene->elipses_values[i][1], gene->elipses_values[i][2], gene->function->get(current_point), gene->resolution));
    }
    
    for (const std::future<Elipse>& ftr : futures)
    {
        ftr.wait();
        Elipse some = ftr._Get_value();
        container.elipses.emplace_back(make_unique<Elipse>(some));
    }

    current_point = FindEndPoint(gene->function, gene->distance, current_point);
    container.end = gene->function->get(current_point);
    container.end_point = current_point;

    container.function = gene->function;
    
    return container;
}

double Length(double a1, double a2, double a3) 
{
    return sqrt(pow(a1, 2) + pow(a2, 2) + pow(a3, 2));
}

vector<unique_ptr<Elipse>> GetElipseConnect(ProkaryotePartContainer pp1, ProkaryotePartContainer pp2, double t1, double t2, double t3) 
{
    double a_length = (Length( pp1.elipses.back()->a[0], pp1.elipses.back()->a[1], pp1.elipses.back()->a[2])+ Length(pp2.elipses[0]->a[0], pp2.elipses[0]->a[1], pp2.elipses[0]->a[2])) / 2;
    double b_length = (Length(pp1.elipses.back()->b[0], pp1.elipses.back()->b[1], pp1.elipses.back()->b[2]) + Length(pp2.elipses[0]->b[0], pp2.elipses[0]->b[1], pp2.elipses[0]->b[2])) / 2;
    int n_angles = (int) (pp1.elipses.back()->n_angles + pp2.elipses[0]->n_angles) / 2;
    double a_angle = (pp1.elipses.back()->a_angle + pp2.elipses[0]->a_angle) / 2;
    
    Elipse elipse1 = ElipsePoints(pp1.function, pp2.origin_point, a_length, a_angle, b_length, pp2.origin, n_angles);
    Elipse elipse2 = ElipsePoints(pp2.function, pp2.origin_point, a_length, a_angle, b_length, pp2.origin, n_angles);

    vector<vector<double>> points;

    for (int i = 0; i < elipse1.points.size(); i++)
    {
        points.push_back({ ((elipse1.points[i][0] + elipse2.points[i][0]) / 2) + t1, ((elipse1.points[i][1] + elipse2.points[i][1]) / 2) + t2 ,
            ((elipse1.points[i][2] + elipse2.points[i][2]) / 2) + t3 });
    }

    unique_ptr<Elipse> inter_elipse = make_unique<Elipse>(Elipse(points, elipse1.lines, elipse1.a, elipse1.b, pp1.end, elipse1.a_angle, n_angles));

    vector<unique_ptr<Elipse>> ConnectedElipses;

    for (int i = 0; i < pp1.elipses.size(); i++)
    {
        ConnectedElipses.emplace_back(move(pp1.elipses[i]));
    }

    ConnectedElipses.emplace_back(move(inter_elipse));

    for (int i = 0; i < pp2.elipses.size(); i++)
    {
        ConnectedElipses.emplace_back(move(pp2.elipses[i]));
    }

    return ConnectedElipses;
}

ProkaryotePartContainer ConnectParts(ProkaryotePartContainer prokaryote_part_1, ProkaryotePartContainer prokaryote_part_2, int counter)
{
    
    //Get move transformators to use
    double t1 = (prokaryote_part_1.end[0]) - (prokaryote_part_2.origin[0]);
    double t2 = (prokaryote_part_1.end[1]) - (prokaryote_part_2.origin[1]);
    double t3 = (prokaryote_part_1.end[2]) - (prokaryote_part_2.origin[2]);
    
    //Change every position in second part
    for (int i = 0; i < prokaryote_part_2.elipses.size(); i++)
    {
        for (int j = 0; j < prokaryote_part_2.elipses[i]->points.size(); j++)
        {
            prokaryote_part_2.elipses[i]->points[j][0] += t1;
            prokaryote_part_2.elipses[i]->points[j][1] += t2;
            prokaryote_part_2.elipses[i]->points[j][2] += t3;
        }
    }
    
    vector<double> end = { prokaryote_part_2.end[0] + t1, prokaryote_part_2.end[1] + t2, prokaryote_part_2.end[2] + t3 };
    
    shared_ptr<VectorFunction> function = prokaryote_part_2.function;
    double origin_point = prokaryote_part_1.origin_point;
    vector<double> origin = prokaryote_part_1.origin;
    double end_point = prokaryote_part_2.end_point;
    
    vector<unique_ptr<Elipse>> ConnectedElipses = GetElipseConnect(move(prokaryote_part_1), move(prokaryote_part_2), t1, t2, t3);

    ProkaryotePartContainer ConnectedPart = ProkaryotePartContainer(function, move(ConnectedElipses), origin_point, origin, end_point, end);
    
    return ConnectedPart;
}

/*In future, add elipses connecting parts. Should work on principle of two planes and plane coming out of their intersection, the same principle.*/

ProkaryoteBodyContainer AssembleParts(vector<ProkaryotePartContainer> prokaryote_parts)
{
    if (prokaryote_parts.size() == 1) { ProkaryoteBodyContainer body = ProkaryoteBodyContainer(move(prokaryote_parts[0])); return body; }
    else if (prokaryote_parts.size() == 2) { 
        ProkaryotePartContainer running_container = ConnectParts(move(prokaryote_parts[0]), move(prokaryote_parts[1]), 0);
        ProkaryoteBodyContainer body = ProkaryoteBodyContainer(move(running_container)); 
        return body; 
    }
    else {
        
        int counter = prokaryote_parts[0].elipses.size() + prokaryote_parts[1].elipses.size();
        ProkaryotePartContainer running_container = ConnectParts(move(prokaryote_parts[0]), move(prokaryote_parts[1]), 0);
        
        for (int i = 2; i < prokaryote_parts.size(); i++)
        {
            running_container = ConnectParts(move(running_container), move(prokaryote_parts[i]), counter);
        }

        ProkaryoteBodyContainer body = ProkaryoteBodyContainer(move(running_container));
        return body;
    }
}


