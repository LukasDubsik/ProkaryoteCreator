
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

double Distance_here(vector<double> a, vector<double> b)
{
    return sqrt(pow((a[0] - b[0]), 2) + pow((a[1] - b[1]), 2) + pow((a[2] - b[2]), 2));
}

vector<double> DecideBase(vector<double> d) {
    double base2 = sqrt(pow(d[0], 2) / (pow(d[0], 2) + pow(d[1], 2)));
    double base1 = -(base2 * d[1]) / d[0]; double base3 = 0;
    if (base2 - (d[1] / d[0]) * base1 <= 0) { return { base1, base2, base3 }; }
    else { return { -base1, -base2, base3 }; };
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

vector<double> RotateEnd(vector<double> end_2, vector<double> end_1, vector<double> origin_2, vector<double> d1, vector<double> d2, vector<double> b1, vector<double> b2, double gamma, double delta, double theta, double testa1, double testa2, double testa3)
{
    end_2[0] -= origin_2[0]; end_2[1] -= origin_2[1]; end_2[2] -= origin_2[2];

    double holder1; double holder2;

    if (b2[0] >= 0) {
        holder1 = end_2[0] * cos(gamma) - end_2[1] * sin(gamma);
        holder2 = end_2[0] * sin(gamma) + end_2[1] * cos(gamma);
        end_2[0] = holder1;
        end_2[1] = holder2;
    }
    else {
        holder1 = end_2[0] * cos(2 * M_PI - gamma) - end_2[1] * sin(2 * M_PI - gamma);
        holder2 = end_2[0] * sin(2 * M_PI - gamma) + end_2[1] * cos(2 * M_PI - gamma);
        end_2[0] = holder1;
        end_2[1] = holder2;
    }

    if (d2[2] <= 0) {
        holder1 = end_2[0] * cos(testa1) + end_2[2] * sin(testa1);
        holder2 = -end_2[0] * sin(testa1) + end_2[2] * cos(testa1);
        end_2[0] = holder1;
        end_2[2] = holder2;
    }
    else {
        holder1 = end_2[0] * cos(2 * M_PI - testa1) + end_2[2] * sin(2 * M_PI - testa1);
        holder2 = -end_2[0] * sin(2 * M_PI - testa1) + end_2[2] * cos(2 * M_PI - testa1);
        end_2[0] = holder1;
        end_2[2] = holder2;
    }

    if (true) {
        holder1 = end_2[1] * cos(testa2) - end_2[2] * sin(testa2);
        holder2 = end_2[1] * sin(testa2) + end_2[2] * cos(testa2);
        end_2[1] = holder1;
        end_2[2] = holder2;
    }
    /*else {
        holder1 = end_2[0] * cos(2 * M_PI - testa2) + end_2[2] * sin(2 * M_PI - testa2);
        holder2 = end_2[0] * sin(2 * M_PI - testa2) + end_2[2] * cos(2 * M_PI - testa2);
        end_2[1] = holder1;
        end_2[2] = holder2;
    }*/

    if (d1[2] >= 0) {
        holder1 = end_2[0] * cos(testa3) + end_2[2] * sin(testa3);
        holder2 = -end_2[0] * sin(testa3) + end_2[2] * cos(testa3);
        end_2[0] = holder1;
        end_2[2] = holder2;
    }
    else {
        holder1 = end_2[0] * cos(2 * M_PI - testa3) + end_2[2] * sin(2 * M_PI - testa3);
        holder2 = -end_2[0] * sin(2 * M_PI - testa3) + end_2[2] * cos(2 * M_PI - testa3);
        end_2[0] = holder1;
        end_2[2] = holder2;
    }

    if (b1[0] <= 0) {
        holder1 = end_2[0] * cos(theta) - end_2[1] * sin(theta);
        holder2 = end_2[0] * sin(theta) + end_2[1] * cos(theta);
        end_2[0] = holder1;
        end_2[1] = holder2;
    }
    else {
        holder1 = end_2[0] * cos(2 * M_PI - theta) - end_2[1] * sin(2 * M_PI - theta);
        holder2 = end_2[0] * sin(2 * M_PI - theta) + end_2[1] * cos(2 * M_PI - theta);
        end_2[0] = holder1;
        end_2[1] = holder2;
    }

    end_2[0] += end_1[0]; end_2[1] += end_1[1]; end_2[2] += end_1[2];

    return end_2;

}

vector<double> RotateEndVector(vector<double> end_2, vector<double> d1, vector<double> d2, vector<double> b1, vector<double> b2, double gamma, double delta, double theta, double testa1, double testa2, double testa3)
{
    double holder1; double holder2;

    if (b2[0] >= 0) {
        holder1 = end_2[0] * cos(gamma) - end_2[1] * sin(gamma);
        holder2 = end_2[0] * sin(gamma) + end_2[1] * cos(gamma);
        end_2[0] = holder1;
        end_2[1] = holder2;
    }
    else {
        holder1 = end_2[0] * cos(2 * M_PI - gamma) - end_2[1] * sin(2 * M_PI - gamma);
        holder2 = end_2[0] * sin(2 * M_PI - gamma) + end_2[1] * cos(2 * M_PI - gamma);
        end_2[0] = holder1;
        end_2[1] = holder2;
    }

    if (d2[2] <= 0) {
        holder1 = end_2[0] * cos(testa1) + end_2[2] * sin(testa1);
        holder2 = -end_2[0] * sin(testa1) + end_2[2] * cos(testa1);
        end_2[0] = holder1;
        end_2[2] = holder2;
    }
    else {
        holder1 = end_2[0] * cos(2 * M_PI - testa1) + end_2[2] * sin(2 * M_PI - testa1);
        holder2 = -end_2[0] * sin(2 * M_PI - testa1) + end_2[2] * cos(2 * M_PI - testa1);
        end_2[0] = holder1;
        end_2[2] = holder2;
    }

    if (true) {
        holder1 = end_2[1] * cos(testa2) - end_2[2] * sin(testa2);
        holder2 = end_2[1] * sin(testa2) + end_2[2] * cos(testa2);
        end_2[1] = holder1;
        end_2[2] = holder2;
    }
    /*else {
        holder1 = end_2[0] * cos(2 * M_PI - testa2) + end_2[2] * sin(2 * M_PI - testa2);
        holder2 = end_2[0] * sin(2 * M_PI - testa2) + end_2[2] * cos(2 * M_PI - testa2);
        end_2[1] = holder1;
        end_2[2] = holder2;
    }*/

    if (d1[2] >= 0) {
        holder1 = end_2[0] * cos(testa3) + end_2[2] * sin(testa3);
        holder2 = -end_2[0] * sin(testa3) + end_2[2] * cos(testa3);
        end_2[0] = holder1;
        end_2[2] = holder2;
    }
    else {
        holder1 = end_2[0] * cos(2 * M_PI - testa3) + end_2[2] * sin(2 * M_PI - testa3);
        holder2 = -end_2[0] * sin(2 * M_PI - testa3) + end_2[2] * cos(2 * M_PI - testa3);
        end_2[0] = holder1;
        end_2[2] = holder2;
    }

    if (b1[0] <= 0) {
        holder1 = end_2[0] * cos(theta) - end_2[1] * sin(theta);
        holder2 = end_2[0] * sin(theta) + end_2[1] * cos(theta);
        end_2[0] = holder1;
        end_2[1] = holder2;
    }
    else {
        holder1 = end_2[0] * cos(2 * M_PI - theta) - end_2[1] * sin(2 * M_PI - theta);
        holder2 = end_2[0] * sin(2 * M_PI - theta) + end_2[1] * cos(2 * M_PI - theta);
        end_2[0] = holder1;
        end_2[1] = holder2;
    }
    
    return end_2;
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

vector<double> VectorOnPlane2(vector<double> derivative, vector<double> b, double init_angle, double angle)
{
    double d_length = sqrt(pow(derivative[0], 2) + pow(derivative[1], 2) + pow(derivative[2], 2));
    vector<double> d = { derivative[0] / d_length, derivative[1] / d_length, derivative[2] / d_length };
    double rot_d_angle = acos(d[0] / (sqrt(pow(d[0], 2) + pow(d[1], 2)))); vector<double> rot_d = { 0,0,0 };

    if (d[1] < 0) { rot_d[0] = d[0] * cos(rot_d_angle) - d[1] * sin(rot_d_angle); rot_d[1] = d[0] * sin(rot_d_angle) + d[1] * cos(rot_d_angle); rot_d[2] = d[2]; }
    else { rot_d[0] = d[0] * cos(2 * M_PI - rot_d_angle) - d[1] * sin(2 * M_PI - rot_d_angle); rot_d[1] = d[0] * sin(2 * M_PI - rot_d_angle) + d[1] * cos(2 * M_PI - rot_d_angle); rot_d[2] = d[2]; }
    
    double angle_rot_d = acos(rot_d[0]);
    vector<double> a = { 0,cos(init_angle + angle),sin(init_angle + angle)}; double holder1 = 0; double holder2 = 0;
    
    if (rot_d[0] <= 0) {
        holder1 = a[0] * cos(angle_rot_d) + a[2] * sin(angle_rot_d);
        holder2 = -a[0] * sin(angle_rot_d) + a[2] * cos(angle_rot_d);
        a[0] = holder1; a[2] = holder2;
    }
    else {
        holder1 = a[0] * cos(2 * M_PI - angle_rot_d) + a[2] * sin(2 * M_PI - angle_rot_d);
        holder2 = -a[0] * sin(2 * M_PI - angle_rot_d) + a[2] * cos(2 * M_PI - angle_rot_d);
        a[0] = holder1; a[2] = holder2;
    }

    if (d[1] > 0) {
        holder1 = a[0] * cos(rot_d_angle) - a[1] * sin(rot_d_angle);
        holder2 = a[0] * sin(rot_d_angle) + a[1] * cos(rot_d_angle);
        a[0] = holder1; a[1] = holder2;
    }
    else {
        holder1 = a[0] * cos(2 * M_PI - rot_d_angle) - a[1] * sin(2 * M_PI - rot_d_angle);
        holder2 = a[0] * sin(2 * M_PI - rot_d_angle) + a[1] * cos(2 * M_PI - rot_d_angle);
        a[0] = holder1; a[1] = holder2;
    }
    
    return a;
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
    
    vector<double> d = { function->functions[0].PartialDerivative(values, function->variable_name)+limit, 
        function->functions[1].PartialDerivative(values, function->variable_name)+limit,
        function->functions[2].PartialDerivative(values, function->variable_name) +limit };
    double d_length = sqrt(pow(d[0], 2) + pow(d[1], 2) + pow(d[2], 2));

    vector<double> point_on_line = { function->functions[0].get(values), function->functions[1].get(values), function->functions[2].get(values)};

    vector<double> base = DecideBase(d);
    
    //Getting axes of elipse as having given angles from base and then a axis
    vector<double> a = VectorOnPlane2(d, { base[0], base[1], base[2]}, 0, a_angle);
    vector<double> b = VectorOnPlane2(d, { base[0], base[1], base[2]}, a_angle, M_PI / 2);
    
    a = { a[0] * a_length, a[1] * a_length, a[2] * a_length }; b = { b[0] * b_length, b[1] * b_length, b[2] * b_length };
    
    vector<double> run_vector = a;

    double run_angle = 0; vector<double> run_point{ 0,0,0 }; vector<vector<double>> points; vector<vector<int>> lines;

    while (run_angle < 2 * M_PI) {
        run_point[0] = a[0] * cos(run_angle) + b[0] * sin(run_angle) + point_on_line[0];
        run_point[1] = a[1] * cos(run_angle) + b[1] * sin(run_angle) + point_on_line[1];
        run_point[2] = a[2] * cos(run_angle) + b[2] * sin(run_angle) + point_on_line[2];

        points.emplace_back(run_point);

        run_angle += (2 * M_PI) / number_angles;
    }

    vector<int> found_points = {0};
    int current_point;
    double distance_point;
    int found_point;

    //Connecting individual points with lines to get connected elipse
    for (int i = 0; i < points.size()-1; i++) 
    { 
        current_point = found_points[0];
        found_point = current_point;
        distance_point = 1e10;
        for (int j = 0; j < points.size(); j++) 
        {
            if (Distance_here(points[current_point], points[j]) < distance_point) {
                if (std::find(found_points.begin(), found_points.end(), j) == found_points.end())
                {
                    distance_point = Distance_here(points[current_point], points[j]);
                    found_point = j;
                }
            }
        }
        found_points.insert(found_points.begin(), found_point);
        lines.push_back({ current_point, found_point });
    }
    lines.push_back({ found_points[0], 0 });

    Elipse elipse = Elipse(points, lines, a, b, point_on_line, a_angle, number_angles);
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
        if (i!=0){ current_point = FindEndPoint(gene->function, gene->distance, current_point); }
        futures.emplace_back(std::async(std::launch::deferred, ElipsePoints, gene->function, current_point, gene->elipses_values[i][0],
            gene->elipses_values[i][1], gene->elipses_values[i][2], gene->function->get(current_point), gene->resolution));
    }
    
    for (const std::future<Elipse>& ftr : futures)
    {
        ftr.wait();
        Elipse some = ftr._Get_value();
        container.elipses.emplace_back(make_unique<Elipse>(some));
    }

    container.end = gene->function->get(current_point);
    container.end_point = current_point;

    container.function = gene->function;
    
    return container;
}

double Length(double a1, double a2, double a3) 
{
    return sqrt(pow(a1, 2) + pow(a2, 2) + pow(a3, 2));
}

vector<unique_ptr<Elipse>> GetElipseConnect(ProkaryotePartContainer pp1, ProkaryotePartContainer pp2) 
{
    vector<unique_ptr<Elipse>> ConnectedElipses;

    for (int i = 0; i < pp1.elipses.size(); i++)
    {
        ConnectedElipses.emplace_back(move(pp1.elipses[i]));
    }

    for (int i = 0; i < pp2.elipses.size(); i++)
    {
        ConnectedElipses.emplace_back(move(pp2.elipses[i]));
    }

    return ConnectedElipses;
}

double GetEndAngle(vector<double> end, vector<double> base, vector<double> d, double angle,  vector<vector<double>> elipse_points) {
    vector<double> furthest_point; double distance = 0; double run_distance;
    double holder1; double holder2;

    for (int i = 0; i < elipse_points.size(); i++) {
        run_distance = Distance_here(end, elipse_points[i]);
        if (run_distance > distance) { distance = run_distance; furthest_point = elipse_points[i]; }
    }

    vector<double> vec = { furthest_point[0] - end[0],furthest_point[1] - end[1], furthest_point[2] - end[2] };
    double vec_length = sqrt(pow(vec[0], 2) + pow(vec[1], 2) + pow(vec[2], 2));
    vec[0] = vec[0] / vec_length; vec[1] = vec[1] / vec_length; vec[2] = vec[2] / vec_length;
    return acos(vec[0] * base[0] + vec[1] * base[1] + vec[2] * base[2]) - angle;
}

ProkaryotePartContainer ConnectParts(ProkaryotePartContainer prokaryote_part_1, ProkaryotePartContainer prokaryote_part_2, int counter, int run)
{
    double limit = 1e-10;

    //function to b modified during the runtime and serve as connected part new one
    VectorFunction function_run = *prokaryote_part_2.function;
    vector<double> d1 = { 0,0,0 };
    double holder1; double holder2; bool holder3 = true;
    
    if (run == 0) {
        map<string, double> values1 = { {prokaryote_part_1.function->variable_name, prokaryote_part_1.end_point} };
        d1[0] = prokaryote_part_1.function->functions[0].PartialDerivative(values1, prokaryote_part_1.function->variable_name) + limit;
        d1[1] = prokaryote_part_1.function->functions[1].PartialDerivative(values1, prokaryote_part_1.function->variable_name) + limit;
        d1[2] = prokaryote_part_1.function->functions[2].PartialDerivative(values1, prokaryote_part_1.function->variable_name) + limit;
    }
    else {
        d1[0] = prokaryote_part_1.derivative_end[0];
        d1[1] = prokaryote_part_1.derivative_end[1];
        d1[2] = prokaryote_part_1.derivative_end[2];
    }
    
    map<string, double> values2 = { {prokaryote_part_2.function->variable_name, prokaryote_part_2.origin_point} };
    vector<double> d2 = { prokaryote_part_2.function->functions[0].PartialDerivative(values2, prokaryote_part_2.function->variable_name) + limit,
        prokaryote_part_2.function->functions[1].PartialDerivative(values2, prokaryote_part_2.function->variable_name) + limit,
        prokaryote_part_2.function->functions[2].PartialDerivative(values2, prokaryote_part_2.function->variable_name) + limit };
    
    double d1_length = sqrt(pow(d1[0], 2) + pow(d1[1], 2) + pow(d1[2], 2));
    double d2_length = sqrt(pow(d2[0], 2) + pow(d2[1], 2) + pow(d2[2], 2));
    
    d1[0] = d1[0] / d1_length; d1[1] = d1[1] / d1_length; d1[2] = d1[2] / d1_length;
    d2[0] = d2[0] / d2_length; d2[1] = d2[1] / d2_length; d2[2] = d2[2] / d2_length;
    
    vector<double> b1 = DecideBase(d1);
    vector<double> b2 = DecideBase(d2);

    vector<double> rotated_d1 = { 0,0,0 }; vector<double> rotated_d2 = { 0,0,0 }; 
    double alpha = acos(d1[0] / sqrt(pow(d1[0], 2) + pow(d1[1], 2))); 
    double beta = acos(d2[0] / sqrt(pow(d2[0], 2) + pow(d2[1], 2)));
    
    if (d1[1] <= 0) {
        rotated_d1[0] = d1[0] * cos(alpha) - d1[1] * sin(alpha);
        rotated_d1[1] = d1[0] * sin(alpha) + d1[1] * cos(alpha);
        rotated_d1[2] = d1[2];
    }
    else {
        rotated_d1[0] = d1[0] * cos(2 * M_PI - alpha) - d1[1] * sin(2 * M_PI - alpha);
        rotated_d1[1] = d1[0] * sin(2 * M_PI - alpha) + d1[1] * cos(2 * M_PI - alpha);
        rotated_d1[2] = d1[2];
    }
    
    if (d2[1] <= 0) {
        rotated_d2[0] = d2[0] * cos(beta) - d2[1] * sin(beta);
        rotated_d2[1] = d2[0] * sin(beta) + d2[1] * cos(beta);
        rotated_d2[2] = d2[2];
    }
    else {
        rotated_d2[0] = d2[0] * cos(2 * M_PI - beta) - d2[1] * sin(2 * M_PI - beta);
        rotated_d2[1] = d2[0] * sin(2 * M_PI - beta) + d2[1] * cos(2 * M_PI - beta);
        rotated_d2[2] = d2[2];
    }
    
    double gamma = acos(b2[1]);
    double delta = acos((rotated_d1[0] * rotated_d2[0] + rotated_d1[2] * rotated_d2[2]));
    double theta = acos(b1[1]);
    
    //Get second part upon origin
    for (int i = 0; i < prokaryote_part_2.elipses.size(); i++) {
        for (int j = 0; j < prokaryote_part_2.elipses[i]->points.size(); j++)
        {
            prokaryote_part_2.elipses[i]->points[j][0] -= prokaryote_part_2.origin[0];
            prokaryote_part_2.elipses[i]->points[j][1] -= prokaryote_part_2.origin[1];
            prokaryote_part_2.elipses[i]->points[j][2] -= prokaryote_part_2.origin[2];
        }
    }

    //Rotate second part's base to y axis
    for (int i = 0; i < prokaryote_part_2.elipses.size(); i++) {
        for (int j = 0; j < prokaryote_part_2.elipses[i]->points.size(); j++)
        {
            if (b2[0] >= 0) {
                holder1 = prokaryote_part_2.elipses[i]->points[j][0] * cos(gamma) - prokaryote_part_2.elipses[i]->points[j][1] * sin(gamma);
                holder2 = prokaryote_part_2.elipses[i]->points[j][0] * sin(gamma) + prokaryote_part_2.elipses[i]->points[j][1] * cos(gamma);
                prokaryote_part_2.elipses[i]->points[j][0] = holder1;
                prokaryote_part_2.elipses[i]->points[j][1] = holder2;
            }
            else {
                holder1 = prokaryote_part_2.elipses[i]->points[j][0] * cos(2 * M_PI - gamma) - prokaryote_part_2.elipses[i]->points[j][1] * sin(2 * M_PI - gamma);
                holder2 = prokaryote_part_2.elipses[i]->points[j][0] * sin(2 * M_PI - gamma) + prokaryote_part_2.elipses[i]->points[j][1] * cos(2 * M_PI - gamma);
                prokaryote_part_2.elipses[i]->points[j][0] = holder1;
                prokaryote_part_2.elipses[i]->points[j][1] = holder2;
            }
        }
    }

    double testa1 = acos(rotated_d2[0]); double testa2 = prokaryote_part_1.end_angle; double testa3 = acos(rotated_d1[0]);
    
    for (int i = 0; i < prokaryote_part_2.elipses.size(); i++) {
        for (int j = 0; j < prokaryote_part_2.elipses[i]->points.size(); j++)
        {
            if (rotated_d2[2] <= 0) {
                holder1 = prokaryote_part_2.elipses[i]->points[j][0] * cos(testa1) + prokaryote_part_2.elipses[i]->points[j][2] * sin(testa1);
                holder2 = -prokaryote_part_2.elipses[i]->points[j][0] * sin(testa1) + prokaryote_part_2.elipses[i]->points[j][2] * cos(testa1);
                prokaryote_part_2.elipses[i]->points[j][0] = holder1;
                prokaryote_part_2.elipses[i]->points[j][2] = holder2;
            }
            else {
                holder1 = prokaryote_part_2.elipses[i]->points[j][0] * cos(2 * M_PI - testa1) + prokaryote_part_2.elipses[i]->points[j][2] * sin(2 * M_PI - testa1);
                holder2 = -prokaryote_part_2.elipses[i]->points[j][0] * sin(2 * M_PI - testa1) + prokaryote_part_2.elipses[i]->points[j][2] * cos(2 * M_PI - testa1);
                prokaryote_part_2.elipses[i]->points[j][0] = holder1;
                prokaryote_part_2.elipses[i]->points[j][2] = holder2;
            }
        }
    }
    
    for (int i = 0; i < prokaryote_part_2.elipses.size(); i++) {
        for (int j = 0; j < prokaryote_part_2.elipses[i]->points.size(); j++)
        {
            if (true) {
                holder1 = prokaryote_part_2.elipses[i]->points[j][1] * cos(testa2) - prokaryote_part_2.elipses[i]->points[j][2] * sin(testa2);
                holder2 = prokaryote_part_2.elipses[i]->points[j][1] * sin(testa2) + prokaryote_part_2.elipses[i]->points[j][2] * cos(testa2);
                prokaryote_part_2.elipses[i]->points[j][1] = holder1;
                prokaryote_part_2.elipses[i]->points[j][2] = holder2;
            }
            else {
                holder1 = prokaryote_part_2.elipses[i]->points[j][0] * cos(2 * M_PI - testa2) + prokaryote_part_2.elipses[i]->points[j][2] * sin(2 * M_PI - testa2);
                holder2 = prokaryote_part_2.elipses[i]->points[j][0] * sin(2 * M_PI - testa2) + prokaryote_part_2.elipses[i]->points[j][2] * cos(2 * M_PI - testa2);
                prokaryote_part_2.elipses[i]->points[j][1] = holder1;
                prokaryote_part_2.elipses[i]->points[j][2] = holder2;
            }
        }
    }

    for (int i = 0; i < prokaryote_part_2.elipses.size(); i++) {
        for (int j = 0; j < prokaryote_part_2.elipses[i]->points.size(); j++)
        {
            if (rotated_d1[2] >= 0) {
                holder1 = prokaryote_part_2.elipses[i]->points[j][0] * cos(testa3) + prokaryote_part_2.elipses[i]->points[j][2] * sin(testa3);
                holder2 = -prokaryote_part_2.elipses[i]->points[j][0] * sin(testa3) + prokaryote_part_2.elipses[i]->points[j][2] * cos(testa3);
                prokaryote_part_2.elipses[i]->points[j][0] = holder1;
                prokaryote_part_2.elipses[i]->points[j][2] = holder2;
            }
            else {
                holder1 = prokaryote_part_2.elipses[i]->points[j][0] * cos(2 * M_PI - testa3) + prokaryote_part_2.elipses[i]->points[j][2] * sin(2 * M_PI - testa3);
                holder2 = -prokaryote_part_2.elipses[i]->points[j][0] * sin(2 * M_PI - testa3) + prokaryote_part_2.elipses[i]->points[j][2] * cos(2 * M_PI - testa3);
                prokaryote_part_2.elipses[i]->points[j][0] = holder1;
                prokaryote_part_2.elipses[i]->points[j][2] = holder2;
            }
        }
    }

    //Rotate second part's base to first part's base
    for (int i = 0; i < prokaryote_part_2.elipses.size(); i++) {
         for (int j = 0; j < prokaryote_part_2.elipses[i]->points.size(); j++)
         {
             if (b1[0] <= 0) {
                 holder1 = prokaryote_part_2.elipses[i]->points[j][0] * cos(theta) - prokaryote_part_2.elipses[i]->points[j][1] * sin(theta);
                 holder2 = prokaryote_part_2.elipses[i]->points[j][0] * sin(theta) + prokaryote_part_2.elipses[i]->points[j][1] * cos(theta);
                 prokaryote_part_2.elipses[i]->points[j][0] = holder1;
                 prokaryote_part_2.elipses[i]->points[j][1] = holder2;
             }
             else {
                 holder1 = prokaryote_part_2.elipses[i]->points[j][0] * cos(2 * M_PI - theta) - prokaryote_part_2.elipses[i]->points[j][1] * sin(2 * M_PI - theta);
                 holder2 = prokaryote_part_2.elipses[i]->points[j][0] * sin(2 * M_PI - theta) + prokaryote_part_2.elipses[i]->points[j][1] * cos(2 * M_PI - theta);
                 prokaryote_part_2.elipses[i]->points[j][0] = holder1;
                 prokaryote_part_2.elipses[i]->points[j][1] = holder2;
             }
         }
    }

    //Get second part upon origin
    for (int i = 0; i < prokaryote_part_2.elipses.size(); i++) {
         for (int j = 0; j < prokaryote_part_2.elipses[i]->points.size(); j++)
         {
             prokaryote_part_2.elipses[i]->points[j][0] += prokaryote_part_1.end[0];
             prokaryote_part_2.elipses[i]->points[j][1] += prokaryote_part_1.end[1];
             prokaryote_part_2.elipses[i]->points[j][2] += prokaryote_part_1.end[2];
         }
    }

    vector<double> d3 = { 0,0,0 };
    map<string, double> values3 = { {prokaryote_part_2.function->variable_name, prokaryote_part_2.end_point } };
    d3[0] = prokaryote_part_2.function->functions[0].PartialDerivative(values3, prokaryote_part_2.function->variable_name) + limit;
    d3[1] = prokaryote_part_2.function->functions[1].PartialDerivative(values3, prokaryote_part_2.function->variable_name) + limit;
    d3[2] = prokaryote_part_2.function->functions[2].PartialDerivative(values3, prokaryote_part_2.function->variable_name) + limit;

    double d3_length = sqrt(pow(d3[0], 2) + pow(d3[1], 2) + pow(d3[2], 2));

    d3[0] = d3[0]/d3_length; d3[1] = d3[1] / d3_length; d3[2] = d3[2] / d3_length;
    
    d3 = RotateEndVector(d3, rotated_d1, rotated_d2, b1, b2, gamma, delta, theta, testa1, testa2, testa3);
    vector<double> b3 = DecideBase(d3);

    vector<double> end = RotateEnd(prokaryote_part_2.end, prokaryote_part_1.end, prokaryote_part_2.origin, rotated_d1, rotated_d2, b1, b2, gamma, delta, theta, testa1, testa2, testa3);
    
    shared_ptr<VectorFunction> function = make_shared<VectorFunction>(function_run);
    double origin_point = prokaryote_part_2.origin_point;
    vector<double> origin = prokaryote_part_1.end;
    double end_point = prokaryote_part_2.end_point;

    double end_angle = GetEndAngle(end, b3, d3, prokaryote_part_2.elipses[prokaryote_part_2.elipses.size() - 1]->a_angle,  prokaryote_part_2.elipses[prokaryote_part_2.elipses.size() - 1]->points);
    
    vector<unique_ptr<Elipse>> ConnectedElipses = GetElipseConnect(move(prokaryote_part_1), move(prokaryote_part_2));

    ProkaryotePartContainer ConnectedPart = ProkaryotePartContainer(function, move(ConnectedElipses), origin_point, origin, end_point, end);
    
    ConnectedPart.derivative_end = d3;
    ConnectedPart.end_angle = end_angle;

    return ConnectedPart;
}

ProkaryoteBodyContainer AssembleParts(vector<ProkaryotePartContainer> prokaryote_parts)
{
    if (prokaryote_parts.size() == 1) { ProkaryoteBodyContainer body = ProkaryoteBodyContainer(move(prokaryote_parts[0])); return body; }
    else if (prokaryote_parts.size() == 2) { 
        ProkaryotePartContainer running_container = ConnectParts(move(prokaryote_parts[0]), move(prokaryote_parts[1]), 0, 0);
        ProkaryoteBodyContainer body = ProkaryoteBodyContainer(move(running_container)); 
        return body; 
    }
    else {
        
        int counter = prokaryote_parts[0].elipses.size() + prokaryote_parts[1].elipses.size();
        ProkaryotePartContainer running_container = ConnectParts(move(prokaryote_parts[0]), move(prokaryote_parts[1]), 0, 0);
        int run = 1;
        for (int i = 2; i < prokaryote_parts.size(); i++)
        {
            running_container = ConnectParts(move(running_container), move(prokaryote_parts[i]), counter, run);
            run += 1;
        }

        ProkaryoteBodyContainer body = ProkaryoteBodyContainer(move(running_container));
        return body;
    }
}


