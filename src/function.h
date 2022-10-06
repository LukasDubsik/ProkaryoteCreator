#ifndef FUNCTION_H
#define FUNCTION_H

#include <iostream>
#include <string>
#include <stdexcept>
#include <vector>
#include <map>
#include <math.h>
#include <mutex>

using namespace std;

class Function
{
public:
    std::map<string, double> variable_values;

    Function(std::string function_def, std::vector<std::string> variables);

    string getDefinition();
    void setDefinition(string function_definition);

    double get(map<string, double> variable_val);
    double Integrate(map<string, double> variable_val, string variable_name, double a, double b, int n);
    double PartialDerivative(map<string, double> variable_val, string variable_name, double limit = 1e-10);

private:
    std::vector<std::string> _variables;
    std::vector<std::string> _definition;
    std::string _function_definition;
    int _definition_length;

    void Preprocess();

    void CheckVariableValues();
    template <typename type> bool FindInVector(vector<type> vec, type value);
    double EvaluateVariable(string type, string value);

    vector<double> OperatorsBasic(int position, string oper);
    vector<double> OperatorsFunction(int position, string function);

    double Basic(double val1, double val2, string oper);
    double Func(double value, string oper);

protected:
    vector<string> _operators_function{ "sin(", "cos(", "tan(", "asin(", "acos(", "atan(",
        "cosh(", "sinh(", "tanh(", "log(", "abs(" };
    vector<string> _operators_basic{ "sum(", "dec(", "mul(", "div(", "pow(" };
    vector<string> _operators_variable{ "con(", "var(" };
};

class VectorFunction 
{
public:
    vector<Function> functions;
    string variable_name;
    int function_size;

    VectorFunction(vector<Function> funcs, string var_name);

    vector<double> get(double variable_value);
    double SizeFunctionValue(double a, double b, int n = 100);
    double ValueToPlane(int type, double v_1, double v_2, double variable_value);

private:
    double CombineValues(double value);
};

VectorFunction RotateFunction(VectorFunction input_function, double angle, string axis);
VectorFunction SumFunction(VectorFunction input_function, double value_x, double value_y, double value_z);
VectorFunction DecrementFunction(VectorFunction input_function, double value_x, double value_y, double value_z);

#endif

