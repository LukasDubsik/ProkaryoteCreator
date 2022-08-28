
#include <chrono>

#include "function.h"
#include "prokaryote_body_creator.h"
#include "export.h"

using namespace std;

//Create the list of elipses per the part of bacteria

/*Create the end parts of bacterium*/
//Some variation, need to think about

/*Connecting functions together*/
//Here it is needed to think about how to connect different functional parts

int main()
{
    Function func1 = Function("var(;t;)", vector<string> {"t"}); //declaring the class
    Function func2 = Function("div(;var(;t;);con(;4;);)", vector<string> {"t"}); //declaring the class
    Function func3 = Function("dec(;mul(;var(;t;);con(;3;););con(;1;);)", vector<string> {"t"}); //declaring the class
    
    vector<Function> functions1 {func1, func2, func3};

    shared_ptr<VectorFunction> funcv1(new VectorFunction(functions1, "t"));
    
    unique_ptr<EndBodyPartGene> gene_begin(new EndBodyPartGene(funcv1, 0.100, 5, { 1.000, 0.100, 1.000 }, 30, 0.9, "begin"));
    unique_ptr<EndBodyPartGene> gene_end(new EndBodyPartGene(funcv1, 0.100, 4, { 1.000, 0.100, 1.000 }, 30, 0.9, "end"));

    unique_ptr<MainBodyPartGene> gene_body(new MainBodyPartGene(funcv1, 0.500, 7, { {1.000,0.100,1.000}, {1.000,0.100,1.000}, {1.000,0.100,1.000}, {1.000,0.100,1.000}, {1.000,0.100,1.000}, {1.000,0.100,1.000}, {1.000,0.100,1.000} }, 30));
    unique_ptr<MainBodyPartGene> gene_begin_body = ConvertEndToMain(move(gene_begin));
    unique_ptr<MainBodyPartGene> gene_end_body = ConvertEndToMain(move(gene_end));
    
    ProkaryotePartContainer res1 = AssembleContainerBody(std::move(gene_begin_body));
    ProkaryotePartContainer res2 = AssembleContainerBody(std::move(gene_body));
    ProkaryotePartContainer res3 = AssembleContainerBody(std::move(gene_end_body));

    vector<ProkaryotePartContainer> parts;
    parts.emplace_back(move(res1));
    parts.emplace_back(move(res2));
    parts.emplace_back(move(res3));
    ProkaryoteBodyContainer assembly = AssembleParts(move(parts));
    
    ExportToBlender(move(assembly), "C:", "test1.txt");

    return 0;
}
