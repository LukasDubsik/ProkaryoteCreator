
#include <chrono>

#include "function.h"
#include "prokaryote_body_creator.h"
#include "export.h"

using namespace std;

int main()
{
    //These functions define vector function which serves as scafold upon which perpendicular elipses, with additional instructions hidden in genes (below), are created
    //As they are part of vector function, all functions must have only one variable
    //Further instructions, on how to write functions, can be found in this repository: https://github.com/LukasDubsik/FunctionCreator.git
    Function func11 = Function("cos(;var(;t;););", vector<string> {"t"}); //declaring the class
    Function func21 = Function("sin(;var(;t;););", vector<string> {"t"}); //declaring the class
    Function func31 = Function("var(;t;);", vector<string> {"t"}); //declaring the class
    
    vector<Function> functions1 {func11, func21, func31};
    
    //The vector function must be shared vector, this serves for more effective passing
    shared_ptr<VectorFunction> funcv1(new VectorFunction(functions1, "t"));

    Function func12 = Function("var(;t;);", vector<string> {"t"}); //declaring the class
    Function func22 = Function("con(;0;);", vector<string> {"t"}); //declaring the class
    Function func32 = Function("con(;0;);", vector<string> {"t"}); //declaring the class

    vector<Function> functions2{ func12, func22, func32 };
    
    //The vector function must be shared vector, this serves for more effective passing
    shared_ptr<VectorFunction> funcv2(new VectorFunction(functions2, "t"));
    
    //Genes to define begging cap and ending cap of prokaryote
    //The definition, on what should be included, can be found in prokaryote_body_gene header
    unique_ptr<EndBodyPartGene> gene_begin(new EndBodyPartGene(funcv2, 0.100, 10, { 1.000, 0.000, 1.000 }, 20, "begin"));
    unique_ptr<EndBodyPartGene> gene_end(new EndBodyPartGene(funcv2, 0.100, 10, { 1.000, 0.000, 1.000 }, 20, "end"));
    
    //main body gene is defined (there can be more of them, but here only one) with capping (begin/end) genes converted to main, as this format 
    //is the only that can be inserted to connecting functions
    unique_ptr<MainBodyPartGene> gene_body(new MainBodyPartGene(funcv2, 0.200, 40, { 1.000, 0.000, 1.000 }, 20));
    unique_ptr<MainBodyPartGene> gene_begin_body = ConvertEndToMain(move(gene_begin));
    unique_ptr<MainBodyPartGene> gene_end_body = ConvertEndToMain(move(gene_end));

    //Genes are converted to body containers, these contain elipse points and other necessary information for export
    ProkaryotePartContainer res1 = AssembleContainerBody(std::move(gene_begin_body));
    ProkaryotePartContainer res2 = AssembleContainerBody(std::move(gene_body));
    ProkaryotePartContainer res3 = AssembleContainerBody(std::move(gene_end_body));

    //Vector of parts is created, with order determing the order of assembling in export
    vector<ProkaryotePartContainer> parts;
    parts.emplace_back(move(res1)); 
    parts.emplace_back(move(res2));
    parts.emplace_back(move(res3));

    //Parts are assembled creating final body conatiner with all parts present in previosly defined order
    ProkaryoteBodyContainer assembly = AssembleParts(move(parts));
    
    //assembly is exported, creating test1.txt file to be used in blender to visualize the result.
    ExportToBlender(move(assembly), "C:", "test1.txt", 0.025, 20, "CLOUDS");

    return 0;
}
