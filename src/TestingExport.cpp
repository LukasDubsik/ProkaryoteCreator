
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
    //These functions define vector function which serves as scafold upon which perpendicular elipses, with additional instructions hidden in genes (below), are created
    //As they are part of vector function, all functions must have only one variable
    //Further instructions, on how to write functions, can be found in this repository: https://github.com/LukasDubsik/FunctionCreator.git
    Function func1 = Function("var(;t;)", vector<string> {"t"}); //declaring the class
    Function func2 = Function("con(;0;)", vector<string> {"t"}); //declaring the class
    Function func3 = Function("con(;0;)", vector<string> {"t"}); //declaring the class
    
    vector<Function> functions1 {func1, func2, func3};
    
    //The vector function must be shared vector, this serves for more effective passing
    shared_ptr<VectorFunction> funcv1(new VectorFunction(functions1, "t"));
    
    //Genes to define begging cap and ending cap of prokaryote
    //The definition, on what should be included, can be found in prokaryote_body_gene header
    unique_ptr<EndBodyPartGene> gene_begin(new EndBodyPartGene(funcv1, 0.100, 10, { 1.000, 0.100, 1.000 }, 30, "begin"));
    unique_ptr<EndBodyPartGene> gene_end(new EndBodyPartGene(funcv1, 0.100, 10, { 1.000, 0.100, 1.000 }, 30, "end"));
    
    //main body gene is defined (there can be more of them, but here only one) with capping (begin/end) genes converted to main, as this format 
    //is the only that can be inserted to connecting functions
    unique_ptr<MainBodyPartGene> gene_body(new MainBodyPartGene(funcv1, 0.500, 7, { {1.000,0.100,1.000}, {1.000,0.100,1.000}, {1.000,0.100,1.000}, {1.000,0.100,1.000}, {1.000,0.100,1.000}, {1.000,0.100,1.000}, {1.000,0.100,1.000} }, 30));
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
    
    //assembly is exported, creating tesrt1.txt file to be used in blender to visualize the result.
    ExportToBlender(move(assembly), "C:", "test1.txt");

    return 0;
}
