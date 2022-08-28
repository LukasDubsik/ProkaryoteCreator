#ifndef EXPORT_H
#define EXPORT_H

#include <iostream>
#include <fstream>

#include "prokaryote_body_creator.h"

using namespace std;

//Defines functions used to export parts of the genome for assemblies and similar functions
//Export prokaryote structure for blender assembly
void ExportToBlender(ProkaryoteBodyContainer assembly, string export_path, string file_name);

#endif
