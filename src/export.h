#ifndef EXPORT_H
#define EXPORT_H

#include <iostream>
#include <fstream>
#include <windows.h>
#include <string>
#include <stdio.h>

#include "prokaryote_body_creator.h"

using namespace std;

//Defines functions used to export parts of the genome for assemblies and similar functions
//Export prokaryote structure for blender assembly
void ExportToBlender(ProkaryoteBodyContainer assembly, string export_path, string file_name, int dimensions = 512, int r_color = 180, int g_color = 180, int b_color = 180);

#endif
