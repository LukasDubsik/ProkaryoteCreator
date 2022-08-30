
#include "export.h"

//Values to export in object blender creation
vector<string> blender_includes = { "import bpy", "import math"};

vector<string> setting_scene = { "#Set the viewing scene for the displacement export",
"bpy.data.scenes['Scene'].render.engine = 'CYCLES'",
"bpy.data.scenes['Scene'].cycles.feature_set = 'EXPERIMENTAL'"};

vector<string> point_cloud_function = { "#function for construction of the bacterial body",
	"def point_cloud(ob_name, coords, edges=[], faces=[]):",
"    me = bpy.data.meshes.new(ob_name + 'Mesh')",
"    ob = bpy.data.objects.new(ob_name, me)",
"    me.from_pydata(coords, edges, faces)",
"    ob.show_name = True", "    me.update()", "    return ob" };

vector<string> select_object = { "#make the mesh active",
"obj = bpy.data.objects['prokaryote_body']",
"bpy.context.view_layer.objects.active = bpy.data.objects['prokaryote_body']",
" ", "#Set the scene",
"bpy.ops.object.mode_set(mode = 'EDIT')",
"bpy.ops.mesh.select_mode(type = 'VERT')",
"bpy.ops.mesh.select_all(action = 'DESELECT')",
"bpy.ops.object.mode_set(mode = 'OBJECT')"};

vector<string> vertices_merge = { "bpy.ops.object.mode_set(mode = 'EDIT')",
"bpy.ops.mesh.merge(type = 'CENTER')",
"bpy.ops.mesh.select_all(action = 'DESELECT')",
"bpy.ops.object.mode_set(mode = 'OBJECT')"};

vector<string> object_smooth = { "#Smooth the surface",
"obj.select_set(True)",
"obj.data.use_auto_smooth = 1",
"obj.data.auto_smooth_angle = math.pi / 6",
"bpy.ops.object.shade_smooth()"};

string point_cloud_comment = "#Create the collection of points, lines and planes";
string point_cloud = "pc = point_cloud(";

string point_cloud_name = "'prokaryote_body'";

string link_object_comment = "#Link the the collection to the scene";
string link_object = "bpy.context.collection.objects.link(pc)";

double Distance(vector<double> a, vector<double> b)
{
	return sqrt(pow((a[0] - b[0]), 2) + pow((a[1] - b[1]), 2) + pow((a[2] - b[2]), 2));
}

int ChoosePosition(int length, int closest) 
{
	int close;
	if (closest == length) { return 0; }
	else if (closest == -1) { return length - 1; }
	else { return closest; }
}

void ExportToBlender(ProkaryoteBodyContainer assembly, string export_path, string file_name) 
{
	ofstream blender_file(file_name);

	for (string line : blender_includes) { blender_file << line + "\n"; }
	blender_file << "\n";
	for (string line : point_cloud_function) { blender_file << line + "\n"; }
	blender_file << "\n";

	blender_file << point_cloud_comment << "\n";
	string object_creation = point_cloud + point_cloud_name + ", ";
	object_creation += "[";

	for (int i = 0; i < assembly.elipses.size(); i++)
	{
		for (int j = 0; j < assembly.elipses[i]->points.size(); j++)
		{
			object_creation += "(" + to_string(assembly.elipses[i]->points[j][0]) + ", " +
				to_string(assembly.elipses[i]->points[j][1]) + ", " + to_string(assembly.elipses[i]->points[j][2]) + "), ";
		}
	}
	object_creation += "], [";

	int counter = 0;
	for (int i = 0; i < assembly.elipses.size(); i++)
	{
		for (int j = 0; j < assembly.elipses[i]->lines.size(); j++)
		{
			object_creation += "(" + to_string(assembly.elipses[i]->lines[j][0]+counter) + ", " + to_string(assembly.elipses[i]->lines[j][1]+counter) + "), ";
		}
		counter += assembly.elipses[i]->points.size();
	}
	object_creation += "], [";
	
	counter = 0;
	double closest_distance;
	int closest;
	int closest_next;
	int closestt;
	double distance;
	
	for (int i = 0; i < assembly.elipses.size() - 1; i++)
	{
		closest_distance = numeric_limits<float>::max();
		closest = 0;

		if (i == 0){
			object_creation += "(";
			for (int s = 0; s < assembly.elipses[i]->points.size(); s++) 
			{
				object_creation += to_string(s) + ", ";
			}
			object_creation += "), ";
		}

		for (int j = 0; j < assembly.elipses[i + 1]->lines.size(); j++) 
		{
			distance = Distance(assembly.elipses[i]->points[0], assembly.elipses[i + 1]->points[j]);
			if (distance < closest_distance) { closest_distance = distance; closest = j; }

		}
		
		closestt = closest;
		if (Distance(assembly.elipses[i]->points[1], assembly.elipses[i + 1]->points[ChoosePosition(assembly.elipses[i + 1]->points.size(), closest + 1)]) <=
			Distance(assembly.elipses[i]->points[1], assembly.elipses[i + 1]->points[ChoosePosition(assembly.elipses[i + 1]->points.size(), closest - 1)])) {
			if (closest == assembly.elipses[i + 1]->points.size() - 1) {closest_next = 0;}
			else { closest_next = closest + 1; }
			for (int k = 0; k < assembly.elipses[i]->points.size(); k++)
			{
				if (closest == assembly.elipses[i + 1]->points.size()) (closest = 0);
				if (closest == assembly.elipses[i + 1]->points.size() - 1) (closest_next = 0);
				object_creation += "(" + to_string(assembly.elipses[i]->lines[k][0] + counter) + ", " + to_string(assembly.elipses[i]->lines[k][1] + counter) + ", " + to_string(closest_next + counter + assembly.elipses[i]->points.size()) + ", " + to_string(closest + counter + assembly.elipses[i]->points.size()) + "), ";
				closest++;
				closest_next++;
			}
		}
		else {
			if (closest == 0) (closest_next = assembly.elipses[i + 1]->points.size() - 1);
			else (closest_next = closest - 1);
			for (int k = 0; k < assembly.elipses[i]->points.size(); k++)
			{
				if (closest == -1) (closest = assembly.elipses[i + 1]->points.size() - 1);
				if (closest == 0) (closest_next = assembly.elipses[i + 1]->points.size() - 1);
				object_creation += "(" + to_string(assembly.elipses[i]->lines[k][0] + counter) + ", " + to_string(assembly.elipses[i]->lines[k][1] + counter) + ", " + to_string(closest_next + counter + assembly.elipses[i]->points.size()) + ", " + to_string(closest + counter + assembly.elipses[i]->points.size()) + "), ";
				closest--;
				closest_next--;
			}
		}
		
		counter += assembly.elipses[i]->points.size();

	}

	object_creation += "(";
	for (int s = 0; s < assembly.elipses[assembly.elipses.size()-1]->points.size(); s++)
	{
		object_creation += to_string(s+counter) + ", ";
	}
	object_creation += ") ";

	object_creation += "])";

	blender_file << object_creation + "\n" << "\n";
	blender_file << link_object_comment << "\n";
	blender_file << link_object + "\n" << "\n";

	for (string line : select_object) { blender_file << line + "\n"; }
	blender_file << "\n";

	blender_file << "#Select all points on one side of cap and merge them" << "\n";
	string for_loop = "for i in range(0, " + to_string(assembly.elipses[0]->points.size())  +"):";
	blender_file << for_loop + "\n";
	blender_file << "    obj.data.vertices[i].select = True" + std::string("\n");

	for (string line : vertices_merge) { blender_file << line + "\n"; }
	blender_file << "\n";

	blender_file << "#Select all points on the other side and merge them" << "\n";
	for_loop = "for i in range(" + to_string(assembly.elipses.size() * assembly.elipses[assembly.elipses.size() - 1]->points.size() - assembly.elipses[0]->points.size() - 59) + "," + to_string(assembly.elipses.size() * assembly.elipses[assembly.elipses.size() - 1]->points.size() - 59) + "):";
	blender_file << for_loop + "\n";
	blender_file << "    obj.data.vertices[i].select = True" + std::string("\n");

	for (string line : vertices_merge) { blender_file << line + "\n"; }
	blender_file << "\n";

	for (string line : object_smooth) { blender_file << line + "\n"; }
	blender_file << "\n";

	blender_file.close();

	//Saved coe for further polishing

	//system('move'+' ' + file_name.c_str() + ' ' + 'C:\\U'+'sers' + '\\ - ' + '\\Des'+ 'ktop');

	/*fstream blender_file;
	blender_file.open(export_path + "\\prokaryote.py", ios::out);

	for (string line : blender_includes) { blender_file << line + "\n"; }
	blender_file << "\n";
	for (string line : point_cloud_function) { blender_file << line + "\n"; }
	blender_file << "\n";

	string object_creation = point_cloud + point_cloud_name + ", ";
	object_creation += "[";

    for (int i = 0; i < assembly.elipses.size(); i++)
    {
        for (int j = 0; j < assembly.elipses[i]->points.size(); j++)
        {
            object_creation += "(" + to_string(assembly.elipses[i]->points[j][0]) + ", " + 
				to_string(assembly.elipses[i]->points[j][1]) + ", " + to_string(assembly.elipses[i]->points[j][2]) + "), ";
        }
    }
	object_creation += "])";

	blender_file << object_creation + "\n";
	blender_file << "\n";
	blender_file << link_object + "\n";*/


	//double val = FindEndPoint(functions, value, 1, 1);

	//vector<double> points = GetSeparators(funcv, 1, 1, 10);

	//for (double i : points) { cout << i << endl; }

	//Elipse elipse = GetElipsePoints(funcv, 1, 1, 2, { 2.07944, 3, 0.693147 }, 0.5, 3);

	//for (vector<double> i : elipse.points) { cout << i[0] << ", "<< i[1] << ", " << i[2] << endl; }

	//Be wary and use 1.0 instead of 1 and similiraly, as there then divisional error (divisions by zero etc.)

}