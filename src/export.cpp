
#include "export.h"

//Values to export in object blender creation

vector<string> blender_includes = { "import bpy", "import math" };

vector<string> setting_scene = { "#Set the viewing scene for the displacement export",
"bpy.data.scenes['Scene'].render.engine = 'CYCLES'",
"bpy.data.scenes['Scene'].cycles.feature_set = 'EXPERIMENTAL'" };

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
"bpy.ops.object.mode_set(mode = 'OBJECT')" };

vector<string> connect_elipses = {
	"bpy.ops.object.mode_set(mode = 'EDIT')",
	"bpy.ops.mesh.bridge_edge_loops()",
	"bpy.ops.mesh.select_all(action = 'DESELECT')",
	"bpy.ops.object.mode_set(mode = 'OBJECT')"};

vector<string> vertices_merge = { "bpy.ops.object.mode_set(mode = 'EDIT')",
"bpy.ops.mesh.merge(type = 'CENTER')",
"bpy.ops.mesh.select_all(action = 'DESELECT')",
"bpy.ops.object.mode_set(mode = 'OBJECT')" };

vector<string> object_smooth = { 
"obj.select_set(True)",
"obj.data.use_auto_smooth = 1",
"obj.data.auto_smooth_angle = math.pi / 6",
"bpy.ops.object.shade_smooth()" };

vector<string> prepare_mesh = { "#Prepare the mesh for matereial addition",
"bpy.ops.mesh.uv_texture_add()",
"#Add Subdivision modifier",
"bpy.ops.object.modifier_add(type = 'SUBSURF')",
"obj.modifiers['Subdivision'].levels = 3",
"obj.modifiers['Subdivision'].render_levels = 3"};

vector<string> add_material = { "#Add the material",
"mat = bpy.data.materials.new('prokaryote_material')",
"mat.use_nodes = True",
"mat.cycles.displacement_method = 'BOTH'",
"obj.data.materials.append(mat)" };

vector<string> set_nodes = { "#Get the current nodesand links of the added material",
"context = bpy.context",
"nodes = context.active_object.active_material.node_tree.nodes",
"links = context.active_object.active_material.node_tree.links",
"base_node = nodes.get('Principled BSDF')",
"material_output = nodes.get('Material Output')" };

vector<string> additional_nodes = { 
"# Add textue coordinate node",
"texture_coordinate = nodes.new('ShaderNodeTexCoord')",
"texture_coordinate.location = (-1100, 100)",
"texture_coordinate.label = 'Texture Coordinate'",
"# Add mapping node",
"mapping = nodes.new('ShaderNodeMapping')",
"mapping.location = (-900, 100)",
"mapping.label = 'Mapping'" };

vector<string> link_nodes = { "#Link the mapping to all input nodes",
"links.new(texture_coordinate.outputs[2], mapping.inputs[0])",
"links.new(mapping.outputs[0], base_color.inputs[0])",
"links.new(mapping.outputs[0], roughness.inputs[0])",
"#Link the individual nodes",
"links.new(base_color.outputs[0], base_node.inputs[0])",
"links.new(roughness.outputs[0], base_node.inputs[9])" };



string point_cloud_comment = "#Create the collection of points, lines and planes";
string point_cloud = "pc = point_cloud(";

string point_cloud_name = "'prokaryote_body'";

string link_object_comment = "#Link the the collection to the scene";
string link_object = "bpy.context.collection.objects.link(pc)";

std::wstring ExePath() {
	TCHAR buffer[MAX_PATH] = { 0 };
	GetModuleFileName(NULL, buffer, MAX_PATH);
	std::wstring::size_type pos = std::wstring(buffer).find_last_of(L"\\/");
	return std::wstring(buffer).substr(0, pos);
}

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

vector<string> CreateShaderNode(int x, int y, string name, string path)
{
	string comment = "# Add " + name + " node";
	string node = name + " = nodes.new('ShaderNodeTexImage')";
	string location = name + ".location = (" + to_string(x) + ", " + to_string(y) + ")";
	string label = name + ".label = '" + name + "'";
	string image = name + ".image = bpy.data.images.load('" + path + "')";

	return { comment, node, location, label, image };
}

void ExportToBlender(ProkaryoteBodyContainer assembly, string export_path, string file_name, double displace, double resolution, bool quadriflow, string displace_type, int dimensions, int r_color, int g_color, int b_color)
{
	ofstream blender_file(file_name);

	for (string line : blender_includes) { blender_file << line + "\n"; }
	blender_file << "\n";
	for (string line : setting_scene) { blender_file << line + "\n"; }
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

	object_creation += "])";

	blender_file << object_creation + "\n" << "\n";
	blender_file << link_object_comment << "\n";
	blender_file << link_object + "\n" << "\n";

	for (string line : select_object) { blender_file << line + "\n"; }
	blender_file << "\n";

	blender_file << "#Connect individual elipses" << "\n";
	blender_file << "counter = 0" << "\n";
	string for_loop_1 = "for i in range(0, " + to_string(assembly.elipses.size()-1) + "):";
	blender_file << for_loop_1 << "\n";
	blender_file << "	for j in range(counter," << 2*resolution << " + counter):" << "\n";
	blender_file << "		obj.data.edges[j].select = True" << "\n";
	blender_file << "	counter += " << resolution << "\n";
	for (string line : connect_elipses) { blender_file << line + "\n"; }
	blender_file << "\n";

	for (string line : select_object) { blender_file << line + "\n"; }
	blender_file << "\n";

	blender_file << "#Select all points on one side of cap and merge them" << "\n";
	string for_loop_2 = "for i in range(0, " + to_string(assembly.elipses[0]->points.size())  +"):";
	blender_file << for_loop_2 + "\n";
	blender_file << "    obj.data.vertices[i].select = True" + std::string("\n");

	for (string line : vertices_merge) { blender_file << line + "\n"; }
	blender_file << "\n";

	blender_file << "#Select all points on the other side and merge them" << "\n";
	for_loop_2 = "for i in range(" + to_string(int(assembly.elipses.size() * assembly.elipses[assembly.elipses.size() - 1]->points.size() - assembly.elipses[0]->points.size() - resolution + 1)) + "," + to_string(int(assembly.elipses.size() * assembly.elipses[assembly.elipses.size() - 1]->points.size() - resolution + 1)) + "):";
	blender_file << for_loop_2 + "\n";
	blender_file << "    obj.data.vertices[i].select = True" + std::string("\n");

	for (string line : vertices_merge) { blender_file << line + "\n"; }
	blender_file << "\n";

	blender_file << "#Smooth the surface" << "\n";
	blender_file << "bpy.ops.object.voxel_remesh()" << "\n";
	if (quadriflow == true) { blender_file << "bpy.ops.object.quadriflow_remesh()" << "\n"; }

	for (string line : object_smooth) { blender_file << line + "\n"; }
	blender_file << "\n";

	//Create materail files (displacement, roughness and base color)
	auto wstr = ExePath();

	std::string current_directory(wstr.begin(), wstr.end());
	current_directory.erase(current_directory.length() - 9);

	std::remove((current_directory + (string)"\\TestingMath\\MaterialCreation\\generated_maps\\base_color_map.jpg").c_str());
	std::remove((current_directory + (string)"\\TestingMath\\MaterialCreation\\generated_maps\\roughness_map.jpg").c_str());
	std::remove((current_directory + (string)"\\TestingMath\\MaterialCreation\\generated_maps\\displacement_map.jpg").c_str());

	//System Path
	string displacement = (string)"cd " + current_directory + (string)"\\TestingMath\\MaterialCreation && python3 displacement.py " + to_string(dimensions);
	string roughness = (string)"cd " + current_directory + (string)"\\TestingMath\\MaterialCreation && python3 roughness.py " + to_string(dimensions);
	string base_color = (string)"cd " + current_directory + (string)"\\TestingMath\\MaterialCreation && python3 base_color.py " + 
		to_string(dimensions) + " " + to_string(r_color) + " " + to_string(g_color) + " " + to_string(b_color);

	//Convert to const char pointer
	const char* pointer_displacement = (displacement).c_str(); 
	const char* pointer_roughness = (roughness).c_str(); 
	const char* pointer_base_color = (base_color).c_str();
	
	std::system(pointer_displacement); std::system(pointer_roughness); std::system(pointer_base_color);

	//Prepare the blender for material addition
	for (string line : prepare_mesh) { blender_file << line + "\n"; }
	blender_file << "\n";

	//Adding the displacement 
	blender_file << "#Add Displace modifier" << "\n";
	blender_file << "bpy.ops.object.modifier_add(type = 'DISPLACE')" << "\n";
	blender_file << "obj.modifiers['Displace'].strength = " << displace << "\n";
	blender_file << "obj.modifiers['Displace'].mid_level = 0.000" << "\n";
	if (displace_type == "IMAGE") {
		blender_file << "texture = bpy.data.textures.new('Image', 'IMAGE')" << "\n";
		blender_file << "texture.image = bpy.data.images.load(" << current_directory + (string)"/TestingMath/MaterialCreation/generated_maps/displacement_map.jpg" << ")" << "\n";
		blender_file << "obj.modifiers['Displace'].texture = bpy.data.textures['Image']" << "\n";
	}
	else if (displace_type == "CLOUDS") {
		blender_file << "texture = bpy.data.textures.new('Cloud', 'CLOUDS')" << "\n";
		blender_file << "obj.modifiers['Displace'].texture = bpy.data.textures['Cloud']" << "\n";
	}

	//Add the material
	for (string line : add_material) { blender_file << line + "\n"; }
	blender_file << "\n";

	//Prepare for shading
	for (string line : set_nodes) { blender_file << line + "\n"; }
	blender_file << "\n";

	//Add the individual nodes
	//base color
	string base_color_path = current_directory + (string)"/TestingMath/MaterialCreation/generated_maps/base_color_map.jpg";
	std::replace(base_color_path.begin(), base_color_path.end(), '\\', '/');
	for (string line : CreateShaderNode(-600, 210, "base_color", base_color_path)) { blender_file << line + "\n"; }
	blender_file << "\n";
	//roughness
	string roughness_path = current_directory + (string)"/TestingMath/MaterialCreation/generated_maps/roughness_map.jpg";
	std::replace(roughness_path.begin(), roughness_path.end(), '\\', '/');
	for (string line : CreateShaderNode(-600, -80, "roughness", roughness_path)) { blender_file << line + "\n"; }
	blender_file << "\n";

	for (string line : additional_nodes) { blender_file << line + "\n"; }
	blender_file << "\n";

	for (string line : link_nodes) { blender_file << line + "\n"; }
	blender_file << "\n";

	blender_file.close();

}