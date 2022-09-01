# Simulating evolution
This is a longreaching project trying to simulate, to some extent, evolution of organisms. The project is currently in the begginning of demostration phase, which aim is to present the main idea in the form of prokaryotic evolution. The finished product should run on UE5 and so is mainly written in c++.

With this having been said, the main purpose of this repository is to post and share main parts of the code, which are generally not directly linked to UE5, for discussion and improvement. As poreviously stated, this project was only started recently and so only one part is fully finished: generating prokaryote main body from gene for blender (and there can still be some undetected errors in there). The main focus is currently on full body design, with motion and entatling physics, with them being added as sonn as finished. Both the readme and source code is only temporali and will be frequently changed with coming changes.

The full vision of the project will be soon presented in scientific paper, to which link will be provided here.

## Testing Prokaryote body creation
How to use and analyze the current body generation functions is described in src under TestingExport.cpp file. The generated file is named test1.txt and should be present in the same directory as the TestingExport.cpp file (one generated example also present in src directory).

The text from this file can then be copied and pasted into blender's scripting workspace, where after hitting run the body is generated and can be exported. Currently, the result is not smoothed, but it can be done by adding shade smooth and clicking auto smooth in normals.

If you find any error, contacting me would be very appreciated.
