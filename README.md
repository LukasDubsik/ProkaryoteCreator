# Simulating evolution
This is a longreaching project trying to simulate, to some extent, evolution of organisms. The project is currently in the begginning of demostration phase, which aim is to present the main idea in the form of prokaryotic evolution. The finished product should run on UE5 and so is mainly written in c++.

With this having been said, the main purpose of this repository is to post and share main parts of the code, which are generally not directly linked to UE5, for discussion and improvement. As poreviously stated, this project was only started recently and so only one part is fully finished: generating prokaryote main body from gene for blender (and there can still be some undetected errors in there). The main focus is currently on full body design, with motion and entatling physics, with them being added as sonn as finished. Both the readme and source code is only temporali and will be frequently changed with coming changes.

The full vision of the project will be soon presented in scientific paper, to which link will be provided here.

## Testing Prokaryote body creation
How to use and analyze the current body generation functions is described in src under TestingExport.cpp file. The generated file is named test1.txt and should be present in the same directory as the TestingExport.cpp file (one generated example also present in src directory).

The text from this file can then be copied and pasted into blender's scripting workspace, where after hitting run the body is generated and can be exported. Currently, the result is not smoothed, but it can be done by adding shade smooth and clicking auto smooth in normals.

If you find any error, contacting me would be very appreciated.

# Udacity Review (Temporaly)
(Excuse my english as it is not my first language)

## Dependencies
The project uses only standart c++ library, so no additional instalations are necessary.

## Rubric: Temporary Udacity submission
This project should include all, at least to some extent, that was covered during udacity c++ nanodegree. Excluding the basic functionalities, like loops or work with types, project also work with file structure when outputing instructions to blender. User input is not defined by std::cin but by changing function or gene definition in TestingExport.cpp. Object oriented programming is present throughout the program, with examples of initialization, encapsulation, abstraction, composition and more. Templates are also present. 

Memory management is present throughout the functions, as both unique (used to pass vectors of elipses) and shared (mainly VectorFunction but also other) pointers are used, with locking and guarding mechanisms also there. Asynchronous programming is present only once, as by testing I determined it to be effective only there. Its function is to create elipses by their definition from gene.

One central part is creating functions during runtime. On how to define runtime functions further informations and commented code can be found in https://github.com/LukasDubsik/FunctionCreator.git where everything is explained. The project is build in the same manner as this one.

I understand that code is lacking comments, and also some parts from the rubric, but as I state above in introduction to this project, it is in early phases of development and the parts, including these requirements, are simply not finished yet. I planned to present more entitling glimpse at my work but the work took longer than expected (not that much the proggramming, but mainly the biological, mathematical and physical aspects regarding this project) and as my deadline was drawing closer I had to scramble the finished parts to have something presentable. In the near future (in month or so) I plan to finish the full body construction with movement and physics, the work that I originally planned to post, and add it here, with the project being then presentable.

## Building the project
The project includes all its cpp and header files in directory called src. I wrote this project in Visual studion 2022, and if you own it, I recommend running the project from there. It can be, of course, also run from other compiler. All that needs to be done is to just download all the files present, sort them by header and cpp, and open the TestingExport.cpp which includes the main function. All the files include necessary dependencies to communicate between themselves, and so does of course TestingExport.cpp.

In this cpp file, example of vector functions and gene is present and you could run or modify it. The result is text file whose content can be copied and pasted to blender (further instructions above). If you do not want to test the result in blender, you just need to check that no nan values are present there, if so, please contact me, the error can be in wrong function/gene definition, or in my math (that happened many, many times before to the point of near insanity) and I would try to fix it. Or if you want to check in blender, the error can be found by the result looking "weird". For comparison, one generated bacteria is present in src directory.
