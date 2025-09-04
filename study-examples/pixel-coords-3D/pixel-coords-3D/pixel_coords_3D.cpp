//[header]
// This program renders a wireframe image of a 3D object whose description is stored
// in the program itself. The result is stored to a SVG file. To draw an image
// of each triangle making up that object, we project the vertices of each triangle
// onto the screen using perspective projection (effectively transforming the vertices
// world coordinates to 2D pixel coordinates). Triangles are stored in the SVG
// file by connecting their respective vertices to each other with lines.
//[/header]
//[compile]
// Download the pixel_coords_3D.cpp and geometry.h files to the same folder.
// Open a shell/terminal, and run the following command where the files are saved:
//
// c++ pixel_coords_3D.cpp  -o pixel_coords_3D -std=c++11
//
// Run with: ./pixel_coords_3D. Open the file ./proj.svg in any Internet browser to see
// the result.
//[/compile]

#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <cmath>

#include "geometry.h"

const Vec3f = verts[] = {} // a list of vertices making up the object.
const uint32_t numTris = 10; //???