// add your includes to this file instead of to individual .cpp files
// to enjoy the benefits of precompiled headers:
// - fast compilation
// - solve issues with the order of header files once (here)
// do not include headers in header files (ever).

// Prevent expansion clashes (when using std::min and std::max):
#define NOMINMAX

#define DEBUG_PERFORMANCE

#define SCRWIDTH 1440
#define SCRHEIGHT 810
// #define FULLSCREEN
// #define ADVANCEDGL	// faster if your system supports it

// Glew should be included first
#include <GL/glew.h>
// Comment for autoformatters: prevent reordering these two.
#include <GL/gl.h>

#ifdef _WIN32
// Followed by the Windows header
#include <Windows.h>

// Then import wglext: This library tries to include the Windows
// header WIN32_LEAN_AND_MEAN, unless it was already imported.
#include <GL/wglext.h>

// Extra definitions for redirectIO
#include <fcntl.h>
#include <io.h>
#endif

// External dependencies:
#include <FreeImage.h>
#include <SDL2/SDL.h>

// C++ headers
#include <algorithm>
#include <chrono>
#include <fstream>
#include <iostream>
#include <memory>
#include <random>
#include <vector>

// Namespaced C headers:
#include <cassert>
#include <cinttypes>
#include <cmath>
#include <cstdio>
#include <cstdlib>

// Header for AVX, and every technology before it.
// If your CPU does not support this, include the appropriate header instead.
// See: https://stackoverflow.com/a/11228864/2844473
#include <immintrin.h>

//imgui lib
#include "lib/imgui/imgui.h"
#include "lib/imgui/examples/imgui_impl_sdl.h"
#include "lib/imgui/imgui_sdl.h"


// clang-format off

// "Leak" common namespaces to all compilation units. This is not standard
// C++ practice but a mere simplification for this small project.
using namespace std;

#include "Bucket.h"
#include "BucketPool.h"

#include "Statistics.h"

#include "template.h"
#include "surface.h"

using namespace Tmpl8;
#include "Graph.h"
#include "swarmz.h"

#include "Scenario.h"


#include "game.h"
// clang-format on
