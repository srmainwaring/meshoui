// ==========================================================================
// MeshOui
//
// Copyright (c) D-ICE Engineering.
// All rights reserved.
// ==========================================================================

#ifndef MESHOUI_MESHOUI_H
#define MESHOUI_MESHOUI_H

#include "maths.h"
#include "mesh.h"
#include "integration.h"
#include "integrand_on_face.h"
#include "clipping_surfaces.h"
#include "clipper.h"
#include "plane.h"
#include "polygon.h"
#include "version.h"

#ifdef USE_VTK
#include "vtkmesh.h"
#endif

#ifdef USE_MMG
#include "remesher.h"
#endif

#endif // MESHOUI_MESHOUI_H
