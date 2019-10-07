// ==========================================================================
// MeshOui
//
// Copyright (c) D-ICE Engineering.
// All rights reserved.
// ==========================================================================

#ifndef MESHOUI_MATH_H
#define MESHOUI_MATH_H

#include "MathUtils/Constants.h"

#include "MathUtils/Matrix33.h"
#include "MathUtils/Matrix66.h"
#include "MathUtils/Matrix.h"

#include "MathUtils/Vector3d.h"
#include "MathUtils/Vector6d.h"
#include "MathUtils/VectorN.h"

namespace meshoui {

    // Vectors.
    using Vector3d = mathutils::Vector3d<double>;
    using Vector6d = mathutils::Vector6d<double>;
    using VectorN = mathutils::VectorN<double>;

    // Matrices.
    using Matrix33 = mathutils::Matrix33<double>;
    using Matrix66 = mathutils::Matrix66<double>;
    using MatrixMN = mathutils::MatrixMN<double>;

}

#endif // MESHOUI_MATH_H
