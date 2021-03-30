// ==========================================================================
// MeshOui
//
// Copyright (c) D-ICE Engineering.
// All rights reserved.
// ==========================================================================

#ifndef MESHOUI_INTEGRAND_ON_FACE_H
#define MESHOUI_INTEGRAND_ON_FACE_H

#include "MathUtils/Integrand.h"

namespace meshoui {

  /**
  * Class for handling the function (integrand) to be integrated on a face.
  */
  template<typename T>
  class IntegrandOnFace : public mathutils::Integrand<T> {

   public:

    /// This function evalutes the integrand at the point x.
    virtual T Evaluate(const Vector3d &x) const = 0;

  };


}

#endif //MESHOUI_INTEGRAND_ON_FACE_H
