// ==========================================================================
// MeshOui
//
// Copyright (c) D-ICE Engineering.
// All rights reserved.
// ==========================================================================

#ifndef MESHOUI_INTEGRANDONFACE_H
#define MESHOUI_INTEGRANDONFACE_H

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

  private:

   /// Face handle to the face where the integrand is integrated.
   FaceHandle m_fh;

   /// Meshoui mesh.
   meshoui::Mesh* m_mesh;

  };


}

#endif //MESHOUI_INTEGRANDONFACE_H
