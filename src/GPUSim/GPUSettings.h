//
// Created by pbahr on 22/06/2022.
//

#ifndef MOLFLOW_PROJ_GPUSETTINGS_H
#define MOLFLOW_PROJ_GPUSETTINGS_H

namespace flowgpu {
    struct MolflowGlobal {
        MolflowGlobal() = default;

        float gasMass;
        bool useMaxwellDistribution;
        size_t recursiveMaxDepth;
        size_t cyclesRNG{1};
        bool randomNumberMethod; /*! 0=bulked, 1=ad hoc */
        /*bool	 lowFluxMode;
        double	 lowFluxCutoff;*/

        float offsetMagnitude{0.0f}; // adaptive offset towards center
        float offsetMagnitudeN{1.0f}; // adaptive offset towards normal
    };
}
#endif //MOLFLOW_PROJ_GPUSETTINGS_H
