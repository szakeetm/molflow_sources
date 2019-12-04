//
// Created by pbahr on 08/10/2019.
//

#ifndef MOLFLOW_SAMPLEWINDOW_H
#define MOLFLOW_SAMPLEWINDOW_H

#include "OptixControl.h"

// our helper library for window handling
#include "glfWindow/GLFWindow.h"
#include <GL/gl.h>

namespace flowgpu {

    struct SampleWindow : public GLFCameraWindow {
        SampleWindow(const std::string &title,
                     const Model *model,
                     const Camera &camera,
                     const float worldScale)
                : GLFCameraWindow(title,camera.from,camera.at,camera.up,worldScale),
                  sample(model)
                  {
            sample.setCamera(camera);
        }

        virtual void render() override;

        virtual void draw() override;

        virtual void run();

        virtual void resize(const vec2i &newSize);

        vec2i                 fbSize;
        GLuint                fbTexture {0};
        OptixControl        sample;
        std::vector<uint32_t> pixels;
    };

}
#endif //MOLFLOW_SAMPLEWINDOW_H
