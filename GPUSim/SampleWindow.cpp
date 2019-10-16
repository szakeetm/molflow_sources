//
// Created by pbahr on 08/10/2019.
//

#include "SampleWindow.h"

void osc::SampleWindow::render(){
    if (cameraFrame.modified) {
        sample.setCamera(Camera{ cameraFrame.get_from(),
                                 cameraFrame.get_at(),
                                 cameraFrame.get_up() });
        cameraFrame.modified = false;
    }
    sample.render();
}

void osc::SampleWindow::draw(){
    sample.downloadPixels(pixels.data());
    if (fbTexture == 0)
        glGenTextures(1, &fbTexture);

    glBindTexture(GL_TEXTURE_2D, fbTexture);
    GLenum texFormat = GL_RGBA;
    GLenum texelType = GL_UNSIGNED_BYTE;
    glTexImage2D(GL_TEXTURE_2D, 0, texFormat, fbSize.x, fbSize.y, 0, GL_RGBA,
                 texelType, pixels.data());

    glDisable(GL_LIGHTING);
    glColor3f(1, 1, 1);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, fbTexture);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

    glDisable(GL_DEPTH_TEST);

    glViewport(0, 0, fbSize.x, fbSize.y);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0.f, (float)fbSize.x, 0.f, (float)fbSize.y, -1.f, 1.f);

    glBegin(GL_QUADS);
    {
        glTexCoord2f(0.f, 0.f);
        glVertex3f(0.f, 0.f, 0.f);

        glTexCoord2f(0.f, 1.f);
        glVertex3f(0.f, (float)fbSize.y, 0.f);

        glTexCoord2f(1.f, 1.f);
        glVertex3f((float)fbSize.x, (float)fbSize.y, 0.f);

        glTexCoord2f(1.f, 0.f);
        glVertex3f((float)fbSize.x, 0.f, 0.f);
    }
    glEnd();
}

//------------------------------------------------------------------------------
//
// GLFW callbacks
//
//------------------------------------------------------------------------------

/*! callback for _moving_ the mouse to a new position */
static void glfwindow_mouseMotion_cb(GLFWwindow *window, double x, double y)
{
    osc::SampleWindow *gw = static_cast<osc::SampleWindow*>(glfwGetWindowUserPointer(window));
    assert(gw);
    gw->mouseMotion(gdt::vec2i((int)x, (int)y));
}

/*! callback for a key press */
static void glfwindow_key_cb(GLFWwindow *window, int key, int scancode, int action, int mods)
{
    osc::SampleWindow *gw = static_cast<osc::SampleWindow*>(glfwGetWindowUserPointer(window));
    assert(gw);
    if (action == GLFW_PRESS) {
        gw->key(key,mods);
    }
}

/*! callback for pressing _or_ releasing a mouse button*/
static void glfwindow_mouseButton_cb(GLFWwindow *window, int button, int action, int mods)
{
    osc::SampleWindow *gw = static_cast<osc::SampleWindow*>(glfwGetWindowUserPointer(window));
    assert(gw);
    // double x, y;
    // glfwGetCursorPos(window,&x,&y);
    gw->mouseButton(button,action,mods);
}

void osc::SampleWindow::run()
{
    int width, height;
    glfwGetFramebufferSize(handle, &width, &height);
    resize(vec2i(width,height));

    // glfwSetWindowUserPointer(window, GLFWindow::current);
    //glfwSetFramebufferSizeCallback(handle, glfwindow_reshape_cb);
    //glfwSetMouseButtonCallback(handle, glfwindow_mouseButton_cb);
    //glfwSetKeyCallback(handle, glfwindow_key_cb);
    //glfwSetCursorPosCallback(handle, glfwindow_mouseMotion_cb);

    //glfwSetMouseButtonCallback( handle, mouseButtonCallback );
    //glfwSetKeyCallback        ( handle, keyCallback         );
    //glfwSetCursorPosCallback  ( handle, cursorPosCallback   );

    //while (!glfwWindowShouldClose(handle)) {
        render();
        draw();

        glfwSwapBuffers(handle);
        glfwPollEvents();
    //}
    std::cout << "#osc: camera at " << cameraFrame.position << std::endl;

}

void osc::SampleWindow::resize(const vec2i &newSize){
    fbSize = newSize;
    sample.resize(newSize);
    pixels.resize(newSize.x*newSize.y);
}