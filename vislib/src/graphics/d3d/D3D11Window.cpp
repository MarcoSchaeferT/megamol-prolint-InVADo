/*
 * D3D11Window.cpp
 *
 * Copyright (C) 2006 - 2012 by Visualisierungsinstitut Universitaet Stuttgart.
 * Alle Rechte vorbehalten.
 */

#include "vislib/graphics/d3d/D3D11Window.h"

#include "vislib/graphics/d3d/D3DException.h"
#include "vislib/Trace.h"


/*
 * vislib::graphics::d3d::D3D11Window::~D3D11Window
 */
vislib::graphics::d3d::D3D11Window::~D3D11Window(void) {
}


/*
 * vislib::graphics::d3d::D3D11Window::onCreated
 */
void vislib::graphics::d3d::D3D11Window::onCreated(HWND hWnd) {
    AbstractWindow::onCreated(hWnd);
    this->initialise(hWnd);
}


/*
 * vislib::graphics::d3d::D3D11Window::onResized
 */
void vislib::graphics::d3d::D3D11Window::onResized(const int width, 
        const int height) {
    this->resizeSwapChain(width, height);
}
