//
// View3DMouse.cpp
//
// Copyright (C) 2012 by University of Stuttgart (VISUS).
// Copyright (C) 2019 by University of Tuebingen (BDVA).
// All rights reserved.
//

#include "stdafx.h"
#define _USE_MATH_DEFINES

#include "View3DMouse.h"
#include "mmcore/factories/CallAutoDescription.h"
#include "MouseInputCall.h"
#include "mmcore/param/ButtonParam.h"
#include "mmcore/view/MouseFlags.h"


using namespace megamol;


prolint::View3DMouse::View3DMouse(void) : core::view::View3D(),
         mouseSlot("mouse", "Slot to send mouse information to the renderer."),
         enableSelectingSlot("enableSelecting", "Enable selecting and picking with the mouse."),
         toggleSelect(true) {

    // Slot for mouse input call
	this->mouseSlot.SetCompatibleCall<core::factories::CallAutoDescription<prolint::MouseInputCall> >();
    this->MakeSlotAvailable(&this->mouseSlot);

    // Slot for key modifier
    this->enableSelectingSlot << new core::param::ButtonParam(core::view::Key::KEY_TAB);
    this->enableSelectingSlot.SetUpdateCallback(&View3DMouse::OnButtonEvent);
    this->MakeSlotAvailable(&this->enableSelectingSlot);
}

using namespace megamol::core;
#include "mmcore/view/CallRender3D.h"
#ifdef tester
bool prolint::View3DMouse::MouseEvent(int x, int y, core::view::MouseFlags flags) { return true; }

//#######

#endif
bool prolint::View3DMouse::OnMouseButton(
    core::view::MouseButton button, core::view::MouseButtonAction action, core::view::Modifiers mods) {
    /*
	if (mods.test(core::view::Modifier::CTRL)) {
        // These clicks go to the view.
        printf("CTRL click\n");
        return false;
    }
	*/
    auto* cr = this->rendererSlot.CallAs<view::CallRender3D>();
    if (cr) {
        view::InputEvent evt;
        evt.tag = view::InputEvent::Tag::MouseButton;
        evt.mouseButtonData.button = button;
        evt.mouseButtonData.action = action;
        evt.mouseButtonData.mods = mods;
        cr->SetInputEvent(evt);
        if ((*cr)(view::CallRender3D::FnOnMouseButton)) return true;
    }
	
    auto down = action == view::MouseButtonAction::PRESS;
    if (mods.test(view::Modifier::SHIFT)) {
        this->modkeys.SetModifierState(vislib::graphics::InputModifiers::MODIFIER_SHIFT, down);
    } else if (mods.test(view::Modifier::CTRL)) {
        this->modkeys.SetModifierState(vislib::graphics::InputModifiers::MODIFIER_CTRL, down);
    } else if (mods.test(view::Modifier::ALT)) {
        this->modkeys.SetModifierState(vislib::graphics::InputModifiers::MODIFIER_ALT, down);
    }

	
    switch (button) {
    case view::MouseButton::BUTTON_LEFT:
        this->cursor2d.SetButtonState(0, down);
        if (down == (int)view::MouseButtonAction::RELEASE) {
            this->mouseFlags = 0;
        } else {
            this->mouseFlags = -1;
		}

        break;
    case view::MouseButton::BUTTON_RIGHT:
        this->cursor2d.SetButtonState(1, down);
        this->mouseFlags = 1;
        break;
    case view::MouseButton::BUTTON_MIDDLE:
        this->cursor2d.SetButtonState(2, down);
        this->mouseFlags = 2;
        break;
    default:
        break;
    }

    SetCursor2DPosition(this->mouseX, this->mouseY);
    
    return false;
}

//#######
bool prolint::View3DMouse::OnMouseMove(double x, double y) {
   
	this->mouseFlags = -1;
    this->mouseX = static_cast<float>(x);
    this->mouseY = static_cast<float>(y);

   auto* cr = this->rendererSlot.CallAs<view::CallRender3D>();
    if (cr) {
        view::InputEvent evt;
        evt.tag = view::InputEvent::Tag::MouseMove;
        evt.mouseMoveData.x = x;
        evt.mouseMoveData.y = y;
        cr->SetInputEvent(evt);
        if ((*cr)(view::CallRender3D::FnOnMouseMove)) return true;
    }

    this->cursor2d.SetPosition(mouseX, mouseY, true);
    SetCursor2DPosition(mouseX, mouseY);
   
    return true;
}



void prolint::View3DMouse::SetCursor2DButtonState(unsigned int btn, bool down) {

    //if(!this->toggleSelect) {
    //    View3D::SetCursor2DButtonState(btn, down); // Keep camera movement functional
    //}
    //else {
        switch (btn) {
        case 0 : // left
            megamol::core::view::MouseFlagsSetFlag(this->mouseFlags,
                    core::view::MOUSEFLAG_BUTTON_LEFT_DOWN, down);
            break;
        case 1 : // right
            megamol::core::view::MouseFlagsSetFlag(this->mouseFlags,
                    core::view::MOUSEFLAG_BUTTON_RIGHT_DOWN, down);
            break;
        case 2 : // middle
            megamol::core::view::MouseFlagsSetFlag(this->mouseFlags,
                    core::view::MOUSEFLAG_BUTTON_MIDDLE_DOWN, down);
            break;
        }
        //}
}


void prolint::View3DMouse::SetCursor2DPosition(float x, float y) {

    if(!this->toggleSelect) {
    //    View3D::SetCursor2DPosition(x, y); // Keep camera movement functional
       
        core::view::MouseFlagsResetAllChanged(this->mouseFlags);
    }
    else {
       
		prolint::MouseInputCall *cm = this->mouseSlot.CallAs<prolint::MouseInputCall>();
        if (cm) {
            cm->SetCam(cam);
            cm->SetMouseInfo(static_cast<int>(x), static_cast<int>(y), this->mouseFlags);
            if ((*cm)(0)) {
                this->mouseX = (float)static_cast<int>(x);
                this->mouseY = (float)static_cast<int>(y);
                megamol::core::view::MouseFlagsResetAllChanged(this->mouseFlags);
                // mouse event consumed
                return;
            }
            megamol::core::view::MouseFlagsResetAllChanged(this->mouseFlags);
        }
    }
}


prolint::View3DMouse::~View3DMouse() {
	// ToDo
}


bool prolint::View3DMouse::OnButtonEvent(core::param::ParamSlot& p) {
    printf("TAB pressed: \n"); // DEBUG

	//printf("x: %f\n", this->cursor2d.X());
	
	this->toggleSelect = !this->toggleSelect;
    if (!this->toggleSelect) {
        printf("mouse picking / hovering: OFF \n"); 
	} else {
        printf("mouse picking / hovering: ON \n");
	}

    return true;
}


