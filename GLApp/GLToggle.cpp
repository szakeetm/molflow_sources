/*
  File:        GLToggle.cpp
  Description: Toggle Button class (SDL/OpenGL OpenGL application framework)
  Author:      J-L PONS (2007)

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/
#include "GLWindow.h"
#include "GLToggle.h"
#include "GLToolkit.h"

// ---------------------------------------------------------------------

GLToggle::GLToggle(int compId,char *text):GLComponent(compId) {
  if(text) 
    strcpy(this->text,text);
  else
    strcpy(this->text,"");
  state=0;
  SetBorder(BORDER_NONE);
  SetTextColor(0,0,0);
}

// ---------------------------------------------------------------------

BOOL GLToggle::IsChecked() {
  return state;
}

// ---------------------------------------------------------------------

void GLToggle::SetCheck(BOOL checked) {
  state=checked;
}

// ---------------------------------------------------------------------

void GLToggle::Paint() {

  if(!parent) return;
  GLFont2D *font = GLToolkit::GetDialogFont();

  GLComponent::Paint();


  font->SetTextColor(rText,gText,bText);
  font->DrawText(posX+16,posY+2,text,FALSE);

  GLToolkit::DrawToggle(posX+2,posY+3);
  if(state) {
    font->SetTextColor(0.0f,0.0f,0.0f);
    font->DrawText(posX+5,posY+1,"\215",FALSE);
  }

}

// ---------------------------------------------------------------------

void GLToggle::SetTextColor(int r,int g,int b) {
  rText = r/255.0f;
  gText = g/255.0f;
  bText = b/255.0f;
}

// ---------------------------------------------------------------------

void GLToggle::ManageEvent(SDL_Event *evt) {

  if(!parent) return;

  int x,y,w,h;
  GetWindow()->GetBounds(&x,&y,&w,&h);

  if( evt->type == SDL_MOUSEBUTTONDOWN ) {
    if( evt->button.button == SDL_BUTTON_LEFT ) {
      state = !state;
      parent->ProcessMessage(this,MSG_TOGGLE);
    }
  }

}
