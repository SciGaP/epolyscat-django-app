// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "timeCritical.h"
namespace timeCritical{

static bool _on=false;
static bool _onOld=false;

bool coefficients=false;
void setOn(){_on=true;_onOld=true;      coefficients=_on;}
void setOff(){_on=false;_onOld=false;   coefficients=_on;}
void suspend(){_onOld=_on;_on=false;    coefficients=_on;}
void resume(){_on=_onOld;               coefficients=_on;}
bool isOn(){return _on;}
}
