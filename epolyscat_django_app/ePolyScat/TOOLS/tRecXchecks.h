// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef _TRECXCHECKS_H_
#define _TRECXCHECKS_H_

#include <vector>
#include <string>
#include <iostream>


/** @defgroup Checks
 *  @ingroup Tools
  * \brief global control of checks and verifications
  * @{
  */
/// \file

namespace tRecX{

/// \ingroup Checks
/// \fn for globally switching on- and off checks
/// \brief usage: off("someName") returns false, if "someName" is in Checks: off=someName someOther yetSome
bool off(const std::string mess="global");
void read(); // create the info
void print(); // display status
void setOff(std::string Off);

}
/** @} */

#endif
