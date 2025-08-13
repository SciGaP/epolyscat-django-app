// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef UNITS_H
#define UNITS_H
#include "toolsHeader.h"

/** @defgroup Units Units and constants
 *  @ingroup Tools
  * \brief Physical, mathematical, chemical constants, unit conversion
  * @{
  */

/// @brief convert between units
///
/// convert to full precisions as available at NIST (see constants.h)
///
/// values_in_new_units=Units::convert(values_in_input_units,input_unit_name,[output_units])
///
/// Examples:
///
/// Units::convert(0.5,"au|energy","eV"): returns ~13.56 <br>
/// Units::convert(speed_of_light,"SI","au|velocity"): returns   ~1/137 <br>
/// Units::setDefault("SI"): convert to these, if no output units are specified <br>
/// Units::convert(1./137,"au|veloctiy"): returns  ~ 3.e8 (convert to default, which was set "SI") <br>
/// Units::convert(1,"au|length"): returns  ~0.5e-10 (Bohr radius in SI) <br>
class Units
{
public:
    static std::string sep; ///< separator between name and dimension, e.g. au|length, W/cm2|intensity, Bohr|length, etc.
    static std::map<std::string,std::string> systems; ///< table of available system names
    static std::map<std::string,double> uniTab;       ///<table of all available units
    static std::map<std::string,std::string> aka;     ///< map from alternate unit name to standard name

    static void setDefault(std::string Def); ///< define default units
    static void standardUnits();             ///< set up frequently used units, if not default units are defined, set default to "SI"

    /// construct a new unit system by specifying basic units
    static void addUnitSystem(std::string System, /**< unit system name, e.g. SI,au,cgs,... */
                              double Length,/**< SI value of unit, e.g. bohr_radius */
                              double Mass,  /**< SI value of unit mass, e.g. electron_mass */
                              double Time,  /**< SI value of unit time */
                              double Charge,/**< SI value of unit charge, e.g. proton_mass */
                              double Mu0,   /**< SI value of \f$\mu_0\f$  */
                              std::string Names=""/**< list of alternate names, e.g. "charge|C,velocity|m/s" */);

    ///< add a freely defined unit
    static void addUnit(std::string NameDimension, /**< unit name and dimension e.g. "W/cm2|intensity, OpticalCycle|time" */
                        double Value               /**< value expressed SI units */);

    static bool isDefined(std::string Name); ///< true if unit name is known

    /// return Val converted from InUnits to OutUnits
    static double convert(double Val,std::string InUnits,std::string OutUnits="DEFAULT_SYSTEM");
    /// vector version (see scalar)
    static std::vector<double> convert(std::vector<double> Val,std::string InUnits,std::string OutUnits="DEFAULT_SYSTEM");
    static std::vector<std::complex<double> > convert(std::vector<std::complex<double> > Val,std::string InUnits,std::string OutUnits="DEFAULT_SYSTEM");

    /// sanity check for a few known units
    static void test();

    /// check wether 1 of Unit equals Val up to Digits digits
    static void check(const std::string Unit, double Val, int Digits);
};
/** @} */
#endif // UNITS_H
