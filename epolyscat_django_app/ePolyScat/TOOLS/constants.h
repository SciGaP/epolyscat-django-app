// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef CONSTANTS_H
#define CONSTANTS_H

/// \ingroup Units
/// \brief some (very few) math constants
namespace math {
    const double pi=3.141592653589793238;           //!< pi with 19 correct digits
}

/// \ingroup Units
/// \brief the most important physics constants
namespace physics {

	const double a_finestructure = 1.0 / 137.035999679;							//!< alpha 2009 (NIST)
	const double electron_mass = 9.10938215e-31;								//!< m_e [kg] 2009 (NIST)
	const double bohr_radius = 5.2917720859e-11;								//!< a_0 [m]  2009 (NIST)
    const double speed_of_light = 299792458.0;									//!< c_0 [m/s]
//    const double speed_of_light_in_au = 137.0359990738557;                      //!< c_0 [a.u.] == 1/a_finestructure
    const double vacuum_permeability = 4.*math::pi*1.e-7;						//!< µ_0 [H/m]
	const double mu0 = vacuum_permeability;										//!< shorthand µ_0
    const double vacuum_permittivity = 1/(mu0*pow(speed_of_light,2));			//!< eps_0 [As/Vm]
	const double eps0 = vacuum_permittivity;									//!< shorthand eps_0
	const double proton_charge = 1.602176487e-19;								//!< e   [C]  2009 (NIST)
	const double planck_constant = 6.62606896e-34;								//!< h   [Js] 2009 (NIST)
    const double h_bar = planck_constant / (2.*math::pi);						//!< h_bar [Js]
    const double hartree = std::pow(speed_of_light*a_finestructure, 2)*electron_mass;//!< E_H [J]
    const double avogadro = 6.02214076e23; //!< 1/mol

//	const static double mc2 = 0.510998928e6; ///< electron energy mass*c^2 in eV
//	const static double h_bar_c = 197.3269718e-9; ///< hbar c in eV*m

	
//	const static double e = 1.602176565e-19; ///< e in As
//	const static double me = 9.10938188e-31; ///< electron mass in kg

//	const static double c0inAU = 137.0359990738557;			// speed of ligh

}
#endif // CONSTANTS_H
