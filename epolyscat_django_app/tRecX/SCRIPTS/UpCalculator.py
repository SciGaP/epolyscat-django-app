#!/usr/bin/env python

# tRecX = tSurff+irECS - a universal Schroedinger solver
# Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
# 
# This program is free software; you can redistribute it and/or modify it 
# under the terms of the GNU General Public License as published by the Free Software Foundation; 
# either version 2 of the License, or (at your option) any later version.
# End of license
 
# -*- coding: utf-8 -*-


from tools import *
from units import *

from plotTools import *

(args,flags) = argsAndFlags();

if len(args)<2:
    print "Usage: "
    print "   UpCalculator.py lambda(nm) Intensity(W/cm2)"
    sys.exit(0)

Units(au())

lambdaNM=float(args[0])
intenWcm2=float(args[1])


print 'U_p (eV):', au(eV()).energy((SI().intensity(intenWcm2*1e4)/nm_inv().energy(1/lambdaNM)**2)/4.)
