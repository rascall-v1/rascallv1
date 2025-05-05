# -*- coding: utf-8 -*-
#!/usr/bin/env python
#
# Copyright (C) 2017 - Massachusetts Institute of Technology (MIT)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""
ATMOS_1 testing cases

"""
import numpy as np

from .xsec import ATMOS_1_Simulator

import matplotlib.pyplot as plt



def ATMOS_crosssection(molecule_smile):
    window = "earth2_atmosphere"
    pressure = 10e6
    temperature = 300

    Simulator = ATMOS_1_Simulator(molecule_smile, window, pressure, temperature)
    nu = np.arange(500, 5000, 1)
    xsec = Simulator.get_cross_section(nu)

    return nu, xsec






