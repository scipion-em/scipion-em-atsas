# *****************************************************************************
# *
# * Authors:     Carlos Oscar Sorzano (coss@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# *****************************************************************************

import os
import pwem
from pyworkflow.utils import Environ

from .constants import *

__version__ = '3.1'
_references = ['Manalastas-Cantos2021']
_logo = "atsas_logo.gif"


class Plugin(pwem.Plugin):
    _homeVar = ATSAS_HOME
    _pathVars = [ATSAS_HOME]
    _supportedVersions = [V3_0_1, V3_2_0]
    _url = "https://www.embl-hamburg.de/biosaxs/software.html"

    @classmethod
    def _defineVariables(cls):
        cls._defineEmVar(ATSAS_HOME,
                         f"atsas-{cls._supportedVersions[-1]}")

    @classmethod
    def getEnviron(cls):
        environ = Environ(os.environ)
        environ.update({
            'ATSAS': cls.getHome(),
            'PATH': cls.getHome('bin'),
            'LD_LIBRARY_PATH': cls.getHome('lib')
        }, position=Environ.BEGIN)

        return environ

    @classmethod
    def defineBinaries(cls, env):
        pass

    @classmethod
    def getProgram(cls, binary='crysol'):
        """ Return the program binary that will be used. """
        return cls.getHome('bin', binary)
