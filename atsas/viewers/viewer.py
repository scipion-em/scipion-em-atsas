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

from pyworkflow.viewer import DESKTOP_TKINTER, Viewer
from pyworkflow.gui.plotter import Plotter

from ..protocols.protocol_pdb_to_saxs import AtsasProtConvertPdbToSAXS


class AtsasViewer(Viewer):
    """ Wrapper to visualize Pdb to SAXS. """
    _targets = [AtsasProtConvertPdbToSAXS]
    _environments = [DESKTOP_TKINTER]

    def __init__(self, **args):
        Viewer.__init__(self, **args)

    def visualize(self, obj, **args):
        cls = type(obj)
        if issubclass(cls, AtsasProtConvertPdbToSAXS):
            if obj.experimentalSAXS.empty():
                fnInt = obj._getExtraPath("pseudoatoms.int")
            else:
                fnInt = obj._getExtraPath("pseudoatoms_experimental_SAXS_curve.fit")

            import numpy
            x = numpy.loadtxt(fnInt, skiprows=1)
            xplotter = Plotter(windowTitle="SAXS Curves")
            a = xplotter.createSubPlot('SAXS curves', 'Angstroms^-1',
                                       'log(SAXS)', yformat=False)
            a.plot(x[:, 0], numpy.log(x[:, 1]))
            a.plot(x[:, 0], numpy.log(x[:, 3]))
            if obj.experimentalSAXS.empty():
                xplotter.showLegend(['SAXS in solution', 'SAXS in vacuum'])
            else:
                xplotter.showLegend(['Experimental SAXS', 'SAXS from volume'])
            xplotter.show()
