# *****************************************************************************
# *
# * Authors:  Carlos Oscar Sanchez Sorzano (coss@cnb.csic.es)
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
import math
import numpy

from pwem.protocols import ProtPreprocessVolumes, PointerParam
from pyworkflow.protocol import (IntParam, FloatParam, FileParam,
                                 StringParam, LEVEL_ADVANCED)
import pyworkflow.utils as pwutils

from .. import Plugin


class AtsasProtConvertPdbToSAXS(ProtPreprocessVolumes):
    """ Protocol for converting a PDB file (true atoms or pseudoatoms) into a
    SAXS curve.

    This is actually a wrapper to the program Crysol from Atsas.
    See documentation at:
       http://www.embl-hamburg.de/biosaxs/manuals/crysol.html
    """
    _label = 'convert PDB to SAXS curve'

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputStructure', PointerParam,
                      pointerClass='AtomStruct',
                      label="Input structure", important=True)
        form.addParam('numberOfSamples', IntParam, default=101,
                      expertLevel=LEVEL_ADVANCED,
                      label='Number of data points',
                      help='Number of calculated data points; '
                           'default: 101, maximum = 10001.')
        form.addParam('maximumFrequency', FloatParam, default=0.5,
                      expertLevel=LEVEL_ADVANCED,
                      label='Maximum scattering angle (freq.)',
                      help='Maximum scattering angle in inverse angstroms, '
                           'either for calculating the theoretical curve up '
                           'to sm or for fitting till sm; default: 0.5Å-1, '
                           'maximum: 2.0Å-1')
        form.addParam('numberOfHarmonics', IntParam, default=20,
                      expertLevel=LEVEL_ADVANCED,
                      label='Number of harmonics',
                      help='Maximum order of harmonics; default: 20, minimum: '
                           '1, maximum: 100. This defines the resolution of '
                           'the calculated curve. The default value should be '
                           'sufficient in most of the cases. For large or '
                           'extended particles higher orders could improve '
                           'the results, at the cost of an increased run time. '
                           'This value must be increased whenever the maximum '
                           'scattering angle is increased (smax).')
        form.addParam('experimentalSAXS', FileParam, filter="*.dat",
                      default='',
                      label='Experimental SAXS curve (optional)',
                      help="This parameter is optional. If provided the "
                           "simulated SAXS curve will be compared to the "
                           "experimental one.")
        form.addParam('otherCrysol', StringParam, default='',
                      expertLevel=LEVEL_ADVANCED,
                      label='Other parameters for Crysol',
                      help='See http://www.embl-hamburg.de/biosaxs/manuals/crysol.html')

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.crysolWrapper)

    # --------------------------- STEPS functions -----------------------------
    def crysolWrapper(self):
        inputFn = self.inputStructure.get().getFileName()
        ext = pwutils.getExt(inputFn)
        pwutils.createLink(os.path.abspath(inputFn),
                           self._getTmpPath(f"pseudoatoms{ext}"))

        args = [
            f"--lm {self.numberOfHarmonics}",
            f"--smax {self.maximumFrequency}",
            f"--ns {self.numberOfSamples}",
            f"{self.otherCrysol}",
            f"../tmp/pseudoatoms{ext}"
            ]

        if not self.experimentalSAXS.empty():
            pwutils.createLink(os.path.abspath(self.experimentalSAXS.get()),
                               self._getTmpPath("experimental_SAXS_curve.dat"))
            args.append(f"../tmp/experimental_SAXS_curve.dat")

        self.runJob(Plugin.getProgram('crysol'),
                    " ".join(args), cwd=self._getExtraPath())

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = []
        summary.append('Number of data points: %d' % self.numberOfSamples)
        summary.append('Maximum frequency: %f' % self.maximumFrequency)
        summary.append('Number of harmonics: %d' % self.numberOfHarmonics)

        if not self.experimentalSAXS.empty():
            summary.append(
                'Experimental SAXS curve: %s' % self.experimentalSAXS)

        # Goodness of fit
        fnInt = self._getExtraPath("pseudoatoms.int")

        if os.path.exists(fnInt):
            x = numpy.loadtxt(fnInt, skiprows=1)
            diff = numpy.log(x[:, 1]) - numpy.log(x[:, 3])
            idx = numpy.isfinite(diff)
            RMS = math.sqrt(
                1.0 / numpy.sum(idx) * numpy.dot(diff[idx], diff[idx]))
            summary.append("RMS=%f" % RMS)

        return summary

    def _citations(self):
        return ['Svergun1995']
