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
    """
    Converts an atomic structure or pseudoatomic model into a theoretical SAXS scattering curve using
    the Crysol program from the ATSAS package. The protocol allows users to simulate small-angle X-ray
    scattering profiles directly from structural models in order to evaluate how well a proposed
    structure agrees with solution scattering experiments. This approach is especially useful for
    validating structural models, comparing conformational states, studying flexible systems, and
    integrating cryo-EM or atomic modeling workflows with experimental SAXS measurements. More info:
    https://www.embl-hamburg.de/biosaxs/manuals/crysol.html

    AI Generated:

    Convert PDB to SAXS Curve (AtsasProtConvertPdbToSAXS) - User Manual
        Overview

        The Convert PDB to SAXS Curve protocol generates theoretical small-angle X-ray scattering
        profiles from atomic or pseudoatomic structural models. Its main objective is to predict the
        scattering behavior of a macromolecule in solution and optionally compare the simulated profile
        against experimental SAXS data. This type of analysis is widely used in structural biology to
        evaluate whether a structural model is consistent with experimentally observed molecular shapes
        and conformations in solution.

        In integrative structural biology workflows, SAXS calculations are particularly valuable because
        they provide low-resolution information about the global organization, flexibility, and overall
        dimensions of biological assemblies. The protocol is therefore useful not only for validating
        high-resolution models, but also for studying conformational variability, multidomain proteins,
        flexible complexes, and assemblies that may adopt multiple states in solution.

        Biological Context and Applications

        SAXS experiments measure the average scattering behavior of particles in solution, making them
        highly complementary to cryo-EM, crystallography, and computational modeling. While high-resolution
        methods often capture a single structural state, SAXS reflects the ensemble-averaged behavior of
        molecules under near-native conditions.

        This protocol is commonly used when evaluating whether a reconstructed model agrees with solution
        behavior observed experimentally. For example, a cryo-EM structure may appear rigid and compact,
        whereas the SAXS profile may indicate substantial flexibility or extension in solution. Comparing
        theoretical and experimental curves can therefore reveal biologically important conformational
        differences that are not immediately visible from static structural models alone.

        The protocol is also useful during model refinement, flexible fitting, coarse-grained modeling,
        and pseudoatomic reconstruction workflows. By generating theoretical scattering curves from
        candidate structures, users can identify which conformations best explain the experimental SAXS
        observations.

        Input Structures and General Considerations

        The protocol accepts atomic structures or pseudoatomic models as input. In practical workflows,
        these models may originate from crystallography, cryo-EM fitting, molecular dynamics simulations,
        or hybrid structural modeling approaches.

        For meaningful SAXS interpretation, users should ensure that the structural model represents a
        biologically realistic state in solution. Missing flexible regions, unresolved loops, or incorrect
        oligomeric assemblies can strongly influence the agreement between theoretical and experimental
        scattering profiles.

        In many biological systems, flexibility plays a central role in SAXS interpretation. A rigid model
        may fail to reproduce the experimental curve even if the core structure is correct. Therefore,
        discrepancies between theoretical and experimental profiles should often be interpreted as possible
        evidence of conformational heterogeneity rather than simply modeling failure.

        Scattering Range and Resolution Considerations

        The maximum scattering frequency determines the highest spatial detail represented in the calculated
        SAXS curve. Lower ranges emphasize global molecular shape and envelope information, whereas higher
        frequencies include finer structural details.

        In routine biological applications, moderate scattering ranges are usually sufficient because SAXS
        inherently provides low-resolution structural information. Increasing the maximum frequency may
        improve sensitivity to fine features, but it also increases computational demands and may amplify
        noise or overfitting when comparing to experimental data.

        The number of harmonics controls the complexity and resolution of the theoretical scattering
        calculation. Compact globular proteins generally require fewer harmonics, while elongated,
        multidomain, or flexible assemblies may benefit from higher values. However, excessively large
        values may increase runtime without producing biologically meaningful improvements.

        Number of Data Points and Sampling Density

        The number of calculated data points determines how finely the scattering curve is sampled across
        the selected angular range. Higher sampling densities provide smoother and more detailed curves,
        which can be useful when comparing against high-quality experimental measurements.

        For most biological applications, standard sampling values are sufficient to capture the overall
        behavior of the scattering profile. Extremely dense sampling is generally unnecessary unless users
        are performing detailed quantitative fitting or benchmarking analyses.

        Comparison with Experimental SAXS Curves

        One of the most important capabilities of the protocol is the optional comparison between the
        simulated SAXS profile and an experimental scattering curve. This comparison allows users to assess
        how closely the structural model reproduces the observed experimental behavior.

        In structural biology practice, agreement between theoretical and experimental SAXS curves is often
        interpreted as evidence that the model captures the global dimensions and organization of the
        molecule in solution. Poor agreement may indicate conformational flexibility, incorrect assembly
        stoichiometry, missing domains, unresolved disordered regions, or biologically relevant structural
        heterogeneity.

        The protocol can also provide quantitative estimates of agreement between theoretical and
        experimental profiles. Such values are useful for comparing alternative structural hypotheses,
        evaluating refinement strategies, or selecting among competing models.

        Practical Recommendations

        For routine use, the default calculation settings are typically sufficient for globular proteins
        and moderately sized complexes. Users studying very large assemblies, elongated particles, or
        flexible systems may consider increasing the number of harmonics and extending the scattering
        range cautiously.

        When comparing against experimental data, careful preprocessing of the experimental SAXS curve is
        essential. Background subtraction, buffer matching, and proper normalization strongly influence the
        reliability of the comparison. Users should also remember that SAXS reflects solution averages, so
        biologically relevant flexibility may require ensemble modeling approaches rather than single rigid
        structures.

        In integrative workflows combining cryo-EM and SAXS, this protocol is especially useful for
        evaluating whether high-resolution reconstructions remain consistent with the solution behavior of
        the macromolecule.

        Outputs and Interpretation

        The protocol produces theoretical SAXS scattering profiles that can be visualized and compared with
        experimental measurements. These curves provide insight into the overall molecular dimensions,
        compactness, and conformational state of the structure.

        When experimental curves are included, the agreement between theoretical and observed scattering
        can guide structural interpretation and model validation. Strong agreement generally supports the
        biological plausibility of the model in solution, whereas systematic deviations may indicate
        flexibility, conformational changes, or incomplete structural representation.

        Final Perspective

        The generation of theoretical SAXS curves from structural models represents a powerful bridge
        between computational structural biology and experimental solution scattering. By connecting atomic
        or pseudoatomic models with experimentally measurable observables, this protocol supports model
        validation, conformational analysis, and integrative interpretation of biological macromolecules in
        their native solution environment.
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
