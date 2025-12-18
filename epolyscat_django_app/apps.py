import logging
import os
from django.apps import AppConfig

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
logger = logging.getLogger(__name__)


class Settings:
    WEBPACK_LOADER = {
        "EPOLYSCAT_DJANGO_APP": {
            "BUNDLE_DIR_NAME": "epolyscat_django_app/dashboard/dist",
            "STATS_FILE": os.path.join(
                BASE_DIR,
                "static",
                "epolyscat_django_app",
                "dashboard",
                "dist",
                "webpack-stats.json",
            ),
        }
    }


class epolyscatDjangoAppConfig(AppConfig):
    name = "epolyscat_django_app"
    label = name
    verbose_name = "Epolyscat Django App"
    url_home = "epolyscat_django_app:home"
    fa_icon_class = "fa-atom"
    settings = Settings()
    default_auto_field = "django.db.models.BigAutoField"
    SCRIPTS = os.path.join(BASE_DIR, "ePolyScat", "SCRIPTS")

    APPLICATION_SETTINGS = {
        "EPOLYSCAT_DJANGO_APP": {
            "INPUT_PAGES": {
                "Discretization": [
                    {
                        "Axis": [
                            "name",
                            "lower end",
                            "upper end",
                            "order",
                            "nElem",
                            "functions",
                            "subset",
                            "exactIntegral",
                        ]
                    },
                    {
                        "Absorption": [
                            "axis",
                            "kind",
                            "theta",
                            "upper",
                            "lower",
                            "unitary",
                        ]
                    },
                    {"IndexConstraint": ["axes", "name"]},
                    {"BasisConstraint": ["axes", "kind"]},
                ],
                "Propagation": [
                    {
                        "TimePropagation": [
                            "begin",
                            "end",
                            "print",
                            "accuracy",
                            "store",
                            "fixStep",
                            "method",
                            "cutEnergy",
                            "operatorThreshold",
                        ]
                    }
                ],
                "Spectrum": [
                    {
                        "Spectrum": [
                            "radialPoints",
                            "maxEnergy",
                            "plot",
                            "energyGrid",
                            "minEnergy",
                            "symmetry12",
                        ]
                    },
                    {"Source": ["turnOff", "Region", "end"]},
                    {"Surface": ["points"]},
                ],
                "System": [
                    {
                        "Operator": [
                            "expectationValue",
                            "hamiltonian",
                            "initial",
                            "interaction",
                            "projection",
                            "smooth",
                            "truncate",
                        ]
                    },
                    {
                        "Laser": [
                            "shape",
                            "FWHM",
                            "I(W/cm2)",
                            "lambda(nm)",
                            "phiCEO",
                            "polarAngle",
                        ]
                    },
                    {"Dipole": ["acceleration", "length", "velocity"]},
                ],
            },
            "ALL_INPUTS": [
                {"Absorption": ["axis", "kind", "lower", "theta", "unitary", "upper"]},
                {
                    "Axis": [
                        "exactIntegral",
                        "functions",
                        "lower end",
                        "nCoefficients",
                        "nElem",
                        "name",
                        "order",
                        "subset",
                        "upper end",
                    ]
                },
                {"BasisConstraint": ["axes", "kind"]},
                {"ChannelsSubregion": ["tStore"]},
                {"Columbus": ["path"]},
                {"Constant": ["name"]},
                {"DEBUG": ["sample"]},
                {"Eigensolver": ["method", "shift"]},
                {"Eigen": ["select"]},
                {"FactorMatrix": ["name", "ncol", "nrow"]},
                {"Initial": ["kind", "state"]},
                {
                    "Laser": [
                        "FWHM",
                        "I(W/cm2)",
                        "azimuthAngle",
                        "ePhoton",
                        "ePhoton@I/2",
                        "ellip",
                        "ellipAng",
                        "lambda(nm)",
                        "lambda@I/2(nm)",
                        "peak",
                        "phiCEO",
                        "polarAngle",
                        "shape",
                    ]
                },
                {
                    "Operator": [
                        "expectationValue",
                        "hamiltonian",
                        "initial",
                        "interaction",
                        "parameterTerm",
                        "projection",
                        "smooth",
                        "spectrum",
                        "truncate",
                    ]
                },
                {"Parallel": ["sort"]},
                {"PlotFunction": ["algebra", "lower", "points", "upper"]},
                {"Plot_user": ["axis", "lowerBound", "points", "upperBound", "usage"]},
                {"Polar2D": ["Mmax", "Nmax", "origin", "quad", "radius"]},
                {
                    "PolarOffCenter": [
                        "Lmax",
                        "Mmax",
                        "Nmax",
                        "origin",
                        "quad",
                        "radius",
                    ]
                },
                {"Print": ["noScreen"]},
                {"Pulse": ["check"]},
                {"Source": ["Region", "end", "turnOff"]},
                {
                    "Spectrum": [
                        "amplitude",
                        "benchmark",
                        "energyGrid",
                        "maxEnergy",
                        "minEnergy",
                        "plot",
                        "radialPoints",
                        "spectrum",
                        "symmetry12",
                    ]
                },
                {"SpectrumPlot": ["noOverwrite"]},
                {"Surface": ["points"]},
                {
                    "TimePropagation": [
                        "accuracy",
                        "begin",
                        "cutEnergy",
                        "end",
                        "fixStep",
                        "method",
                        "operatorThreshold",
                        "print",
                        "store",
                    ]
                },
                {
                    "_EXPERT_": [
                        "FDgrid",
                        "FDmatrix",
                        "flatDer",
                        "fuseOperator",
                        "halfPlaneEta",
                        "halfPlaneR",
                        "halfPlaneS",
                        "newFloor",
                        "suppressTensor",
                    ]
                },
            ],
            "MASTER_LINP": os.path.join(BASE_DIR, "data", "eployscat", "linp"),
            "FILE_PLOTABLE": {
                "spec_total": "energy-differential spectrum",
                "spec_partial": "partial-wave spectrum",
                "expec": "expectation values (defined in input)",
                "Laser": "field, vector potential etc. in all directions",
                "eig": "eigenvalues",
            },
            "FILE_SCREENOUTPUT": {
                "outf": "screen-style output",
                "outspec": "screen-style output for amplitude integration",
            },
            "FILE_INPUT": {
                "inpc": "literal copy of input file",
            },
        }
    }

    APPLICATION_SETTINGS["EPOLYSCAT_DJANGO_APP"]["FILE_VIEWABLE"] = {
        "WARNINGS": "list of all warnings with detailed explanations where available",
        **APPLICATION_SETTINGS["EPOLYSCAT_DJANGO_APP"]["FILE_INPUT"],
        **APPLICATION_SETTINGS["EPOLYSCAT_DJANGO_APP"]["FILE_SCREENOUTPUT"],
        **APPLICATION_SETTINGS["EPOLYSCAT_DJANGO_APP"]["FILE_PLOTABLE"],
        "timer": "run-time for sections of the code",
        "mon": "monitor of ongoing run",
        "log": "list of all monitor outputs",
        "linp": "tRecX-processed input, including all defaults",
    }
