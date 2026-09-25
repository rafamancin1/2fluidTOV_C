from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import setup
from glob import glob
import os

__version__ = "0.0.1"

# Default EOS table directory baked into the module (overridable at runtime with $TWOFLUID_EOS_DIR)
EOS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "eos_tables")

# The main interface is through Pybind11Extension.
# * You can add cxx_std=11/14/17, and then build_ext can be removed.
# * You can set include_pybind11=false to add the include directory yourself,
#   say from a submodule.
#
# Note:
#   Sort input source files if you glob sources to ensure bit-for-bit
#   reproducible builds (https://github.com/pybind/python_example/pull/53)

ext_modules = [
    Pybind11Extension(
        "twofluidTOV",
        sorted(glob("src/*.cpp")),
        # Example: passing in the version to the compiled code
        define_macros=[("VERSION_INFO", __version__),
                       ("TWOFLUID_DEFAULT_EOS_DIR", '"{}"'.format(EOS_DIR))],
        extra_link_args=['-lgsl', '-lgslcblas', '-lm']
    ),
]

setup(
    name="twofluidTOV",
    version=__version__,
    author="Rafael Mancini Santos",
    author_email="rmancinisan@gmail.com",
    license="MIT",
    description="A TOV solver using pybind11",
    long_description="",
    ext_modules=ext_modules,
    extras_require={"test": "pytest"},
    # Currently, build_ext only provides an optional "highest supported C++
    # level" feature, but in the future it may provide more features.
    cmdclass={"build_ext": build_ext},
    zip_safe=False,
    python_requires=">=3.7",
)