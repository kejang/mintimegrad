import setuptools
import os
from pathlib import Path
from pybind11.setup_helpers import Pybind11Extension, build_ext

try:
    from importlib.metadata import version
except ImportError:
    from importlib_metadata import version

module_path = Path(os.path.abspath(__file__)).parent.absolute()
package_name = "mintimegrad"

try:
    pkg_version = version(package_name)
except Exception:
    pkg_version = "0.0.8"

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

base_path = os.path.dirname(__file__)

ext_modules = [
    Pybind11Extension(
        "mintimegrad",
        ["mintimegrad/mtg.cpp", "mintimegrad/mtg_c.c", "mintimegrad/spline.c"],
    ),
]

setuptools.setup(
    name=package_name,
    version=pkg_version,
    author="Kwang Eun Jang",
    author_email="kejang@stanford.edu",
    description="Python Wrapper of Miki Lustig's tOptGrad",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/kejang/mintimegrad",
    project_urls={
        "Bug Tracker": "https://github.com/kejang/mintimegrad/issues",
    },
    package_dir={"": "."},
    packages=setuptools.find_packages(where="."),
    python_requires=">=3.10",
    ext_modules=ext_modules,
    cmdclass={"build_ext": build_ext},
    zip_safe=False,
)