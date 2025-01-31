import setuptools
import os
from pathlib import Path
from pybind11.setup_helpers import Pybind11Extension, build_ext

module_path = Path(os.path.abspath(__file__)).parent.absolute()

ver = {}
with open(module_path.joinpath('version.py')) as ver_file:
    exec(ver_file.read(), ver)

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
    name="mintimegrad",
    version=ver['__version__'],
    author="Kwang Eun Jang",
    author_email="kejang@stanford.edu",
    description="Python Wrapper of Miki Lustig's tOptGrad",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/kejang/mintimegrad",
    project_urls={
        "Bug Tracker": "https://github.com/kejang/mintimegrad/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3.10",
        "Operating System :: POSIX :: Linux",
    ],
    package_dir={"": "."},
    packages=setuptools.find_packages(where="."),
    python_requires=">=3.10",
    install_requires=[
        "numpy >= 1.17.3",
        "pybind11 >=2.6.0",
    ],
    ext_modules=ext_modules,
    cmdclass={"build_ext": build_ext},
    zip_safe=False,
)
