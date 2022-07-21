import setuptools
from pybind11.setup_helpers import Pybind11Extension, build_ext

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

import os
base_path = os.path.dirname(__file__)

ext_modules = [
    Pybind11Extension(
        "mintimegrad",
        ["mintimegrad/mtg.cpp"],
    ),
]

setuptools.setup(
    name="mintimegrad",
    version="0.0.2",
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
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Operating System :: POSIX :: Linux",
    ],
    package_dir={"": "."},
    packages=setuptools.find_packages(where="."),
    python_requires=">=3.7",
    install_requires=[
        "numpy >= 1.17.3",
        "pybind11 >=2.6.0",
    ],
    ext_modules=ext_modules,
    cmdclass={"build_ext": build_ext},
    zip_safe=False,
)
