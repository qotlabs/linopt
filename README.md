# Linopt

Linopt is a library for linear optics calculations. It consists of C++17 core library and Python bindings to it.

## Getting started

Linopt has the following directory structure:

* `lib` -- C++ library source code,
* `examples` -- examples of C++ library using,
* `pylib` -- Python bindings to C++ library,
* `pyexamples` -- Python example scripts.

### Prerequisites

Linopt uses Eigen library for matrix operations and pybind11 for Python bindings to C++ code. We use Doxygen to generate documentation from the source code comments.

* [Eigen](http://eigen.tuxfamily.org/index.php) -- library for linear algebra.
* [pybind11](https://github.com/pybind/pybind11) -- library for exposing C++ types in Python.
* [OpenMP](https://www.openmp.org) -- API for parallel programming.
* [Doxygen](http://www.doxygen.org/) -- documentation generator from source code.

### Building
Linopt is likely to compile with any C++17 compatible compiler which is supported by its dependencies (libraries from [prerequisites](#prerequisites) section). We have tested Linopt with the following compilers:

* GCC 7.2 or newer,
* Microsoft Visual Studio 2017.

We develop Linopt under [Ubuntu](https://www.ubuntu.com/) 24.04 and provide Makefiles for building. To build C++ library you should run `make` command from the `lib` directory. It will generate static version of the library `liblinopt.a`. Type `make help` for a list of available commands.

If you want to build Python version of the library go to `pylib` directory and execute `make` in it. If the build is successful Python module `linopt.so` will appear. We also provide [Qt Creator](https://www.qt.io/qt-features-libraries-apis-tools-and-ide/#ide) project files (with extension *.creator) for convenience (however, Qt Creator uses the same Makefiles for building).

Makefile for building Python version of Linopt assumes that [pybind11 Python module](https://pypi.org/project/pybind11/) is installed on your system. The easiest way to get it is via PIP:
```
pip3 install pybind11
```

### Documentation

Documentation is available [here](https://www.qotlabs.org/linopt/doc/index.html). The most recent version of documentation can be generated from the source code using Doxygen. Simply run `make doc` in the `lib` directory.

## About

### Versioning

We use [SemVer](http://semver.org/) for versioning.

Currently this project is under initial development, so the major version is 0.

### Contributing

Contributing to Linopt is done via pull request mechanism. See [CONTRIBUTING.md](CONTRIBUTING.md) for more details.

### Authors

* Struchalin Gleb <struchalin.gleb@physics.msu.ru>
* Dyakonov Ivan <iv.dyakonov@physics.msu.ru>

### License

Copyright © 2018-2025, [Quantum Optical Technologies Laboratories](https://www.qotlabs.org/en/).

Linopt is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Linopt is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See [COPYING](COPYING) and [COPYING.LESSER](COPYING.LESSER) for more details.

