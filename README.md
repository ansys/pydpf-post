# PyDPF-Post
[![PyAnsys](https://img.shields.io/badge/Py-Ansys-ffc107.svg?logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABAAAAAQCAIAAACQkWg2AAABDklEQVQ4jWNgoDfg5mD8vE7q/3bpVyskbW0sMRUwofHD7Dh5OBkZGBgW7/3W2tZpa2tLQEOyOzeEsfumlK2tbVpaGj4N6jIs1lpsDAwMJ278sveMY2BgCA0NFRISwqkhyQ1q/Nyd3zg4OBgYGNjZ2ePi4rB5loGBhZnhxTLJ/9ulv26Q4uVk1NXV/f///////69du4Zdg78lx//t0v+3S88rFISInD59GqIH2esIJ8G9O2/XVwhjzpw5EAam1xkkBJn/bJX+v1365hxxuCAfH9+3b9/+////48cPuNehNsS7cDEzMTAwMMzb+Q2u4dOnT2vWrMHu9ZtzxP9vl/69RVpCkBlZ3N7enoDXBwEAAA+YYitOilMVAAAAAElFTkSuQmCC)](https://docs.pyansys.com/)
[![Python](https://img.shields.io/pypi/pyversions/ansys-dpf-post?logo=pypi)](https://pypi.org/project/ansys-dpf-post/)
[![pypi](https://badge.fury.io/py/ansys-dpf-post.svg?logo=python&logoColor=white)](https://pypi.org/project/ansys-dpf-post)
[![MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

PyDPF-Post is part of the Ansys Data Processing Framework (DPF), which provides reusable operators
that let you access and transform simulation data coming from different Ansys and third-party
result files.

PyDPF-Post leverages [PyDPF-Core](https://github.com/ansys/pydpf-core), a package for building more
advanced and customized workflows using Ansys DPF. After loading a simulation's result file in PyDPF-Post,
you can extract simulation metadata and results and then apply postprocessing operations on them.

PyDPF-Core is physics-agnostic and operator-oriented. It is a direct exposure of the DPF features and
objects with a few helpers. It is not restricted to postprocessing.

PyDPF-Post on the other hand provides a high-level, physics-oriented API for postprocessing. It has a
more Pythonic and user-friendly API dedicated to postprocessing, with new objects meant to provide a specific
interface for each type of physics and analysis (such as mechanics static versus transient versus fluids).
Most importantly, there is no need to manipulate operators or workflows.

## Documentation and issues

Documentation for the latest stable release of PyPDF-Post is hosted at
[PyDPF-Post documentation](https://post.docs.pyansys.com/version/stable/).

The PyDPF-Post documentation has five sections:

- [Getting started](https://post.docs.pyansys.com/version/stable/getting_started/index.html): Learn how to
  install PyDPF-Post in user mode and quickly begin using it.
- [User guide](https://post.docs.pyansys.com/version/stable/user_guide/index.html): Understand key concepts
  for using PyDPF-Post.
- [API reference](https://post.docs.pyansys.com/version/stable/api/index.html): Understand how to use
  Python to interact programmatically with PyDPF-Post.
- [Examples](https://post.docs.pyansys.com/version/stable/examples/index.html): Explore examples
  that show how to use PyDPF-Post to perform nay different types of operations.
- [Contribute](https://post.docs.pyansys.com/version/stable/contributing.html): Learn how to
  contribute to the PyDPF-Post codebase or documentation.

In the upper right corner of the documentation's title bar, there is an option
for switching from viewing the documentation for the latest stable release
to viewing the documentation for the development version or previously
released versions.

You can also [view](https://cheatsheets.docs.pyansys.com/pydpf-post_cheat_sheet.png) or
[download](https://cheatsheets.docs.pyansys.com/pydpf-post_cheat_sheet.pdf) the
PyDPF-Post cheat sheet. This one-page reference provides syntax rules and commands
for using PyDPF-Post.

On the [PyDPF-Post Issues](https://github.com/ansys/pydpf-post/issues) page,
you can create issues to report bugs and request new features. On the
[PyDPF-Post Discussions](https://github.com/ansys/pydpf-post/discussions) page or
the [Discussions](https://discuss.ansys.com/) page on the Ansys Developer portal,
you can post questions, share ideas, and get community feedback. 

To reach the project support team, email [pyansys.core@ansys.com](mailto:pyansys.core@ansys.com).
