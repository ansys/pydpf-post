"""Installation file for ansys-dpf-post module"""
import os
from io import open as io_open
from setuptools import setup

install_requires = ["ansys.dpf.core>=0.3.0", "scooby"]


# Get version from version info
filepath = os.path.dirname(__file__)
__version__ = None
version_file = os.path.join(filepath, "ansys", "dpf", "post", "_version.py")
with io_open(version_file, mode="r") as fd:
    exec(fd.read())  # execute file from raw string


readme_file = os.path.join(filepath, "README.md")

setup(
    name="ansys-dpf-post",
    packages=["ansys.dpf.post", "ansys.dpf.post.examples"],
    version=__version__,
    author='ANSYS',
    author_email='ramdane.lagha@ansys.com',
    maintainer="ANSYS",
    maintainer_email="ramdane.lagha@ansys.com",
    description="DPF-Post Python client",
    url="https://github.com/pyansys/pydpf-post",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Information Analysis",
        "Operating System :: Microsoft :: Windows",
        "Operating System :: POSIX",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
    ],
    python_requires=">=3.7.*,!=3.10.*",
    extras_require={
        "plotting": ["vtk<9.1.0", "pyvista>=0.24.0", "matplotlib"],
    },
    install_requires=install_requires,
    license='MIT',
)
