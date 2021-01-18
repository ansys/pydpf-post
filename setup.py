"""Installation file for ansys-dpf-post module"""
import os
from io import open as io_open
from setuptools import setup

install_requires = ['pyvista>=0.24.0',
                    'matplotlib',
                    'ansys.dpf.core>=0.2.0']


# Get version from version info
filepath = os.path.dirname(__file__)
__version__ = None
version_file = os.path.join(filepath, 'ansys', 'dpf', 'post', '_version.py')
with io_open(version_file, mode='r') as fd:
    exec(fd.read())  # execute file from raw string


readme_file = os.path.join(filepath, 'README.md')

setup(
    name='ansys-dpf-post',
    packages=['ansys.dpf.post', 'ansys.dpf.post.examples'],
    version=__version__,
    maintainer='ANSYS',
    maintainer_email='alexander.kaszynski@ansys.com',
    description='DPF-Post Python gRPC client',
    url='https://github.com/pyansys/DPF-Post',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Information Analysis',
        'Operating System :: Microsoft :: Windows',
        'Operating System :: POSIX',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    python_requires='>=3.6.*',
    install_requires=install_requires,
)
