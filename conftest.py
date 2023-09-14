"""This runs at the init of the doctest pytest session."""
import doctest
from doctest import DocTestRunner
from unittest import mock

from ansys.dpf.core.misc import module_exists
import pytest

from ansys.dpf import core

# enable matplotlib off_screen plotting to avoid test interruption

if module_exists("matplotlib"):
    import matplotlib as mpl

    mpl.use("Agg")


# enable off_screen plotting to avoid test interruption
core.settings.disable_off_screen_rendering()
core.settings.bypass_pv_opengl_osmesa_crash()


class _DPFDocTestRunner(DocTestRunner):
    def run(self, test, compileflags=None, out=None, clear_globs=True):
        try:
            return DocTestRunner.run(self, test, compileflags, out, clear_globs)
        except doctest.UnexpectedException as e:
            feature_str = "Feature not supported. Upgrade the server to"
            if feature_str in str(e.exc_info):
                pass
            else:
                raise e


@pytest.fixture(autouse=True)
def _doctest_runner_dpf():
    with mock.patch("doctest.DocTestRunner", _DPFDocTestRunner):
        yield
