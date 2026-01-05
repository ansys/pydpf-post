# Copyright (C) 2020 - 2026 ANSYS, Inc. and/or its affiliates.
# SPDX-License-Identifier: MIT
#
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

from ansys.dpf.post.index import ResultsIndex


class TestResultsIndex:
    def test_results_index(self):
        result_index = ResultsIndex(values=["U", "V", "A"], units=None)
        ref = "ResultIndex<['U', 'V', 'A']>"
        assert repr(result_index) == ref
        result_index = ResultsIndex(values=["U", "V", "A"], units=["m", "m/s"])
        ref = "ResultIndex<['U (m)', 'V (m/s)', 'A']>"
        assert repr(result_index) == ref
        result_index = ResultsIndex(values=["U", "V", "A"], units=["m", None, "m/s^2"])
        ref = "ResultIndex<['U (m)', 'V', 'A (m/s^2)']>"
        assert repr(result_index) == ref
        result_index = ResultsIndex(values=["U", "V", "A"], units=["m", "m/s", "m/s^2"])
        ref = "ResultIndex<['U (m)', 'V (m/s)', 'A (m/s^2)']>"
        assert repr(result_index) == ref

    def test_results_index_units(self):
        result_index = ResultsIndex(values=["U", "V", "A"], units=None)
        assert result_index.units == [None, None, None]
        result_index = ResultsIndex(values=["U", "V", "A"], units=["m", "m/s"])
        assert result_index.units == ["m", "m/s", None]
        result_index = ResultsIndex(values=["U", "V", "A"], units=["m", None, "m/s^2"])
        assert result_index.units == ["m", None, "m/s^2"]
        result_index = ResultsIndex(values=["U", "V", "A"], units=["m", "m/s", "m/s^2"])
        assert result_index.units == ["m", "m/s", "m/s^2"]
