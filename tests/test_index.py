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
