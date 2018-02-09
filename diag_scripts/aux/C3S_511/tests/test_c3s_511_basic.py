import pytest
import os
import sys
sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.abspath(__file__)), '..'))
from c3s_511_basic import __Diagnostic_skeleton__, Basic_Diagnostic


class TestDiagnosticSkeleton:
    def setup(self):
        self.S = __Diagnostic_skeleton__()

    def test_init(self):
        self.S.__init__()
        assert isinstance(self.S.__project_info__, dict)
        assert os.path.isdir(self.S.__plot_dir__)
        assert os.path.isdir(self.S.__work_dir__)
        assert self.S.__varname__ == "var"
        assert self.S.__output_type__ == "png"
        assert self.S.__regions__ == {"example": (10, 20, -10, -20)}
        assert self.S.__verbosity_level__ == 0
        assert self.S.__debug_info__ == "No debug info"
        assert isinstance(self.S.__config__, dict)
        assert isinstance(self.S.__basetags__, list)
        assert isinstance(self.S.__infiles__, list)
        assert isinstance(self.S.authors, list)
        assert isinstance(self.S.diagname, str)

    def test_read_data(self):
        assert False

    def test_run_diagnostic(self):
        assert False

    def test_file_check(self):
        assert False

    def test_do_overview(self):
        assert False

    def test_do_mean_var(self):
        assert False

    def test_do_trends(self):
        assert False

    def test_do_extremes(self):
        assert False

    def test_do_matureity_matrix(self):
        assert False

    def test_do_prepare_report(self):
        assert False

    def test_do_report(self):
        flist = ["diag_scripts/aux/C3S_511/example_images/albedo_QA4ECV_all_models_regionalized_smean_ts.png"]
        testdict = {"resolution" : {"dx":"1 deg", "dy":"1 deg"}, "time":{"period":{"start":"2000-01-01", "end":"2000-12-31"}}, "attribute1":"none"}

        do_report(flist, "figure list test")
        do_report(testdict, "dictionary test")

