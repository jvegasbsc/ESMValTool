import pytest
import os
import sys
sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.abspath(__file__)),'..'))
from c3s_511_basic import __Diagnostic_skeleton__, Basic_Diagnostic

class TestDiagnosticSkeleton:
    def setup(self):
        pass
    def test_init(self):
        assert False
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
