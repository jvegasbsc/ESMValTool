"""Unit tests for the module :mod:`esmvaltool.diag_scripts.mlr`."""
import os

import mock

from esmvaltool.diag_scripts.mlr import _load_mlr_models


@mock.patch('esmvaltool.diag_scripts.mlr.importlib', autospec=True)
@mock.patch('esmvaltool.diag_scripts.mlr.os.walk', autospec=True)
@mock.patch('esmvaltool.diag_scripts.mlr.os.path.dirname', autospec=True)
def test_load_mlr_models(mock_dirname, mock_walk, mock_importlib):
    """Test for loading mlr models."""
    root_dir = '/root/to/something'
    models_root = os.path.join(root_dir, 'models')
    models = [
        (models_root, ['dir', '__pycache__'], ['test.py', '__init__.py']),
        (os.path.join(models_root, 'root2'), ['d'], ['__init__.py', '42.py']),
        (os.path.join(models_root, 'root3'), [], []),
        (os.path.join(models_root, 'root4'), ['d2'], ['egg.py']),
    ]
    mock_dirname.return_value = root_dir
    mock_walk.return_value = models
    _load_mlr_models()
    modules = [
        'esmvaltool.diag_scripts.mlr.models.{}'.format(mod) for mod in
        ['test', '__init__', 'root2.__init__', 'root2.42', 'root4.egg']
    ]
    calls = [mock.call(module) for module in modules]
    mock_importlib.import_module.assert_has_calls(calls)