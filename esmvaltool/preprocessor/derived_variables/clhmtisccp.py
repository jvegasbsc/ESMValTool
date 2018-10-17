"""Derivation of variable `clhmtisccp`."""


import iris
from iris import Constraint

from .derived_variable import DerivedVariable


class clhmtisccp(DerivedVariable):  # noqa
    """Derivation of variable `clhmtisccp`."""

    def get_required(self, frequency):
        """Get variable `short_name` and `field` pairs required for derivation.

        Parameters
        ----------
        frequency : str
            Frequency of the desired derived variable.

        Returns
        -------
        list of tuples
            List of tuples (`short_name`, `field`) of all variables required
            for derivation.

        """
        return [('clisccp', 'T4' + frequency)]

    def calculate(self, cubes):
        """Compute ISCCP high level medium-thickness cloud area fraction.

        Parameters
        ----------
        cubes : iris.cube.CubeList
            `CubeList` containing `clisccp` (`isccp_cloud_area_fraction`).

        Returns
        -------
        iris.cube.Cube
            `Cube` containing ISCCP high level medium-thickness cloud area
            fraction.

        """
        clisccp_cube = cubes.extract_strict(
            Constraint(name='isccp_cloud_area_fraction'))

        tau = iris.Constraint(
            atmosphere_optical_thickness_due_to_cloud=lambda t: 3.6 < t <= 23.)
        plev = iris.Constraint(air_pressure=lambda p: p <= 44000.)
        clhmtisccp_cube = clisccp_cube
        clhmtisccp_cube = clhmtisccp_cube.extract(tau & plev)
        coord_names = [
            coord.standard_name for coord in clhmtisccp_cube.coords()
            if len(coord.points) > 1
        ]
        if 'atmosphere_optical_thickness_due_to_cloud' in coord_names:
            clhmtisccp_cube = clhmtisccp_cube.collapsed(
                'atmosphere_optical_thickness_due_to_cloud', iris.analysis.SUM)
        if 'air_pressure' in coord_names:
            clhmtisccp_cube = clhmtisccp_cube.collapsed('air_pressure',
                                                        iris.analysis.SUM)

        return clhmtisccp_cube
