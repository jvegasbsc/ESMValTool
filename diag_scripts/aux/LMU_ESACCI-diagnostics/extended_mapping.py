# -*- coding: utf-8 -*-
import numpy as np

installed_backends = []

# note that this import statement needs to come BEFORE basemap, as
# otherwise some cartopy libraries are not found for some strange reasons
# (at least on some machines)
try:
    import cartopy.crs as ccrs
    installed_backends.append('cartopy')
except:
    print(
         'WARNING: CARTOPY seems not to be installed and can therefore not be used as plotting backend')

try:
    from mpl_toolkits.basemap import Basemap
    installed_backends.append('basemap')
except:
    print(
        'WARNING: BASEMAP seems not to be installed and can therefore not be used as plotting backend')

try:
    from matplotlib import pyplot as plt
    installed_backends.append('imshow')
except:
    raise ValueError(
        'Fatal error: You need to have a valid matplotlib installation to be able to work with GEOVAL')


from geoval.core.mapping import MapPlotGeneric, SingleMap

"""fixes object/data return errors"""


class MapPlotGeneric(MapPlotGeneric):
    
    def _draw_basemap(self, proj_prop=None, drawparallels=True,
                      vmin_polygons=None, vmax_polygons=None, **kwargs):
        """
        """
        if proj_prop is None:
            raise ValueError(
                'No projection properties are given! Please modify or ' +\
                'choose a different backend!')

        the_map = Basemap(ax=self.pax, **proj_prop)
        xm = self.x.timmean(return_data=False)

        Z = xm
        lon = self.x.lon
        lat = self.x.lat

        X, Y = the_map(lon, lat)
        self.im = the_map.pcolormesh(X, Y, Z, **kwargs)

        self.__basemap_ancillary(the_map, drawparallels=drawparallels)

        # add polygons to map
        if self.polygons is not None:
            if False:  # individual polygons
                for p in self.polygons:
                    self._add_single_polygon_basemap(the_map, p)
            else:  # plot all polygons at once
                self._add_polygons_as_collection_basemap(
                    the_map, vmin=vmin_polygons, vmax=vmax_polygons)
                
    def _draw_cartopy(self, proj_prop=None, vmin_polygons=None, 
                      vmax_polygons=None, **kwargs):
        if proj_prop is None:
            raise ValueError(
                'No projection properties are given! ' +\
                'Please modify or choose a different backend!')

        if proj_prop['projection'] in ['robin', 'TransverseMercator',
                    'mercator', 'stereo']:
            pass
        else:
            raise ValueError('Unsupported projection type')

        if hasattr(self.x, 'data'):
            plot_data_field = True
        else:
            plot_data_field = False

        if plot_data_field:
            xm = self.x.timmean(return_object=False)
            Z = xm
            lon = self.x.lon
            lat = self.x.lat

            if np.prod(lon.shape) == 0:  # no geometry
                print 'ERROR: invalid shape for plotting!'
                return

        if proj_prop['projection'] == 'robin':
            act_ccrs = ccrs.Robinson()
        elif proj_prop['projection'] == 'stereo':
            act_ccrs = ccrs.Stereographic(central_longitude=proj_prop.pop(
                'central_longitude', 0.),
                central_latitude=proj_prop.pop('central_latitude', 0.))
        elif proj_prop['projection'] == 'TransverseMercator':
            act_ccrs = ccrs.TransverseMercator(central_longitude=proj_prop.pop(
                'central_longitude', 0.),
                central_latitude=proj_prop.pop('central_latitude', 0.))
        elif proj_prop['projection'] == 'mercator':
            if 'extent' in proj_prop.keys():
                ymin = proj_prop['extent']['ymin']
                ymax = proj_prop['extent']['ymax']
            else:
                raise ValueError('Need to specify extent!')
            act_ccrs = ccrs.Mercator(central_longitude=proj_prop.pop(
                'central_longitude', 0.), min_latitude=ymin, max_latitude=ymax)
        else:
            raise ValueError('Unsupported projection')

        self.pax = self._ax2geoax(self.pax, act_ccrs)

        # add cyclic coordinates if possible
        if plot_data_field:
            if self.x._equal_lon():
                try:
                    lon1, lat1, Z1 = self._add_cyclic_to_field(
                        self.x._get_unique_lon(), lat, Z)
                except:
                    lon1 = None
                if lon1 is not None:
                    lon = lon1
                    lat = lat1
                    Z = Z1

        # plot and ancillary plots
        if 'extent' in proj_prop.keys():
            if proj_prop['projection'] == 'mercator':
                pass
            else:
                xmin = proj_prop['extent']['xmin']
                xmax = proj_prop['extent']['xmax']
                ymin = proj_prop['extent']['ymin']
                ymax = proj_prop['extent']['ymax']
                try:
                    # , crs=act_ccrs)  
                    # problem was fixed by explicitely setting CRS
                    self.pax.set_extent([xmin, xmax, ymin, ymax])
                    # NO! the problem can not be fixed by providing the CRS
                    # explicitely! this results in strange results for the
                    # final maps!
                except:
                    print 'ERROR in set_extent. This is a known problem ' +\
                    'for cartopy geoaxes (see documentation in set_extent ' +\
                    'routine). Can not be fixed here.'
                    # try workaround
                    try:
                        # problem might be fixed by explicitely setting CRS
                        self.pax.set_extent(
                            [xmin, xmax, ymin, ymax], crs=act_ccrs)
                        # CAUTION This can result however in weird plots!!!
                        # Caused problems in the past!
                    except:
                        print 'Workaround did also not work, try to ' +\
                        'continue without setting extent!'
        else:
            self.pax.set_global()  # ensure global plot
        self.pax.coastlines()

        if plot_data_field:
            try:
                self.im = self.pax.pcolormesh(
                    lon, lat, Z, transform=ccrs.PlateCarree(), **kwargs)
            except:
                print '*** WARNING: something did not work with pcolormesh ' +\
                'plotting in mapping.py'
                self.im = None
        else:
            self.im = None

        self.pax.gridlines()

        # plot polygons
        if self.polygons is not None:
            if len(self.polygons) > 0:
                if False:  # plot all polygons individually
                    for p in self.polygons:
                        self._add_single_polygon_cartopy(p)
                else:  # all polygons as collection
                    self._add_polygons_as_collection_cartopy(
                        act_ccrs, vmin=vmin_polygons, vmax=vmax_polygons)
                    
    def _draw_imshow(self, **kwargs):
        """
        draw data using imshow command
        """
        if self.pax is None:
            raise ValueError('Fatal Error: no axis for plotting specified')
        if self.pax is not None:
            self._set_axis_invisible(self.pax, frame=True)
        if self.cax is not None:
            self._set_axis_invisible(self.cax, frame=True)
        if self.zax is not None:
            self._set_axis_invisible(self.zax, frame=True)

        self.im = self.pax.imshow(self.x.timmean(return_object=False),
                                  interpolation='nearest', **kwargs)
        
    MapPlotGeneric._draw_basemap = _draw_basemap
    MapPlotGeneric._draw_cartopy = _draw_cartopy
    MapPlotGeneric._draw_imshow = _draw_imshow

        
class SingleMap(SingleMap):
    
    def _get_statistics_str(self):
        tmp_xm = self.x.timmean(return_object=True)  # from temporal mean
        s = ''
        if self.show_statistic:
            if self.stat_type == 'mean':
                me = tmp_xm.fldmean(return_data=False)
                st = tmp_xm.fldstd(return_data=False)
                assert(len(me) == 1)
                assert(len(st) == 1)
                me = me[0]
                st = st[0]
                s = 'mean: $' + str(round(me, 2)) + \
                    ' \pm ' + str(round(st, 2)) + '$'
            elif self.stat_type == 'sum':  # area sum
                me = tmp_xm.areasum()
                assert(len(me) == 1)
                me = me[0]
                s = 'sum: $' + str(round(me, 2)) + '$'
            else:
                me = np.ma.median(tmp_xm.data)
                s = 'median: $' + str(round(me, 2)) + '$'
        return s
    
    SingleMap._get_statistics_str = _get_statistics_str
    