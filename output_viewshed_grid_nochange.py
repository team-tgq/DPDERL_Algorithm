import numpy as np
from osgeo import gdal, osr
import shapefile
from main import analysis_by_spderl_simplified


def generate_viewable_raster_nochange(visibility_array, lon, lat, r_distance, h, start_angle, end_angle, dem_path,
                                      file_name, start_lon, start_lat):
    # 假设你已经有了视域分析结果数组
    # 1代表可视，0或2代表不可视

    # 加载DEM数据import shapefile  # PyShp库
    dem_data = gdal.Open(dem_path)
    geo_transform = dem_data.GetGeoTransform()
    geo_list = list(geo_transform)
    geo_list[0] = start_lon
    geo_list[3] = start_lat

    # 根据需要处理或保存visibility_grid
    # ...
    # 设置输出路径
    output_dir = 'output/grid/'
    output_path = output_dir + file_name

    # 创建新的栅格数据集
    driver = gdal.GetDriverByName('GTiff')
    out_raster = driver.Create(output_path, visibility_array.shape[1], visibility_array.shape[0], 1, gdal.GDT_Float32)

    # 设置地理变换和坐标系统（这里假设使用DEM的原始设置）
    out_raster.SetGeoTransform(geo_list)
    out_raster_srs = osr.SpatialReference()
    out_raster_srs.ImportFromWkt(dem_data.GetProjectionRef())
    out_raster.SetProjection(out_raster_srs.ExportToWkt())

    # 写入数据并关闭文件
    out_band = out_raster.GetRasterBand(1)
    out_band.WriteArray(visibility_array)
    out_band.FlushCache()
    out_raster = None

    # 创建一个点Shapefile
    w = shapefile.Writer(f'output/shape/{file_name}', shapeType=shapefile.POINT)
    w.field('LONGITUDE', 'F', decimal=8)
    w.field('LATITUDE', 'F', decimal=8)
    w.field('radium', 'F', decimal=8)
    w.field('start_angle', 'F', decimal=8)
    w.field('end_angle', 'F', decimal=8)
    w.field('height', 'F', decimal=8)

    # 添加观察点数据
    w.point(lon, lat)
    w.record(lon, lat, r_distance, start_angle, end_angle, h)

    # 保存Shapefile
    w.close()

    print("可视性数据已保存为GeoTIFF")
