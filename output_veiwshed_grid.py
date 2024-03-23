import numpy as np
from osgeo import gdal, osr
import shapefile
from main import analysis_by_spderl_simplified


def generate_viewable_raster(visibility_array, lon, lat, r_distance, h, start_angle, end_angle, dem_path,
                             file_name, x_grid_center, y_grid_center, x_grid_observe, y_grid_observe):
    # 假设你已经有了视域分析结果数组
    # 1代表可视，0或2代表不可视

    # 加载DEM数据import shapefile  # PyShp库
    dem_data = gdal.Open(dem_path)
    geo_transform = dem_data.GetGeoTransform()

    # 创建可视性栅格
    visibility_grid = np.zeros((3600, 3600), dtype=np.float32)

    # 映射视域分析数组到DEM，生成可视性栅格
    for i in range(visibility_array.shape[0]):
        for j in range(visibility_array.shape[1]):
            # 计算对应DEM上的行列号
            dem_row = y_grid_center + i - visibility_array.shape[0] // 2
            dem_col = x_grid_center + j - visibility_array.shape[1] // 2

            # 确保不越界
            if 0 <= dem_row < dem_data.RasterYSize and 0 <= dem_col < dem_data.RasterXSize:
                # 设置可视性
                if visibility_array[i, j] == 1:
                    print(f"{dem_row},{dem_col}")
                visibility_grid[dem_row, dem_col] = visibility_array[i, j]

    # 根据需要处理或保存visibility_grid
    # ...
    # 设置输出路径
    output_dir = 'output/grid/'
    output_path = output_dir + file_name

    # 创建新的栅格数据集
    driver = gdal.GetDriverByName('GTiff')
    out_raster = driver.Create(output_path, visibility_grid.shape[1], visibility_grid.shape[0], 1, gdal.GDT_Float32)

    # 设置地理变换和坐标系统（这里假设使用DEM的原始设置）
    out_raster.SetGeoTransform(geo_transform)
    out_raster_srs = osr.SpatialReference()
    out_raster_srs.ImportFromWkt(dem_data.GetProjectionRef())
    out_raster.SetProjection(out_raster_srs.ExportToWkt())

    # # 设置RGBA颜色通道，将可视性1的区域设为绿色，0为透明
    # red = np.zeros_like(visibility_grid, dtype=np.uint8)
    # green = visibility_grid * 255  # 将1转换为255，实现绿色
    # blue = np.zeros_like(visibility_grid, dtype=np.uint8)
    # alpha = visibility_grid * 255  # 透明度通道，1为不透明，0为完全透明

    # # 写入颜色通道
    # out_raster.GetRasterBand(1).WriteArray(red)
    # out_raster.GetRasterBand(2).WriteArray(green)
    # out_raster.GetRasterBand(3).WriteArray(blue)
    # out_raster.GetRasterBand(4).WriteArray(alpha)

    # 写入数据并关闭文件
    out_band = out_raster.GetRasterBand(1)
    out_band.WriteArray(visibility_grid)
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
