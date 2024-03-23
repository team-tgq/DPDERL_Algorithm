import dem_data
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

plt.rcParams['font.sans-serif'] = ['SimHei']  # 替换为您安装的中文字体名称
plt.rcParams['axes.unicode_minus'] = False  # 解决负号显示异常的问题

# 山区
dem1 = dem_data.Dem("data/Copernicus_DSM_COG_10_N28_00_E097_00_DEM.tif", 30).height
# 平原
dem2 = dem_data.Dem("data/Copernicus_DSM_COG_10_N34_00_E114_00_DEM.tif", 30).height
# 丘陵
dem3 = dem_data.Dem("data/Copernicus_DSM_COG_10_N41_00_E119_00_DEM.tif", 30).height

# 创建x, y坐标网格
x = np.linspace(0, 3599, 3600)
y = np.linspace(0, 3599, 3600)
x, y = np.meshgrid(x, y)

# 选择性降采样，这里以每30个点采样1个点为例
decimation = 100
x = x[::decimation, ::decimation]
y = y[::decimation, ::decimation]
dem1 = dem1[::decimation, ::decimation]
dem2 = dem2[::decimation, ::decimation]
dem3 = dem3[::decimation, ::decimation]

# 创建一个绘图窗口并添加3个子图
fig = plt.figure(figsize=(18, 6))

# 绘制DEM1
ax1 = fig.add_subplot(131, projection='3d')
ax1.plot_surface(x, y, dem1, cmap='Greys', edgecolor='none')
ax1.set_title('DEM-山区地形', fontsize=16)

# 绘制DEM2
ax2 = fig.add_subplot(132, projection='3d')
ax2.plot_surface(x, y, dem2, cmap='Greys', edgecolor='none')
ax2.set_title('DEM-平原地形', fontsize=16)

# 绘制DEM3
ax3 = fig.add_subplot(133, projection='3d')
ax3.plot_surface(x, y, dem3, cmap='Greys', edgecolor='none')
ax3.set_title('DEM-丘陵地形', fontsize=16)

plt.savefig('output/image/三种地形3D地形图')
# 展示图形
plt.show()
