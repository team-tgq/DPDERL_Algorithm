import csv

import pandas as pd
import matplotlib.pyplot as plt


# 定义一个函数来安全地计算比值
def safe_divide(a, b):
    if b == 0:
        return 1  # 如果分母为零，返回NaN
    else:
        return a / b


terrain_type_dict = {
    "data/Copernicus_DSM_COG_10_N28_00_E097_00_DEM.tif": "山区",
    "data/Copernicus_DSM_COG_10_N34_00_E114_00_DEM.tif": "平原",
    "data/Copernicus_DSM_COG_10_N41_00_E119_00_DEM.tif": "丘陵"
}
plt.rcParams['font.sans-serif'] = ['SimHei']  # 替换为您安装的中文字体名称
plt.rcParams['axes.unicode_minus'] = False  # 解决负号显示异常的问题

# 加载数据
df = pd.read_csv('output/参考线优化0127_9.csv', encoding='gbk')
# 在计算平均值之前，将高度调整为50的倍数
# df['视点高度'] = (df['视点高度'] // 30) * 30

# 然后其余代码保持不变，计算平均值并绘图...

# 计算平均值
average_values = df.groupby(['视点高度', 'DEM数据文件'])['精度增量'].mean().reset_index()

# 获取地形类型
terrains = average_values['DEM数据文件'].unique()

# 创建共用y轴的子图
fig, axes = plt.subplots(1, len(terrains), figsize=(15, 5), sharey=True)

for ax, terrain in zip(axes, terrains):
    terrain_data = average_values[average_values['DEM数据文件'] == terrain]
    ax.plot(terrain_data['视点高度'], terrain_data['精度增量'])
    ax.set_title(f'Terrain: {terrain_type_dict[terrain]}')
    ax.set_xlabel('视点高度')

axes[0].set_ylabel('精度增量', fontsize=16)
plt.savefig('output/image/精度增量_0127')
plt.show()
