import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
terrain_type_dict = {
    "data/Copernicus_DSM_COG_10_N28_00_E097_00_DEM.tif": "山区",
    "data/Copernicus_DSM_COG_10_N34_00_E114_00_DEM.tif": "平原",
    "data/Copernicus_DSM_COG_10_N41_00_E119_00_DEM.tif": "丘陵"
}
plt.rcParams['font.sans-serif'] = ['SimHei']  # 替换为您安装的中文字体名称
plt.rcParams['axes.unicode_minus'] = False  # 解决负号显示异常的问题

# 读取CSV文件
df = pd.read_csv('output/精度测试_0109.csv', encoding='gbk')  # 替换为你的CSV文件名

# 分组并计算统计量
grouped = df.groupby(['视点高度', 'DEM数据文件'])
stats = grouped['SPDERL算法整体错误率'].agg(['mean', 'max', lambda x: np.sqrt(np.mean(np.square(x)))])
stats.columns = ['Mean', 'Max', 'RMS']  # 重命名列，便于理解
stats = stats.reset_index()
title = ['错误率均值', '错误率最大值', '错误率均方误差']

# 获取所有唯一的地形类型
terrain_types = stats['DEM数据文件'].unique()
# 创建图形和子图，水平排列
fig, axs = plt.subplots(1, 3, figsize=(18, 5), sharey='all', sharex='all')
# 绘制每种统计量的图
for i, measure in enumerate(['Mean', 'Max', 'RMS']):
    ax = axs[i]
    for terrain_type in terrain_types:
        # 筛选当前地形类型的数据
        terrain_data = stats[stats['DEM数据文件'] == terrain_type]
        ax.plot(terrain_data['视点高度'], terrain_data[measure], label=terrain_type_dict[terrain_type])
        plt.ylim(0, 0.3)  # 设置y轴的范围是0到0.3，这里0.3就是y轴的最大值
    ax.set_title(f'{title[i]}')
    if i == 0:  # 只在第一个子图上标记y轴
        ax.set_ylabel('SPDERL算法整体错误率')
    ax.legend(title='DEM数据文件')
fig.tight_layout(pad=4)
fig.text(0.5, 0.02, '视点高度', ha='center', fontsize=12)  # 调整文本内容和位置

plt.tight_layout(rect=(0, 0.03, 1, 1))  # 调整布局以为底部文字留出空间
plt.savefig('output/image/精度测试_0116')
plt.show()
