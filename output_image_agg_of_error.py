import pandas as pd
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
# 读取Excel文件
df = pd.read_excel('output/误差点聚合0308.csv')

# 将数据从宽格式转换为长格式，这样每行都代表一个邻域和高度的组合
df_long = pd.melt(df, id_vars=['视点高度'], var_name='邻域', value_name='值')

# 对于每个视点高度和邻域组合，计算值的平均
df_avg = df_long.groupby(['视点高度', '邻域']).agg({'值': 'mean'}).reset_index()


# 为邻域创建一个数值标识符（假设邻域列包含如“1邻域（XPDERL）”的字符串）
df_avg['邻域编号'] = df_avg['邻域'].apply(lambda x: int(x.split('邻域')[0]))

# 创建三维图形
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# 绘制三维散点图，x为邻域编号，y为视点高度，z为平均值
ax.scatter(df_avg['邻域编号'], df_avg['视点高度'], df_avg['值'])

# 设置坐标轴标签
ax.set_xlabel('Neighborhood')
ax.set_ylabel('Height')
ax.set_zlabel('Average Value')

# 显示图形
plt.show()
