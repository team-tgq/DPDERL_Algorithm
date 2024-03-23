import pandas as pd
import matplotlib.pyplot as plt

plt.rcParams['font.sans-serif'] = ['SimHei']  # 替换为您安装的中文字体名称
plt.rcParams['axes.unicode_minus'] = False  # 解决负号显示异常的问题

# 加载CSV文件
df = pd.read_csv('output/京都测试_h_1.csv', encoding='gbk')

# 根据列'DEM数据文件'和'视点高度'进行分组，然后计算每个分组的列'时间消耗比'的均值
grouped = df.groupby(['DEM数据文件', '视点高度'])['SPDERL算法整体错误率']

mean = grouped.mean()
max_value = grouped.max()

# 计算均方误差（MSE）
mse = grouped.apply(lambda x: ((x - x.mean()) ** 2).mean())

# 将结果合并为一个新的DataFrame
result = pd.DataFrame({
    'DEM数据文件': mean.index.get_level_values(0),
    '视点高度': mean.index.get_level_values(1),
    '错误率均值': mean.values,  # 均值
    '错误率最大值': max_value.values,  # 最大值
    '错误率均方误差': mse.values  # 均方误差
})

# 将结果保存到新的CSV文件中
result.to_csv('output/预处理_整体精度.csv', index=False, encoding='gbk')
