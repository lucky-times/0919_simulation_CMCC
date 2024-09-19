
import pandas as  pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
# 获取基站的位置
bs_file = r'./7个4.9G站点.xlsx'
df = pd.read_excel(bs_file)

pos = []
row, col = df.shape
for i in range(row):
    pos.append((df.iloc[i]['纬度'], df.iloc[i]['经度']))


# 使用set去除重复项，然后转换回NumPy数组
pos = np.array(list(map(tuple, set(pos))))

file1 = r'.\7.4.2多层窄波束配置--蜂窝三扇区组网\7.4.2多层窄波束配置--蜂窝三扇区组网-1km站间距-飞行高度100米-加扰50%.csv'
df_100 = pd.read_csv(file1)
# print(df_100)
filtered_df = df_100[(df_100['Altitude'] >= 117) & (df_100['Altitude'] <= 125)]
print(filtered_df.shape)

info_ue = filtered_df[['Latitude', 'Longitude', 'SS-SINR']].values

plt.plot(info_ue[:, 0], info_ue[:, 1])




# 假设你已经有一组给定的数组 x 和 y
# 示例数据：随机生成一组数据
x = info_ue[:, 0]  # 生成 0 到 10 之间的 100 个点
y = info_ue[:, 1]  # 计算这些点的正弦值

# 初始化空的图形和轴
fig, ax = plt.subplots()
points, = ax.plot([], [], 'bo')  # 初始时没有点，'bo'表示蓝色圆点

# 初始化函数，用于设置图形的初始状态
def init():
    ax.set_xlim(31.7031, 31.7236)  # x 轴范围
    ax.set_ylim(119.8309, 119.8517)  # y 轴范围
    return points,

# 更新函数，每次更新图形时调用
def update(frame):
    index = frame  # 每次更新按顺序取下一个点
    if index < len(x):
        points.set_data(x[index], y[index])  # 设置点的位置
    return points,

# 创建动画，interval 控制每隔多少毫秒更新一次
ani = FuncAnimation(fig, update, frames=len(x), init_func=init, blit=True, interval=10)
plt.show()