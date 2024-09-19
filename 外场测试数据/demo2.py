import matplotlib.pyplot as plt
import numpy as np
import time

# 示例数据
x = np.linspace(0, 10, 100)
y = np.sin(x)

# 创建画布和子图
fig, ax = plt.subplots()

# 循环遍历数组x和y，每隔100ms画一个点
for i in range(len(x)):
    ax.scatter(x[i], y[i])  # 画点
    plt.pause(0.1)  # 暂停100ms

# 显示图形
plt.show()
