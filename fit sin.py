#!/usr/bin/env python
from matplotlib.pylab import *
from matplotlib import gridspec
import numpy as np
import os
import sys
from astropy.io import fits

#间隔bin越大，效果可能越好?对于此，100>50

# obsid = '7830-50' 
# obsid = '7830-97-QPO' 
# freq1 = 0.000591
# obsid='0761670301pn-50-QPO'
obsid='0761670301pn'
freq1 = 0.000568

freq = float(1 / freq1)
nfreq = 25    # 数据最大分块数量
data0 = np.loadtxt(obsid + '.txt')  # 从文件加载数据

# 时间区间选择
if len(sys.argv) > 3:
    s1 = float(sys.argv[3])
    s2 = float(sys.argv[4])
else:
    s1 = 0
    s2 = data0[-1, 0] - data0[0, 0]

# 数据预处理
n1 = data0[0][0]
for i in range(len(data0)):
    data0[i][0] -= n1  # 时间数据平移到以 0 为起点

n1 = data0[0][0]
n2 = data0[1][0] - n1  # 时间间隔
if nfreq * n2 > freq:
    nfreq = int(freq / n2)

data0 = data0[~np.isnan(data0[:, 1])]  # 去除 NaN 数据

lcv = data0[data0[:, 0] > s1 + n1]  # 筛选时间区间内的数据
lcv = lcv[lcv[:, 0] < s2 + n1]

# 初始化结果数组
res = np.array([[0.0] * 4] * nfreq)
for i in range(nfreq):
    res[i, 0] = i / float(nfreq)

# 数据分配到对应相位
for i in range(0, len(lcv)):
    time2 = int((lcv[i, 0] + (lcv[1, 0] - lcv[0, 0]) / 2) % freq / freq * nfreq)
    time1 = int((lcv[i, 0] - (lcv[1, 0] - lcv[0, 0]) / 2) % freq / freq * nfreq)
    time3 = int(lcv[i, 0] % freq / freq * nfreq)
    if time2 - time1 == 0:
        res[time3, 1] += lcv[i, 1]
        res[time3, 2] += lcv[i, 2]
        res[time3, 3] += 1
    if time2 - time1 == 1 or time2 - time1 == 1 - nfreq:
        a = ((int((lcv[i, 0] + (lcv[1, 0] - lcv[0, 0]) / 2) / freq) +
              float(time2) / nfreq) * freq - lcv[i, 0] +
             (lcv[1, 0] - lcv[0, 0]) / 2) / (lcv[1, 0] - lcv[0, 0])
        res[time3, 1] += lcv[i, 1] * (1 - a)
        res[time3, 2] += lcv[i, 2] * (1 - a)
        res[time3, 3] += (1 - a)
        if time3 == 0:
            res[nfreq - 1, 1] += lcv[i, 1] * a
            res[nfreq - 1, 2] += lcv[i, 2] * a
            res[nfreq - 1, 3] += a
        else:
            res[time3 - 1, 1] += lcv[i, 1] * a
            res[time3 - 1, 2] += lcv[i, 2] * a
            res[time3 - 1, 3] += a
        if a > 1 or a < 0:
            print('??')
            break
    if time2 - time1 == 2 or time2 - time1 == 2 - nfreq:
        print('error')
        break

# 归一化结果并计算误差
for i in range(0, nfreq):
    if np.isnan(res[i, 3]) or res[i, 3] == 0:
        res[i, 1] = np.nan
        res[i, 2] = np.nan
    else:
        res[i, 1] = res[i, 1] / res[i, 3]  # 归一化信号值
        res[i, 2] = res[i, 2] / (res[i, 3] ** 0.5)  / (res[i,3]**1.5)# 修正误差，除以根号下数据点数量

# 绘图部分
res1 = np.vstack((res, res))  # 扩展结果数组，用于绘制两周期数据
for i in range(0, nfreq):
    res1[nfreq + i, 0] += 1
plt.errorbar(res1[:, 0], res1[:, 1], res1[:, 2], c='k', label=str(s1) + "-" + str(s2))

# 设置图像标题
plt.title(f'freq={freq1}   time={1/freq1:.2f}s   num={nfreq}', fontsize=12)
plt.xlim(-0.01, 2.01)

# 第一个周期数据分析
first_period_data = res[res[:, 0] < 1]
sum_y = np.sum(first_period_data[:, 1])
num_points = len(first_period_data)
average_y = sum_y / num_points

diff_square = (first_period_data[:, 1] - average_y) ** 2
sum_diff_square = (np.sum(diff_square)) ** 0.5
variance = sum_diff_square * (num_points - 1)

# 计算周期折叠的 RMS 值、振幅和标准差
rms = (np.nanmean(first_period_data[:, 1]**2) - (np.nanmean(first_period_data[:, 1])**2))**0.5
amplitude = (np.nanmax(first_period_data[:, 1]) - np.nanmin(first_period_data[:, 1])) / 2
std_dev = np.nanstd(first_period_data[:, 1])

# 输出 RMS 值、振幅和标准差
print("第一个周期内的 RMS 值为：", rms)
print("第一个周期内的振幅为：", amplitude, '变化%:', amplitude/average_y)
print("第一个周期内的平均纵坐标值为：", average_y, 'max:', np.nanmax(first_period_data[:, 1]), 'min:', np.nanmin(first_period_data[:, 1]))
print("最大最小偏离平均值幅度：", average_y - max(first_period_data[:, 1]), average_y - min(first_period_data[:, 1]))
print("第一个周期内的样本方差为：", sum_diff_square, '标准差为：', std_dev, '总和/数量-1: ', variance, sum_diff_square, num_points)

from scipy.optimize import curve_fit

def sine_function(x, A, omega, phi, C):
    """
    定义正弦函数模型：
    y = A * sin(omega * x + phi) + C
    参数:
    - A: 振幅
    - omega: 角频率
    - phi: 相位偏移
    - C: 垂直平移量
    """
    return A * np.sin(omega * x + phi) + C

# 对周期折叠结果进行正弦拟合 (基于两个周期的数据)
# 构造用于拟合的 x_data 和 y_data，扩展为两个周期
x_data = np.concatenate((res[:, 0], res[:, 0] + 1))  # 两个周期的相位
y_data = np.concatenate((res[:, 1], res[:, 1]))      # 对应的信号值（重复一份）

# 去除 NaN 数据点
valid_mask = ~np.isnan(y_data)  # 掩码：筛选出非 NaN 的数据点
x_data = x_data[valid_mask]
y_data = y_data[valid_mask]

# 初始拟合参数的估计值 [振幅, 角频率, 相位, 平移量]
initial_guess = [np.ptp(y_data) / 2, 2 * np.pi, 0, np.mean(y_data)]

# 执行非线性最小二乘拟合
popt, pcov = curve_fit(sine_function, x_data, y_data, p0=initial_guess)

# 提取拟合参数
A, omega, phi, C = popt

# 计算拟合优度 R²
residuals = y_data - sine_function(x_data, *popt)
ss_res = np.sum(residuals**2)  # 残差平方和
ss_tot = np.sum((y_data - np.mean(y_data))**2)  # 总体平方和
r_squared = 1 - (ss_res / ss_tot)

# 输出拟合结果
T = 2 * np.pi / omega  # 计算周期（以相位单位）
T_in_seconds = T / freq1  # 转换为秒
print(f"拟合结果: 振幅 A = {A:.6f}, 角频率 ω = {omega:.6f}, 相位 φ = {phi:.6f}, 平移量 C = {C:.6f}")
print(f"拟合优度 R² = {r_squared:.6f}")
print(f"拟合周期（相位单位） T = {T:.6f}")
print(f"拟合周期（时间单位） T = {T_in_seconds:.6f} s")

# 绘制两个周期的拟合曲线
x_fit = np.linspace(0, 2, 1000)  # 两个周期的相位点
y_fit = sine_function(x_fit, *popt)

# 绘图部分
plt.errorbar(res1[:, 0], res1[:, 1], res1[:, 2], c='k', label='Original Data')  # 原始误差棒图
plt.plot(x_fit, y_fit, label='Sine Fit (Two Periods)', color='red', linestyle='--')  # 两周期的拟合曲线
plt.axhline(C, color='blue', linestyle='--', label=f'Offset (C={C:.4f})')  # 水平虚线表示偏移量 C

# 更新标题和图例
# plt.title(f'freq={freq1} (T={T_in_seconds:.2f}s)   num={nfreq}', fontsize=12)
plt.xlim(-0.01, 2.01)
# # plt.legend()
# # 隐藏图例
plt.legend().set_visible(False)

# # plt.savefig(obsid + '-period.png', dpi=100)
plt.show()
