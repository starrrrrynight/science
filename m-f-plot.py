import matplotlib.pyplot as plt
import numpy as np

# 绘制质量(QPO)频率图，对于两个QPO频率的，取较高的那个
# 部分未使用的数据已注释掉

# 数据字典：天体源的质量区间和对应的QPO频率
data = {
    'GRO J1655-40': ([6.3-0.3, 6.3+0.3], [300]),
    'XTE J1550-564': ([8.4, 10.8], [180]),
    'GRS 1915+105': ([10.6, 14.4], [113]),
    'M82 X-1': ([415-63, 415+63], [3.32]),
    'NGC 1313 X-1': ([2524, 6811], [0.45]),
    'NGC 5408 X-1': ([0.7*10**5, 1.7*10**5], [20.2*10**-3]),
    'Swift J1644+57': ([(3.16-2.9)*10**5, (3.16+35.6)*10**5], [4.8*10**-3]),
    'RE J1034+396': ([(4-1.5)*10**6, (4+3)*10**6], [2.68*10**-4]),
    # 'Sgr A*': ([(4.31-0.42)*10**6, (4.31+0.42)*10**6], [1.07*10**-3]),  
    'Mrk 766': ([6.76 * 10**6, 6.87 * 10**6], [1.55 * 10**-4]),
    '1H 0707-495': ([5.2 * 10**6, 5.6 * 10**6], [1.21 * 10**-4, 2.64 * 10**-4]),
    # 'NGC 4151': ([2.5 * 10**6, 3 * 10**7], [5.78 * 10**-4]), 
    # 'ESO 113-G010': ([0.4*10**7,1*10**7], [1.24*10**-4, 6.8*10**-5]),  
    "NGC 1365": ([5*10**6, 10**7], [2.19*10**-4]),
    "NGC 4151": ([0.25*10**7, 5.6*10**7], [5.8*10**-4])
    # "4U 1608-52": ([1.4], [620, 900]),
    # "4U 1728-34": ([1.5], [870, 1150]),
    # "Cyg X-2": ([1.70, 1.72], [650, 750]),
    # "Aql X-1": ([1.4, 2], [1000, 1200])
}

# 创建名称到编号的映射，用于指定某些数据的编号
name_to_index = {
    "4U 1608-52": 1,
    "4U 1728-34": 2,
    "Cyg X-2": 3,
    "Aql X-1": 4
}

# 创建图和坐标轴
fig, ax = plt.subplots(figsize=(8.5,6))

# 遍历数据字典
for name, (mass, frequencies) in data.items():
    if name in name_to_index:
        # 对特定源（如"4U 1608-52", "Cyg X-2"等）处理其频率区间，绘制竖线和中心点
        log_frequencies = np.log10(frequencies)  # 频率取对数
        log_freq_center = np.mean(log_frequencies)  # 计算频率区间的对数中心
        log_mass = np.log10(np.sqrt(np.prod(mass)))  # 计算质量的几何平均值并取对数

        # 在对数坐标下绘制中心点，使用金色表示
        ax.loglog([10**log_mass], [10**log_freq_center], 'o', color='gold')

        # 根据不同的源调整标签位置和颜色
        if name == "4U 1608-52":
            ax.text(10**log_mass * 1.3, 10**log_freq_center * 0.8, name, ha='right', va='top', fontsize=10, color='black')
        elif name == "4U 1728-34":
            ax.text(10**log_mass * 1.2, 10**log_freq_center * 1.2, name, ha='right', va='bottom', fontsize=10, color='black')
        elif name == "Cyg X-2":
            ax.text(10**log_mass * 1.2, 10**log_freq_center * 1.3, name, ha='left', va='top', fontsize=10, color='black')
        elif name == "Aql X-1":
            ax.text(10**log_mass * 1.1, 10**log_freq_center * 1.1, name, ha='left', va='bottom', fontsize=10, color='black')

    else:
        # 处理其他天体源，绘制相应的频率点
        for frequency in frequencies:
            mass = np.where(np.array(mass) <= 0, 1e-10, mass)  # 确保没有负质量值
            log_mass = np.log10(np.sqrt(np.prod(mass)))  # 计算质量的几何平均值并取对数
            log_frequency = np.log10(frequency)  # 频率取对数

            # 根据不同源使用不同的标记颜色
            if name == "NGC 4151":
                ax.loglog([10**log_mass], [10**log_frequency], 'rs', markersize=7)  # 红色方块
            elif name in ['NGC 5408 X-1', 'NGC1313 X-1', 'M82 X-1']:
                ax.loglog([10**log_mass], [10**log_frequency], 'go')  # 紫？色圆点
            elif name in ['GRS 1915+105', 'XTE J1550-564', 'GRO J1655-40']:
                ax.loglog([10**log_mass], [10**log_frequency], 'go')  # 红色圆点
            else:
                ax.loglog([10**log_mass], [10**log_frequency], 'go')  # 绿色圆点

            # 如果质量区间内有不确定性，绘制误差棒
            if len(mass) > 1:
                ax.errorbar(10**log_mass, 10**log_frequency,
                            xerr=[[10**log_mass - 10**np.log10(mass[0])], [10**np.log10(mass[1]) - 10**log_mass]], zorder=1,
                            ecolor='black' if name != "NGC 4151" else 'red')

            # 根据不同天体源调整标签的位置和样式
            if name == "NGC 4151":
                ax.text(10**log_mass * 0.9, 10**log_frequency * 1.25, name, ha='left', va='bottom', fontsize=10, color='red')
            elif name == "1H 0707-495":
                ax.text(10**log_mass, 10**log_frequency * 0.8, name, ha='right', va='top', fontsize=10)
            elif name in ["ESO 113-G010", 'Mrk 766']:
                ax.text(10**log_mass * 1.2, 10**log_frequency * 1.2, name, ha='left', va='top', fontsize=10)
            elif name == 'RE J1034+396':
                ax.text(10**log_mass * 0.8, 10**log_frequency * 1.05, name, ha='right', va='bottom', fontsize=10)
            elif name in ['XTE J1550-564', 'GRO J1655-40']:
                ax.text(10**log_mass * 1.2, 10**log_frequency * 0.8, name, ha='left', va='bottom', fontsize=10)
            elif name == 'GRS 1915+105':
                ax.text(10**log_mass * 1.3, 10**log_frequency * 1.1, name, ha='left', va='top', fontsize=10)
            else:
                ax.text(10**log_mass, 10**log_frequency * 1.1, name, va='bottom', fontsize=10)

# 绘制参考线并放在图层的最底层
x = np.linspace(0.1, 10**8)
y1 = 6 * 10**(2 - np.log10(x))
y2 = 2 * 10**(3 - np.log10(x))
y3 = 2.2 * 10**3 / x
y4 = 16.2 * 10**3 / x
ax.loglog(x, y1, '--', color='grey', zorder=0)
ax.loglog(x, y2, '-', color='grey', zorder=0)
ax.loglog(x, y3, ':', color='grey', zorder=0)
ax.loglog(x, y4, color='c', zorder=0)

# 设置坐标轴标签和对数刻度
ax.set_xlabel('BH Mass (M$_{\\odot}$)', fontsize=12)
ax.set_ylabel('QPO frequency (Hz)', fontsize=12)
ax.set_xscale('log')
ax.set_yscale('log')

# 设置自定义次刻度
from matplotlib.ticker import FixedLocator
minor_ticks = np.concatenate([np.arange(10**i, 10**(i+1), 10**i) for i in range(-1, 8)])
ax.xaxis.set_minor_locator(FixedLocator(minor_ticks))

# 设置刻度线方向
ax.tick_params(which='both', direction='in', top=True, right=True, bottom=True, left=True)

# 设置坐标轴范围
ax.set_xlim(left=5*10**0, right=10**8)
ax.set_ylim(bottom=5 * 10**-5, top=5 * 10**2)

plt.savefig("f-m.png", format='png',dpi=1000)

# 展示图形
# plt.show()
