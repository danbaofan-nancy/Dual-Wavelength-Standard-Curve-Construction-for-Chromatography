import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter

# ============================================================
# 1. 输入已知的摩尔吸光系数 ε (单位: M⁻¹·cm⁻¹)
# ============================================================
# 细胞色素 C（氧化型）
e_C_408 = 85084   # 408 nm
e_C_444 = 12312   # 444 nm 

# 核黄素磷酸钠
e_R_408 = 10273   # 408 nm 
e_R_444 = 16122   # 444 nm

# 光程 (cm)
l = 1.0

# ============================================================
# 2. 构建系数矩阵 A 和其逆矩阵
# ============================================================
# 矩阵形式: [e_C_408, e_R_408; e_C_444, e_R_444] * [c_C; c_R] * l = [A408; A444]
coeff_matrix = np.array([[e_C_408, e_R_408],
                         [e_C_444, e_R_444]])

# 计算逆矩阵
inv_coeff_matrix = np.linalg.inv(coeff_matrix)

# ============================================================
# 3. 读取 Excel 数据
# ============================================================
# 第1列: A408, 第2列: A444
df = pd.read_excel('吸光度数据.xlsx')

# 提取吸光度数据（假设列名分别为 'A408' 和 'A444'）
# 如果你的列名是其他名称，请修改下面的列名
A408 = df['A408'].values
A444 = df['A444'].values

# ============================================================
# 4. 逐个样品计算浓度
# ============================================================
c_C_list = []   # 细胞色素 C 浓度 (mol/L)
c_R_list = []   # 核黄素磷酸钠浓度 (mol/L)

for a408, a444 in zip(A408, A444):
    # 吸光度向量
    A_vec = np.array([a408, a444])
    
    # 解方程组: [c_C; c_R] = inv(coeff_matrix) * A_vec / l
    conc = inv_coeff_matrix @ A_vec / l
    
    c_C_list.append(conc[0])
    c_R_list.append(conc[1])

# ============================================================
# 5. 输出结果
# ============================================================
# 添加结果到 DataFrame
df['c_CytC (mol/L)'] = c_C_list
df['c_FMN (mol/L)'] = c_R_list

# 单位换算 mg/mL 
MW_C = 12400   # 细胞色素 C 分子量 g/mol
MW_R = 478     # 核黄素磷酸钠分子量 g/mol

df['c_CytC (mg/mL)'] = df['c_CytC (mol/L)'] * MW_C
df['c_FMN (mg/mL)'] = df['c_FMN (mol/L)'] * MW_R

# 保存结果
df.to_excel('计算结果.xlsx', index=False)

# 打印前几行
print("===== 计算结果 =====")
print(df.head())


# ============================================================
# 6. 画洗脱曲线
# ============================================================
def elution_plot(df, output_image=None,smooth=True,window_length=7, polyorder=2):
        # 确定洗脱体积列名
    volume_col = 'V'
    volumes = df[volume_col]

    if smooth:
        # 对细胞色素 C 和核黄素磷酸钠的浓度进行平滑处理
        df['c_CytC (mg/mL)'] = savgol_filter(df['c_CytC (mg/mL)'], window_length=window_length, polyorder=polyorder)
        df['c_FMN (mg/mL)'] = savgol_filter(df['c_FMN (mg/mL)'], window_length=window_length, polyorder=polyorder)
    
    # 创建图形
    fig, ax1 = plt.subplots(figsize=(7, 4))
    
    # ========== 左侧纵坐标 ==========
    # 绘制细胞色素c
    ax1.plot(volumes[:], df['c_CytC (mg/mL)'][:], '-o', color='orange', linewidth=3, markersize=3, 
             label='Cytochrome c', alpha=0.8)
    
    ax1.set_xlabel('Elution Volume (µL)', fontsize= 12)
    ax1.set_ylabel('Concentration-cyt c (mg/mL)', fontsize=12, color='black')
    ax1.tick_params(axis='y', labelcolor='black')
   
   # ========== 右侧纵坐标（核黄素磷酸钠）==========
    ax2 = ax1.twinx()
    ax2.plot(volumes[:], df['c_FMN (mg/mL)'][:], '-o', color='#ffd600', linewidth=3, markersize=3, 
             label='FMN', alpha=0.8)
    
    ax2.set_ylabel('Concentration-FMN (mg/mL)', fontsize=12, color='black')
    ax2.tick_params(axis='y', labelcolor='black')
    
    # 设置右侧纵坐标范围：0 到 0.04，步长 0.02
    ax2.set_ylim([0, 0.04])
    ax2.set_yticks(np.arange(0, 0.041, 0.01))  # 步长 0.02

    # ========== 图例合并 ==========
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc='best', fontsize=12)
    
    plt.tight_layout()
    
    if output_image:
        plt.savefig(output_image, dpi=300, bbox_inches='tight')
        print(f"图片已保存到: {output_image}")
    
    plt.show()
elution_plot(df, output_image='洗脱曲线.png')
