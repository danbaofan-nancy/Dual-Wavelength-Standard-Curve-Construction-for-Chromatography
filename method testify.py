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

# 直接读取数据
A408 = np.array([0.792,0.769,1.171])  
A444 = np.array([0.339,0.487,0.395])  

c_C_list = []   # 细胞色素 C 浓度 (mol/L)
c_R_list = []   # 核黄素磷酸钠浓度 (mol/L)


for a408, a444 in zip(A408, A444):
    # 吸光度向量
    A_vec = np.array([a408, a444])
    
    # 解方程组: [c_C; c_R] = inv(coeff_matrix) * A_vec / l
    conc = inv_coeff_matrix @ A_vec / l
    
    c_C_list.append(conc[0])
    c_R_list.append(conc[1])



MW_C = 12400   # 细胞色素 C 分子量 g/mol
MW_R = 478     # 核黄素磷酸钠分子量 g/mol


#打印结果
for i, (c_C, c_R) in enumerate(zip(c_C_list, c_R_list)):
    print(f"样品 {i+1}: 细胞色素 C 浓度 = {c_C:.6e} M, 核黄素磷酸钠浓度 = {c_R:.6e} M")
    print(f"样品 {i+1}: 细胞色素 C 浓度 = {c_C*MW_C:.6f} mg/mL, 核黄素磷酸钠浓度 = {c_R*MW_R:.6f} mg/mL")

