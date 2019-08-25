import os
import numpy
import math
import time

start=time.time()

# -----------------------
# -- initial value     --
# -----------------------
"""
D=1.0,L=1.0に
"""
nt=20                   # 時間ステップ数
dt=1.0e-2                 # 時間刻み幅
nx0=100                   # xのセル数
dx=1.0e-2                 # xの刻み幅
lbound=1                  # 仮想セル数
nx=nx0+2*lbound           # xの総セル数

D=1.0                    # 拡散係数

# 収束判定
norm_ok=1.0e-4

# 初期温度
u_init=0.0

# 出力情報
every_outnum=1
out_file_front="time_"
out_file_back="d-2"
out_ext=".csv"

# 計算準備
cal_c=dx**2/D/dt+2

# ----------------------
# -- function         --
# ----------------------

def setup():                               # 格子生成
    global u,x                           # グローバル変数の指定

    x=[0.0] * nx
    for i in range(lbound,nx-lbound):                    # 位置xの作成
    	x[i]=(i-1)*dx+dx/2
    
    u=[u_init] * nx    # 初期温度uの代入
    for i in range(lbound,nx-lbound):
        if x[i] <= 0.5:
            u[i]=2*x[i]
        else:
            u[i]=2-2*x[i]

def u_cal():                               # メインの計算
    global u
    delta_u=[0.0]*nx

    R=[0.0]*nx
    sum_b=0
    for i in range(lbound,nx-lbound):
        R[i]=u[i-1]-2*u[i]+u[i+1]
        sum_b +=abs(R[i])

    ite=0
    con=0
    while con==0:
        for i in range(lbound,nx-lbound):
            delta_u[i]=(delta_u[i-1]+delta_u[i+1]+R[i])/cal_c

        if (ite+1) % 100 ==0:
            sum_b_Ax=0
            for i in range(lbound,nx-lbound):
                sum_b_Ax += abs(delta_u[i-1]+delta_u[i+1]+R[i]-cal_c*delta_u[i])

            norm2d=sum_b_Ax/sum_b
            print(norm2d)

            if norm2d < norm_ok:
                con=1

        ite += 1

    for i in range(lbound,nx-lbound):
        u[i]=u[i]+delta_u[i]


def output_T(out_num): # 出力について
    outlist=['x[m],T[K]']
    for i in range(lbound,nx-lbound):
        outlist.append(str(x[i])+","+str(u[i]))
    outlist='\n'.join(outlist)

    with open(out_file_front+out_num+out_file_back+out_ext,'wt') as f:
        f.write(outlist)

    
# -- main --
setup()
output_T("0")

for k in range(nt):
    u_cal()

    if (k+1) % every_outnum == 0:
        print("nt_______________________________"+str(int((k+1)/every_outnum)))
        output_T(str(int((k+1)/every_outnum)))

    
# 経過時間の書き込み
end=time.time()-start
elapsed_time=str(end)+" sec"
with open("time",'wt') as f:
    f.write(elapsed_time)
    
