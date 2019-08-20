import os
import numpy
import math

# -----------------------
# -- initial value     --
# -----------------------
"""
D=1.0,L=1.0に
"""
nt=200                   # 時間ステップ数
dt=1.0e-2                 # 時間刻み幅
nx0=100                   # xのセル数
dx=1.0e-2                 # xの刻み幅
lbound=1                  # 仮想セル数
nx=nx0+2*lbound           # xの総セル数

D=1.0                    # 拡散係数

# 初期温度
u_init=0.0

# 出力情報
every_outnum=10
out_dir="1d_plot_x5mm_exact"
out_file_front="time_"
out_file_back="d-2"
out_ext=".csv"

#定数
pi=math.pi

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

def u_cal(time):                               # メインの計算
    global u

    for i in range(lbound,nx-lbound):
        iti_koume=math.sin(pi*x[i])*math.exp(-pi**2*time)
        ni_koume=-math.sin(3*pi*x[i])*math.exp(-3**2*pi**2*time)/3**2
        san_koume=math.sin(5*pi*x[i])*math.exp(-5**2*pi**2*time)/5**2
        u[i]=8/pi**2*(iti_koume+ni_koume+san_koume)

def output_T(out_num): # 出力について
    outlist=['x[m],T[K]']
    for i in range(lbound,nx-lbound):
        outlist.append(str(x[i])+","+str(u[i]))
    outlist='\n'.join(outlist)

    with open(out_dir+"/"+out_file_front+out_num+out_file_back+out_ext,'wt') as f:
        f.write(outlist)

def cre_dir():
    try:
        os.mkdir(out_dir)
    except:
        pass
    
# -- main --
setup()
cre_dir()
output_T("init_")


for k in range(nt):
    u_cal((k+1)*dt)
    
    print("nt_______________________________"+str(int(k+1)))
    output_T(str(k+1))

    
