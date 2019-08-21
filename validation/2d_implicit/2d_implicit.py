import os
import numpy
import time

start=time.time()

# -----------------------
# -- initial value     --
# -----------------------
nt=400                   # 時間ステップ数
dt=5.0e-4                 # 時間刻み幅
nx0=100                    # xのセル数
ny0=100                   # yのセル数
dx=1.0e-2                 # xの刻み幅
dy=1.0e-2                 # yの刻み幅
lbound=1                  # 仮想セル数
nx=nx0+2*lbound           # xの総セル数
ny=ny0+2*lbound           # yの総セル数

D=1.0    # 拡散係数

# 初期温度
u_init=0.0

# norm評価値
norm_ok=1.0e-4

# 境界条件
"""
ディリクレ条件  =[0,温度[K]]
ノイマン条件    =[1,何か]
"""
left_bd=[0,0.0]
right_bd=[0,0.0]
upper_bd=[1,1]
lower_bd=[1,1]

# 出力情報
every_outnum=20
out_file_front="time_"
out_file_back="d-2"

# ----------------------
# -- function         --
# ----------------------

def setup():                               # 格子生成
    global u,x,y                           # グローバル変数の指定

    x=[0.0] * nx
    y=[0.0] * ny
    for i in range(lbound,nx-lbound):                    # 位置xの作成
    	x[i]=(i-1)*dx+dx/2
    	y[i]=(i-1)*dy+dy/2
    
    ux=[u_init] * nx    # 初期温度uの代入
    for i in range(lbound,nx-lbound):
        if x[i] <= 0.5:
            ux[i]=2*x[i]
        else:
            ux[i]=2-2*x[i]
   
    u=(numpy.array([ux]*ny).transpose()).tolist()    # 初期温度uの代入


def u_cal():                               # メインの計算
    global u,ite
    
    delta_u=[[0.0] * ny for i in [1] * nx]   # 初期Δu^nを定義(Δu^n=u^{n+1}-u^{n})
    R=[[0.0] * ny for i in [1] * nx]         # 右辺Rの定義
    sum_b=0                                  # ノルムチェック用b
        
    for i in range(lbound,nx-lbound):
        for j in range(lbound,ny-lbound):
            R[i][j]=(u[i-1][j]-2*u[i][j]+u[i+1][j])/(dx**2)+(u[i][j-1]-2*u[i][j]+u[i][j+1])/(dy**2)   # Rの計算
            sum_b += (R[i][j])**2                                                                 # ノルム用分母sum_bの計算

    con=0       # ループ条件
    roop_num=0  # ループ回数
    
    while con==0:
        
        for i in range(lbound,nx-lbound):
            for j in range(lbound,ny-lbound):
                delta_u[i][j]=((delta_u[i][j-1]+delta_u[i][j+1])/(dy**2)+(delta_u[i-1][j]+delta_u[i+1][j])/(dx**2)+R[i][j])/(1/(D*dt)+2*(1/dy**2+1/dx**2)) # Δuの計算


        if int(roop_num+1) % 100 == 0:
            sum_b_Ax=0 # ノルムチェック用b-Ax
        
            for i in range(lbound,nx-lbound):
                for j in range(lbound,ny-lbound):
                    # ノルムチェック用b-Axの計算
                    sum_b_Ax += ((delta_u[i][j-1]+delta_u[i][j+1])/(dy**2)+(delta_u[i-1][j]+delta_u[i+1][j])/(dx**2)+R[i][j]-delta_u[i][j]*(1/(D*dt)+2*(1/dy**2+1/dx**2)))**2
            
            norm2d=numpy.sqrt(sum_b_Ax/sum_b)    # 2乗ノルムの計算
            print(norm2d)                        # ノルムのプリント

            if norm2d<norm_ok:                   # ノルムが規定以下のときループ終了
                con=1

        roop_num += 1    # ループカウンター

    # 時間の更新
    for i in range(lbound,nx-lbound):
            for j in range(lbound,ny-lbound):
                u[i][j]=u[i][j]+delta_u[i][j]

    u=bound(u) # 境界の更新
    
        
def boundary_left(u,type,temp=None):  # x=bd_leftの境界条件
    if type==0:
        for j in range(1,ny-1):
            u[lbound-1][j]=temp
    if type==1:
        for j in range(1,ny-1):
            u[lbound-1][j]=u[lbound][j]
    return u

def boundary_right(u,type,temp=None): # x=bd_rightの境界条件
    if type==0:
        for j in range(1,ny-1):
            u[nx-lbound][j]=temp
    if type==1:
        for j in range(1,ny-1):
            u[nx-lbound][j]=u[nx-lbound-1][j]

    return u
		    
def boundary_lower(u,type,temp=None): # y=bd_lowerの境界条件
    if type==0:
        for i in range(1,nx-1):
            u[i][lbound-1]=temp
    if type==1:
        for i in range(1,nx-1):
            u[i][lbound-1]=u[i][lbound]
    return u

def boundary_upper(u,type,temp=None): # y=bd_upperの境界条件
    if type==0:
        for i in range(1,nx-1):
            u[i][ny-lbound]=temp
    if type==1:
        for i in range(1,nx-1):
            u[i][ny-lbound]=u[i][ny-lbound-1]
    return u
	
def bound(kai):    # 境界の計算
    kai=boundary_left(kai,left_bd[0],left_bd[1])
    kai=boundary_right(kai,right_bd[0],right_bd[1])
    kai=boundary_upper(kai,upper_bd[0],upper_bd[1])
    kai=boundary_lower(kai,lower_bd[0],lower_bd[1])
    return kai


def output_T(out_num): # 出力について
    outlist=["x[m],T[K]"]
    for i in range(lbound,nx-lbound):
        outlist.append(str(x[i])+","+str(u[i][int(ny/2)]))
    outlist='\n'.join(outlist)

    with open(out_file_front+out_num+out_file_back+'.csv','wt') as f:
        f.write(outlist)

    
# -- main --
setup()
output_T("0")

for k in range(nt):
    u_cal()
    if int(k+1) % every_outnum==0:
        print("nt_______________________________"+str(int((k+1)/every_outnum)))
        output_T(str(int((k+1)/every_outnum)))


# 経過時間の書き込み
end=time.time()-start
elapsed_time=str(end)+" sec"
with open("time",'wt') as f:
    f.write(elapsed_time)
    
