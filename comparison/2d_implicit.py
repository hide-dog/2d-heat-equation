import os
import numpy
import time



# -----------------------
# -- initial value     --
# -----------------------
nt=1000                  # 時間ステップ数
dt=1.0e-2                 # 時間刻み幅
nx0=100                    # xの全体のセル数
ny0=100                   # yのセル数
dx=1.0e-3                 # xの刻み幅
dy=1.0e-3                 # yの刻み幅
lbound=1                  # 仮想セル数
nx=nx0+2*lbound           # xの総セル数
ny=ny0+2*lbound           # yの総セル数

# アルミナ
lam=32.0               # 熱伝導率lambda[W/(mK)]
rho=3.94e3             # 密度[kg/m^3]
c=0.78e3               # 比熱[J/(kgK)]

# 初期温度
u_init=293.0

# norm評価値
norm_ok=1.0e-4            # これ以下に収束したらok
norm_check_num=100        # ノルムチェックを行う反復数

# 境界条件
"""
ディリクレ条件  =[0,温度[K]]
ノイマン条件    =[1,何か]
熱流束          =[2,熱流束[W/m2]]
"""
left_bd=[0,573.0]
right_bd=[0,293.0]
upper_bd=[1,0]
lower_bd=[1,0]

# 出力情報
every_outnum=100
out_file_front="time_"
out_file_back="d-3"
out_ext=".csv"

# 初期発散抑え用
init_small=1.0e-4


# ----------------------
# -- function         --
# ----------------------

def setup():                               # 格子生成
    global u,x,y,rhoc,lam_bdx,lam_bdy           # グローバル変数の指定
   
    u=[[u_init] * ny for i in [1] * nx]    # 初期温度uの代入
    u[lbound][lbound]=u_init+init_small        # 初期発散抑え用(R=0よりsum_b=0になる)
    
    x=[0.0] * nx
    y=[0.0] * ny 
    for i in range(nx):                    # 位置xの作成
    	x[i]=i*dx-dx/2
    for i in range(ny):                    # 位置yの作成
    	y[i]=i*dy-dy/2
   
    # セル中心におけるrhoとc
    rhoc=[[rho*c] * ny for i in [1] * nx]

    # 熱伝導率lambdaを全空間に適用
    lam_bdx=[[lam] * (ny-1) for i in [1] * (nx-1)]
    lam_bdy=[[lam] * (ny-1) for i in [1] * (nx-1)]


def u_cal():                                 # メインの計算
    global u,ite
    
    delta_u=[[0.0] * ny for i in [1] * nx]   # 初期Δu^nを定義(Δu^n=u^{n+1}-u^{n})
    R=[[0.0] * ny for i in [1] * nx]         # 右辺Rの定義
    sum_b=0                                  # ノルムチェック用b
        
    for i in range(lbound,nx-lbound):
        for j in range(lbound,ny-lbound):
            x_numerator=lam_bdx[i][j]*(u[i+1][j]-u[i][j])-lam_bdx[i-1][j]*(u[i][j]-u[i-1][j])     # 右辺Rのxの分子
            y_numerator=lam_bdy[i][j]*(u[i][j+1]-u[i][j])-lam_bdy[i][j-1]*(u[i][j]-u[i][j-1])     # 右辺Rのyの分子
            R[i][j]=x_numerator/(rhoc[i][j]*dx**2)+y_numerator/(rhoc[i][j]*dy**2)                 # 右辺Rの計算
            sum_b += (R[i][j])**2                                                                 # ノルム用分母sum_bの計算

    con=0       # ループ条件
    roop_num=0  # ループ回数
    
    while con==0:
        
        for i in range(lbound,nx-lbound):
            for j in range(lbound,ny-lbound):
                keisu = 1/dt+(lam_bdx[i][j]+lam_bdx[i-1][j])/(rhoc[i][j]*dx**2)+(lam_bdy[i][j]+lam_bdy[i][j-1])/(rhoc[i][j]*dy**2)    # Δu_i,jの係数
                Rx=lam_bdx[i][j]*delta_u[i+1][j]+lam_bdx[i-1][j]*delta_u[i-1][j]
                Ry=lam_bdy[i][j]*delta_u[i][j+1]+lam_bdy[i][j-1]*delta_u[i][j-1]
                delta_u[i][j]=(Rx/(rhoc[i][j]*dx**2)+Ry/(rhoc[i][j]*dy**2)+R[i][j])/keisu # Δuの計算

        if int(roop_num+1) % norm_check_num == 0:
            sum_b_Ax=0 # ノルムチェック用b-Ax
        
            for i in range(lbound,nx-lbound):
                for j in range(lbound,ny-lbound):
                    # ノルムチェック用b-Axの計算
                    keisu = 1/dt+(lam_bdx[i][j]+lam_bdx[i-1][j])/(rhoc[i][j]*dx**2)+(lam_bdy[i][j]+lam_bdy[i][j-1])/(rhoc[i][j]*dy**2)    # Δu_i,jの係数
                    Rx=lam_bdx[i][j]*delta_u[i+1][j]+lam_bdx[i-1][j]*delta_u[i-1][j]
                    Ry=lam_bdy[i][j]*delta_u[i][j+1]+lam_bdy[i][j-1]*delta_u[i][j-1]
                    sum_b_Ax += ((Rx/(rhoc[i][j]*dx**2)+Ry/(rhoc[i][j]*dy**2)+R[i][j])/keisu-delta_u[i][j])**2
            
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
    if type==2:
        for j in range(1,ny-1):
            u[lbound-1][j]=u[lbound][j]+temp*dx/lam
            
    return u

def boundary_right(u,type,temp=None): # x=bd_rightの境界条件
    if type==0:
        for j in range(1,ny-1):
            u[nx-lbound][j]=temp
    if type==1:
        for j in range(1,ny-1):
            u[nx-lbound][j]=u[nx-lbound-1][j]
    if type==2:
        for j in range(1,ny-1):
            u[nx-lbound][j]=u[nx-lbound-1][j]+temp*dx/lam

    return u
		    
def boundary_lower(u,type,temp=None): # y=bd_lowerの境界条件
    if type==0:
        for i in range(1,nx-1):
            u[i][lbound-1]=temp
    if type==1:
        for i in range(1,nx-1):
            u[i][lbound-1]=u[i][lbound]
    if type==2:
        for i in range(1,nx-1):
            u[i][lbound-1]=u[i][lbound]+temp*dx/lam    
    return u

def boundary_upper(u,type,temp=None): # y=bd_upperの境界条件
    if type==0:
        for i in range(1,nx-1):
            u[i][ny-lbound]=temp[i]
    if type==1:
        for i in range(1,nx-1):
            u[i][ny-lbound]=u[i][ny-lbound-1]
    if type==2:
        for i in range(1,nx-1):
            u[i][ny-lbound]=u[i][ny-lbound-1]+temp*dx/lam    
    return u
	
def bound(kai):    # 境界の計算
    kai=boundary_left(kai,left_bd[0],left_bd[1])
    kai=boundary_right(kai,right_bd[0],right_bd[1])
    kai=boundary_upper(kai,upper_bd[0],upper_bd[1])
    kai=boundary_lower(kai,lower_bd[0],lower_bd[1])
    return kai


def output_T(out_num): # 出力について
    outlist=["x[m],y[m],T[K]"]
    
    for i in range(nx):
        for j in range(ny):
            outlist.append(str(x[i])+","+str(y[j])+","+str(u[i][j]))
    outlist='\n'.join(outlist)

    out_num=out_num.zfill(4)
    with open(out_file_front+out_num+out_file_back+out_ext,'wt') as f:
        f.write(outlist)

    
# -- main --
start=time.time()
setup()
output_T("0")

for k in range(nt):
    u_cal()
    if int(k+1) % every_outnum==0:
        print("nt_______________________________"+str(int(k+1)))
        output_T(str(k+1))

end=time.time()-start
print(end)

with open("time",'wt') as f:
    f.write(str(end)+" sec")
