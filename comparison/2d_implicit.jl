
# ----------------------
# -- function         --
# ----------------------
function u_cal(u,rhoc,lam_bdx,lam_bdy,nx,ny,dx,dy,dt,lbound,norm_check_num,norm_ok)                                 # メインの計算
    delta_u=zeros(nx,ny)   # 初期Δu^nを定義(Δu^n=u^{n+1}-u^{n})
    R=zeros(nx,ny)         # 右辺Rの定義
    sum_b=0                                  # ノルムチェック用b
        
    for i in 1+lbound:nx-lbound
        for j in 1+lbound:ny-lbound
            x_numerator=lam_bdx[i,j]*(u[i+1,j]-u[i,j])-lam_bdx[i-1,j]*(u[i,j]-u[i-1,j])     # 右辺Rのxの分子
            y_numerator=lam_bdy[i,j]*(u[i,j+1]-u[i,j])-lam_bdy[i,j-1]*(u[i,j]-u[i,j-1])     # 右辺Rのyの分子
            R[i,j]=x_numerator/(rhoc[i,j]*dx^2)+y_numerator/(rhoc[i,j]*dy^2)                 # 右辺Rの計算
            sum_b += (R[i,j])^2 
        end
    end                                                                # ノルム用分母sum_bの計算

    con=0       # ループ条件
    roop_num=0  # ループ回数

    while con==0
        
        for i in 1+lbound:nx-lbound
            for j in 1+lbound:ny-lbound
                keisu = 1/dt+(lam_bdx[i,j]+lam_bdx[i-1,j])/(rhoc[i,j]*dx^2)+(lam_bdy[i,j]+lam_bdy[i,j-1])/(rhoc[i,j]*dy^2)    # Δu_i,jの係数
                Rx=lam_bdx[i,j]*delta_u[i+1,j]+lam_bdx[i-1,j]*delta_u[i-1,j]
                Ry=lam_bdy[i,j]*delta_u[i,j+1]+lam_bdy[i,j-1]*delta_u[i,j-1]
                delta_u[i,j]=(Rx/(rhoc[i,j]*dx^2)+Ry/(rhoc[i,j]*dy^2)+R[i,j])/keisu # Δuの計算
            end
        end

        if round(roop_num+1) % norm_check_num == 0
            sum_b_Ax=0 # ノルムチェック用b-Ax
        
            for i in 1+lbound:nx-lbound
                for j in 1+lbound:ny-lbound
                    # ノルムチェック用b-Axの計算
                    keisu = 1/dt+(lam_bdx[i,j]+lam_bdx[i-1,j])/(rhoc[i,j]*dx^2)+(lam_bdy[i,j]+lam_bdy[i,j-1])/(rhoc[i,j]*dy^2)    # Δu_i,jの係数
                    Rx=lam_bdx[i,j]*delta_u[i+1,j]+lam_bdx[i-1,j]*delta_u[i-1,j]
                    Ry=lam_bdy[i,j]*delta_u[i,j+1]+lam_bdy[i,j-1]*delta_u[i,j-1]
                    sum_b_Ax += ((Rx/(rhoc[i,j]*dx^2)+Ry/(rhoc[i,j]*dy^2)+R[i,j])/keisu-delta_u[i,j])^2
                end
            end
            
            norm2d=(sum_b_Ax/sum_b)^0.5    # 2乗ノルムの計算
            println(norm2d)                        # ノルムのプリント

            if norm2d < norm_ok                   # ノルムが規定以下のときループ終了
                con=1
            end

        end
        roop_num += 1    # ループカウンター
    end


    # 時間の更新
    for i in 1+lbound:nx-lbound
        for j in 1+lbound:ny-lbound
            u[i,j]=u[i,j]+delta_u[i,j]
        end
    end
    return u
end 
        
function boundary_left(u,type,temp,lam,dx,dy,nx,ny,lbound)  # x=bd_leftの境界条件
    if type==0
        for j in 1+lbound:ny-1
            u[lbound,j]=temp
        end
    elseif type==1
        for j in 1+lbound:ny-1
            u[lbound,j]=u[lbound+1,j]
        end
    elseif type==2 || type==3
        for j in 1+lbound:ny-1
            u[lbound,j]=u[lbound+1,j]+temp*dx/lam
        end
    else
        println("boundary_left miss")
    end
            
    return u
end

function boundary_right(u,type,temp,lam,dx,dy,nx,ny,lbound) # x=bd_rightの境界条件
    if type==0
        for j in 1+lbound:ny-1
            u[nx,j]=temp
        end
    elseif type==1
        for j in 1+lbound:ny-1
            u[nx,j]=u[nx-lbound,j]
        end
    elseif type==2 || type==3
        for j in 1+lbound:ny-1
            u[nx,j]=u[nx-lbound,j]+temp*dx/lam
        end
    else
        println("boundary_right miss")
    end

    return u
end
    
function boundary_lower(u,type,temp,lam,dx,dy,nx,ny,lbound) # y=bd_lowerの境界条件
    if type==0
        for i in 1+lbound:nx-1
            u[i,lbound]=temp
        end
    elseif type==1
        for i in 1+lbound:nx-1
            u[i,lbound]=u[i,lbound+1]
        end
    elseif type==2 || type==3
        for i in 1+lbound:nx-1
            u[i,lbound]=u[i,lbound+1]+temp*dx/lam
        end
    end
            
    return u
end

function boundary_upper(u,type,temp,lam,dx,dy,nx,ny,lbound) # y=bd_upperの境界条件
    if type==0
        for i in 1+lbound:nx-1
            u[i,ny]=temp
        end
    elseif type==1
        for i in 1+lbound:nx-1
            u[i,ny]=u[i,ny-lbound]
        end
    elseif type==2 || type==3
        for i in 1+lbound:nx-1
            u[i,ny]=u[i,ny-lbound]+temp*dx/lam
        end
    end
            
    return u
end

function output_T(out_num,u,x,y,lam,dx,dy,nx,ny,out_file_front,out_file_back,out_ext) # 出力について
    outlist=["x[m],y[m],T[K]\n"]
    while length(out_num) < 4
        out_num="0"*out_num
    end
    
    for i in 1:nx
        for j in 1:ny
            push!(outlist,string(x[i])*","*string(y[j])*","*string(u[i,j])*"\n")
        end
    end

    open(out_file_front*out_num*out_file_back*out_ext,"w") do f
        for i in 1:length(outlist)
            write(f,outlist[i])
        end
    end
end

function main()
    # -----------------------
    # -- initial value     --
    # -----------------------
    nt=1000                    # 時間ステップ数
    dt=1.0e-2                    # 時間刻み幅
    nx0=100                  # xの全体のセル数
    ny0=100                  # yのセル数
    dx=1.0e-3                    # xの刻み幅
    dy=1.0e-3                    # yの刻み幅
    lbound=1            # 仮想セル数
    nx=nx0+2*lbound                      # xの総セル数
    ny=ny0+2*lbound                      # yの総セル数

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
    #= 
    ディリクレ条件  =[0,温度[K]]
    ノイマン条件    =[1,何か]
    熱流束          =[2,熱流束[W/m2]]
    =# 
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
    # -- setup         --
    # ----------------------

    @time begin

        u=ones(nx,ny) *u_init       # 初期温度uの代入

        u[1+lbound,1+lbound]=u_init+init_small        # 初期発散抑え用(R=0よりsum_b=0になる)

        x=Array{Float64}(undef,nx)
        y=Array{Float64}(undef,ny)
        for i in 1:nx                    # 位置xの作成
            x[i]=(i-1)*dx-dx/2
        end
        for i in 1:ny                    # 位置yの作成
            y[i]=(i-1)*dy-dy/2
        end

        # セル中心におけるrhoc
        rhoc=ones(nx,ny)*rho*c

        # 銅のlambdaを全空間に適用,物体セル（セル(2,2)）の左面をlam_bdx[1],上面をlam_bdy[1]とする
        lam_bdx=ones(nx-1,ny-1)*lam
        lam_bdy=ones(nx-1,ny-1)*lam
        # ----------------------
        # -- main loop        --
        # ----------------------        
        output_T("0",u,x,y,lam,dx,dy,nx,ny,out_file_front,out_file_back,out_ext)

        for k in 1:nt
            u=u_cal(u,rhoc,lam_bdx,lam_bdy,nx,ny,dx,dy,dt,lbound,norm_check_num,norm_ok)
            u=boundary_left(u,left_bd[1],left_bd[2],lam,dx,dy,nx,ny,lbound)
            u=boundary_right(u,right_bd[1],right_bd[2],lam,dx,dy,nx,ny,lbound)
            u=boundary_upper(u,upper_bd[1],upper_bd[2],lam,dx,dy,nx,ny,lbound)
            u=boundary_lower(u,lower_bd[1],lower_bd[2],lam,dx,dy,nx,ny,lbound)
            if round(k+1) % every_outnum == 0
                println("nt_______________________________"*string(round(k+1)))
                output_T(string(k+1),u,x,y,lam,dx,dy,nx,ny,out_file_front,out_file_back,out_ext)
            end
            
        end
    end
end

# -- main --
main()
