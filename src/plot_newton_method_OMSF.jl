using PyPlot: plt
using LinearAlgebra
include("OSC.jl")
include("vech.jl")
include("data_file_input.jl")
using PyCall
@pyimport matplotlib.patches as patches
plt.matplotlib[:rc]("text",usetex=true)

###------------------------------------------------------------------------------------------###
function label_func(idx)
    if idx == 1
        return raw"control gain $\Pi_{xx}$"
    elseif idx == 2
        return raw"control gain $\Pi_{zx}$"
    elseif idx == 3
        return raw"control gain $\Pi_{zz}$"
    elseif idx == 4
        return raw"$\Lambda_{xx}$"
    elseif idx == 5
        return raw"$\Lambda_{zx}$"
    elseif idx == 6
        return raw"$\Lambda_{zz}$"
    elseif idx == 7
        return raw"objective function $J$"
    elseif idx == 8
        return raw"estimation error $J_{Q}$"
    elseif idx == 9
        return raw"control cost $J_{M}$"
    end  
end

###------------------------------------------------------------------------------------------###
### objective function with memory
function J_func(x,m::OSC)
    Pi = inv_vech(x[1:Int(length(x)/2)])
    La = inv_vech(x[Int(length(x)/2)+1:end])
    Si = pinv(La)
    K  = Si*transpose(m.H)*pinv(m.E+m.H*Si*transpose(m.H))
    RRR = transpose(K)*(Pi*m.B-m.S)*pinv(m.R)*transpose(Pi*m.B-m.S)*K
    QQQ = m.Q + 2.0*m.S*pinv(m.R)*transpose(Pi*m.B-m.S)*K*m.H + transpose(m.H)*RRR*m.H
    return tr(QQQ*Si + RRR*m.E)
end
### objective function without memory
function J_func_OSF(m::OSC)
    A = m.A[1,1]
    D = m.D[1,1]
    E = m.E[1,1]
    Q = m.Q[1,1]
    return Q*D*E/(D-2.0*A*E)
end

###------------------------------------------------------------------------------------------###
function plot_newton_method_OMSF_idx(idx_y)
    ### parameters
    A = -1.0
    D = 100.0
    F = 1.0
    E = 1.0
    Q = 10.0
    M = 1.0

    ### plot
    y_max = -1000.0
    y_min = +1000.0
    x_flag = true
    x_boun = 0.0
    x_list = collect(1.0:0.05:7.0)
    for E in x_list
        dirname  = "data_newton_method_OMSF_positive_clustering"
        filename = dirname*"/data_A="*string(Int(round(-1000.0*A)))*"_D="*string(Int(round(1000.0*D)))*"_F="*string(Int(round(1000.0*F)))*"_E="*string(Int(round(1000.0*E)))*"_Q="*string(Int(round(1000.0*Q)))*"_M="*string(Int(round(1000.0*M)))*".txt"
        y = input_data_file_matrix(filename)
        m   = set_OMSF(A,D,F,E,Q,M)
        m_Q = set_OMSF(A,D,F,E,Q,0.0)
        m_M = set_OMSF(A,D,F,E,0.0,M)
        J_list   = []
        J_Q_list = []
        J_M_list = []
        for i in 1:size(y,2)
            J   = J_func(y[:,i],m)
            J_Q = J_func(y[:,i],m_Q)
            J_M = J_func(y[:,i],m_M)
            push!(J_list,J)
            push!(J_Q_list,J_Q)
            push!(J_M_list,J_M)
        end
        y = vcat(y,reshape(J_list,1,length(J_list)))
        y = vcat(y,reshape(J_Q_list,1,length(J_Q_list)))
        y = vcat(y,reshape(J_M_list,1,length(J_M_list)))
        J_OSC  = J_func_OSF(m)
        for i in 1:size(y,2)
            if y[6,i] < 10.0^(-3.0)
                cn = "tab:blue"
                mk = "."
                ms = 4.0
                zo = 1
            elseif y[6,i] < 1.5
                cn = "tab:green"
                mk = "."
                ms = 3.0
                zo = 2
            else
                cn = "tab:orange"
                mk = "."
                ms = 5.0
                zo = 3
            end
            plt.plot(E,y[idx_y,i],marker=mk,markersize=ms,color=cn,zorder=zo)
            y_max = ifelse(y_max < y[idx_y,i],y[idx_y,i],y_max)
            y_min = ifelse(y_min > y[idx_y,i],y[idx_y,i],y_min)
            if y[6,i] > 10.0^(-3.0) && x_flag
                if y[7,i] <= J_OSC
                    x_boun = E - (x_list[2] - x_list[1])/2.0
                    x_flag = false
                end
            end
        end
    end
    x_mod = 0.05*(x_list[end]-x_list[1])
    y_mod = 0.05*(y_max-y_min)
    plt.fill_between([x_list[1]-x_mod,x_boun],y_min-y_mod,y_max+y_mod,fc="tab:blue",alpha=0.1,zorder=0)
    plt.fill_between([x_boun,x_list[end]+x_mod],y_min-y_mod,y_max+y_mod,fc="tab:orange",alpha=0.1,zorder=0)
    plt.xlim(x_list[1]-x_mod,x_list[end]+x_mod)
    plt.ylim(y_min-y_mod,y_max+y_mod)
    
    ### save
    # plt.title(raw"memory noise $F="*string(F)*raw"$",fontsize=30)
    plt.xlabel(raw"observation noise $E$",fontsize=30)
    plt.ylabel(label_func(idx_y),fontsize=30)
    plt.tick_params(labelsize=15)
    plt.tight_layout()
    
    dirname = "fig_newton_method_OMSF"
    try mkdir(dirname) catch end
    plt.savefig(dirname*"/fig_newton_method_OMSF_x=E_F="*string(Int(round(1000.0*log10(F))))*"_y="*string(idx_y)*".pdf",transparent=true)
    plt.show()
    plt.close()
end
function plot_newton_method_OMSF()
    for idx_y in 1:9
        plot_newton_method_OMSF_idx(idx_y)
    end
end
# plot_newton_method_OMSF()