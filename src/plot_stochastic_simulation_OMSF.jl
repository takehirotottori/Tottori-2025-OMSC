using PyPlot: plt
include("data_file_input.jl")
plt.matplotlib[:rc]("text",usetex=true)

###------------------------------------------------------------------------------------------###
### plot state & optimal state estimator
### OSC:  optimal state estimator w/o memory
### OMSC: optimal state estimator w/  memory
### idx = 0: OSC and OMSC
### idx = 1: OSC
### idx = 2: OMSC
function plot_stochastic_simulation_OSC_OMSC(E,F,Q,M,idx)
    ### Set parameters
    A = -1.0
    D = 100.0
 
    ### Plot time course
    dirname = "data_stochastic_simulation_OMSF"
    t = input_data_file_vector(dirname*"/data_time.txt")
    i = 1
    x = input_data_file_matrix(dirname*"/data_A="*string(Int(round(-1000.0*A)))*"_D="*string(Int(round(1000.0*D)))*"_F="*string(Int(round(1000.0*F)))*"_E="*string(Int(round(1000.0*E)))*"_Q="*string(Int(round(1000.0*Q)))*"_M="*string(Int(round(1000.0*M)))*"_x_memory="*string(false)*"/data_x_i="*string(i)*".txt")
    z = input_data_file_matrix(dirname*"/data_A="*string(Int(round(-1000.0*A)))*"_D="*string(Int(round(1000.0*D)))*"_F="*string(Int(round(1000.0*F)))*"_E="*string(Int(round(1000.0*E)))*"_Q="*string(Int(round(1000.0*Q)))*"_M="*string(Int(round(1000.0*M)))*"_x_memory="*string(true)*"/data_x_i="*string(i)*".txt")
    plt.plot(t,x[1,:],linewidth=2.0,"k",label=raw"$x$",zorder=2)
    if idx == 0
        plt.plot(t,x[5,:],linewidth=1.0,"tab:blue",label=raw"$\hat{x}^{*}(t,y)$",zorder=1)
        plt.plot(t,z[5,:],linewidth=1.0,"tab:orange",label=raw"$\hat{x}^{*}(t,y,z)$",zorder=3)
    elseif idx == 1
        plt.plot(t,x[5,:],linewidth=1.0,"tab:blue",label=raw"$\hat{x}^{*}$",zorder=4)
    elseif idx == 2
        plt.plot(t,z[5,:],linewidth=1.0,"tab:orange",label=raw"$\hat{x}^{*}$",zorder=4)
    end
    
    ### Plot setting
    plt.xlabel(raw"time $t$",fontsize=30)
    plt.ylabel(raw"state estimator $\hat{x}^{*}$",fontsize=30)
    plt.legend(fontsize=15)
    plt.tick_params(labelsize=15)
    plt.tight_layout()

    dirname = "fig_stochastic_simulation_OMSF"
    try mkdir(dirname) catch end
    if idx == 0
        plt.savefig(dirname*"/fig_stochastic_simulation_x_OSC_OMSC.pdf",transparent=true)
    elseif idx == 1
        plt.savefig(dirname*"/fig_stochastic_simulation_x_OSC.pdf",transparent=true)
    else
        plt.savefig(dirname*"/fig_stochastic_simulation_x_OMSC.pdf",transparent=true)
    end
    plt.show()
    plt.close()
end

###------------------------------------------------------------------------------------------###
### plot state, observation, memory, memory control
### idx = 0: state
### idx = 1: observation
### idx = 2: memory
### idx = 3: memory control
function plot_stochastic_simulation_xyz(E,F,Q,M,idx)
    ### Set parameters
    A = -1.0
    D = 100.0

    ### Plot time course
    dirname = "data_stochastic_simulation_OMSF"
    t = input_data_file_vector(dirname*"/data_time.txt")
    i = 1
    z = input_data_file_matrix(dirname*"/data_A="*string(Int(round(-1000.0*A)))*"_D="*string(Int(round(1000.0*D)))*"_F="*string(Int(round(1000.0*F)))*"_E="*string(Int(round(1000.0*E)))*"_Q="*string(Int(round(1000.0*Q)))*"_M="*string(Int(round(1000.0*M)))*"_x_memory="*string(true)*"/data_x_i="*string(i)*".txt")
    if idx == 0
        plt.plot(t,z[1,:],linewidth=2.0,"k")
    elseif idx == 1
        plt.plot(t,z[3,:],linewidth=2.0,"tab:green")
    elseif idx == 2
        plt.plot(t,z[2,:],linewidth=2.0,"tab:red")
    elseif idx == 3
        plt.plot(t,z[6,:],linewidth=2.0,"tab:blue")
    end
    
    ### Plot setting
    plt.xlabel(raw"time $t$",fontsize=30)
    if idx == 0
        plt.ylabel(raw"state $x$",fontsize=30)
    elseif idx == 1
        plt.ylabel(raw"observation $y$",fontsize=30)
    elseif idx == 2
        plt.ylabel(raw"memory $z$",fontsize=30)
    elseif idx == 3
        plt.ylabel(raw"memory control $v^{*}$",fontsize=30)
    end
    plt.tick_params(labelsize=15)
    plt.tight_layout()

    dirname = "fig_stochastic_simulation_OMSF"
    try mkdir(dirname) catch end
    if idx == 0
        plt.savefig(dirname*"/fig_stochastic_simulation_x.pdf",transparent=true)
    elseif idx == 1
        plt.savefig(dirname*"/fig_stochastic_simulation_y.pdf",transparent=true)
    elseif idx == 2
        plt.savefig(dirname*"/fig_stochastic_simulation_z.pdf",transparent=true)
    elseif idx == 3
        plt.savefig(dirname*"/fig_stochastic_simulation_v.pdf",transparent=true)
    end
    plt.show()
    plt.close()
end

###------------------------------------------------------------------------------------------###
### plot objective function
function plot_stochastic_simulation_J(E,F,Q,M)
    ### Set parameters
    A = -1.0
    D = 100.0
    dt = 0.01

    ### Plot time course
    dirname = "data_stochastic_simulation_OMSF"
    t = input_data_file_vector(dirname*"/data_time.txt")

    j = 1
    x = input_data_file_matrix(dirname*"/data_A="*string(Int(round(-1000.0*A)))*"_D="*string(Int(round(1000.0*D)))*"_F="*string(Int(round(1000.0*F)))*"_E="*string(Int(round(1000.0*E)))*"_Q="*string(Int(round(1000.0*Q)))*"_M="*string(Int(round(1000.0*M)))*"_x_memory="*string(false)*"/data_x_i="*string(j)*".txt")
    z = input_data_file_matrix(dirname*"/data_A="*string(Int(round(-1000.0*A)))*"_D="*string(Int(round(1000.0*D)))*"_F="*string(Int(round(1000.0*F)))*"_E="*string(Int(round(1000.0*E)))*"_Q="*string(Int(round(1000.0*Q)))*"_M="*string(Int(round(1000.0*M)))*"_x_memory="*string(true)*"/data_x_i="*string(j)*".txt")
    x_cum = [0.0]
    z_cum = [0.0]
    for i in 2:size(x,2)
        x_cum = vcat(x_cum,x_cum[i-1]+Q*((x[5,i]-x[1,i])^2)*dt+M*x[6,i]^2*dt)
        z_cum = vcat(z_cum,z_cum[i-1]+Q*((z[5,i]-z[1,i])^2)*dt+M*z[6,i]^2*dt)
    end
    alpha_SE = 0.1
    plt.plot(t,x_cum,linewidth=1.0,"tab:blue",alpha=alpha_SE,zorder=1)
    plt.plot(t,z_cum,linewidth=1.0,"tab:orange",alpha=alpha_SE,zorder=3)
    x_MSE = x_cum
    z_MSE = z_cum

    for j in 2:100
        x = input_data_file_matrix(dirname*"/data_A="*string(Int(round(-1000.0*A)))*"_D="*string(Int(round(1000.0*D)))*"_F="*string(Int(round(1000.0*F)))*"_E="*string(Int(round(1000.0*E)))*"_Q="*string(Int(round(1000.0*Q)))*"_M="*string(Int(round(1000.0*M)))*"_x_memory="*string(false)*"/data_x_i="*string(j)*".txt")
        z = input_data_file_matrix(dirname*"/data_A="*string(Int(round(-1000.0*A)))*"_D="*string(Int(round(1000.0*D)))*"_F="*string(Int(round(1000.0*F)))*"_E="*string(Int(round(1000.0*E)))*"_Q="*string(Int(round(1000.0*Q)))*"_M="*string(Int(round(1000.0*M)))*"_x_memory="*string(true)*"/data_x_i="*string(j)*".txt")
        
        x_cum = [0.0]
        z_cum = [0.0]
        for i in 2:size(x,2)
            x_cum = vcat(x_cum,x_cum[i-1]+Q*((x[5,i]-x[1,i])^2)*dt+M*x[6,i]^2*dt)
            z_cum = vcat(z_cum,z_cum[i-1]+Q*((z[5,i]-z[1,i])^2)*dt+M*z[6,i]^2*dt)
                end
        plt.plot(t,x_cum,linewidth=1.0,"tab:blue",alpha=alpha_SE,zorder=1)
        plt.plot(t,z_cum,linewidth=1.0,"tab:orange",alpha=alpha_SE,zorder=3)
        x_MSE = x_MSE + x_cum
        z_MSE = z_MSE + z_cum    
    end
    x_MSE = x_MSE/100.0
    z_MSE = z_MSE/100.0
    plt.plot(t,x_MSE,linewidth=3.0,"tab:blue",label=raw"$\hat{x}^{*}(t,y)$",zorder=2)
    plt.plot(t,z_MSE,linewidth=3.0,"tab:orange",label=raw"$\hat{x}^{*}(t,y,z)$",zorder=4)

    ### Plot setting
    plt.xlabel(raw"time $t$",fontsize=30)
    plt.ylabel(raw"cumulative cost",fontsize=30)
    plt.legend(fontsize=15)
    plt.tick_params(labelsize=15)
    plt.tight_layout()

    dirname = "fig_stochastic_simulation_OMSF"
    try mkdir(dirname) catch end
    plt.savefig(dirname*"/fig_stochastic_simulation_J.pdf",transparent=true)
    plt.show()
    plt.close()
end

### main
function plot_stochastic_simulation_OMSF()
    Q = 10.0
    M = 1.0
    E = 500.0
    F = 1.0
    for idx in 0:2
        plot_stochastic_simulation_OSC_OMSC(E,F,Q,M,idx)
    end
    for idx in 0:3
        plot_stochastic_simulation_xyz(E,F,Q,M,idx)
    end
    plot_stochastic_simulation_J(E,F,Q,M)
end
# plot_stochastic_simulation_OMSF()
