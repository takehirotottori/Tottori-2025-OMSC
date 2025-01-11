include("do_newton_method.jl")

###------------------------------------------------------------------------------------------###
### obtain optimal solutions of target estimation problem
function do_newton_method_OMSF()
    ### folder
    dirname = "data_newton_method_OMSF"
    dirname_positive = dirname*"_positive"
    dirname_clustering = dirname_positive*"_clustering"
    # try mkdir(dirname) catch end
    # try mkdir(dirname_positive) catch end
    try mkdir(dirname_clustering) catch end

    ### parameters
    A = -1.0
    D = 100.0
    F = 1.0
    E = 1.0
    Q = 10.0
    M = 1.0

    ### for stochastic simulation
    E_list = [500.0]
    for E in E_list
        println("OMSF: E="*string(Int(round(1000.0*E))))
        m = set_OMSF(A,D,F,E,Q,M)
        x_sol = newton_method_OSC(m)
        y_sol = newton_method_positive(x_sol,m)
        z_sol = newton_method_clustering(y_sol)
        # output_data_file(x_sol,dirname*"/data_A="*string(Int(round(-1000.0*A)))*"_D="*string(Int(round(1000.0*D)))*"_F="*string(Int(round(1000.0*F)))*"_E="*string(Int(round(1000.0*E)))*"_Q="*string(Int(round(1000.0*Q)))*"_M="*string(Int(round(1000.0*M)))*".txt")
        # output_data_file(y_sol,dirname_positive*"/data_A="*string(Int(round(-1000.0*A)))*"_D="*string(Int(round(1000.0*D)))*"_F="*string(Int(round(1000.0*F)))*"_E="*string(Int(round(1000.0*E)))*"_Q="*string(Int(round(1000.0*Q)))*"_M="*string(Int(round(1000.0*M)))*".txt")
        output_data_file(z_sol,dirname_clustering*"/data_A="*string(Int(round(-1000.0*A)))*"_D="*string(Int(round(1000.0*D)))*"_F="*string(Int(round(1000.0*F)))*"_E="*string(Int(round(1000.0*E)))*"_Q="*string(Int(round(1000.0*Q)))*"_M="*string(Int(round(1000.0*M)))*".txt")
    end

    ### for phase transition
    E_list = collect(1.0:0.05:7.0)
    for E in E_list
        println("OMSF: E="*string(Int(round(1000.0*E))))
        m = set_OMSF(A,D,F,E,Q,M)
        x_sol = newton_method_OSC(m)
        y_sol = newton_method_positive(x_sol,m)
        z_sol = newton_method_clustering(y_sol)
        # output_data_file(x_sol,dirname*"/data_A="*string(Int(round(-1000.0*A)))*"_D="*string(Int(round(1000.0*D)))*"_F="*string(Int(round(1000.0*F)))*"_E="*string(Int(round(1000.0*E)))*"_Q="*string(Int(round(1000.0*Q)))*"_M="*string(Int(round(1000.0*M)))*".txt")
        # output_data_file(y_sol,dirname_positive*"/data_A="*string(Int(round(-1000.0*A)))*"_D="*string(Int(round(1000.0*D)))*"_F="*string(Int(round(1000.0*F)))*"_E="*string(Int(round(1000.0*E)))*"_Q="*string(Int(round(1000.0*Q)))*"_M="*string(Int(round(1000.0*M)))*".txt")
        output_data_file(z_sol,dirname_clustering*"/data_A="*string(Int(round(-1000.0*A)))*"_D="*string(Int(round(1000.0*D)))*"_F="*string(Int(round(1000.0*F)))*"_E="*string(Int(round(1000.0*E)))*"_Q="*string(Int(round(1000.0*Q)))*"_M="*string(Int(round(1000.0*M)))*".txt")
    end
end
# do_newton_method_OMSF()