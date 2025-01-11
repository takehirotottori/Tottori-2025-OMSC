using Random
include("OSC.jl")
include("vech.jl")
include("data_file_input.jl")
include("data_file_output.jl")

###------------------------------------------------------------------------------------------###
### Rate function
function rate_func(m::OSC,Pi_vec,Si_vec,x,y)
    Pi = inv_vech(Pi_vec)
    Si = inv_vech(Si_vec)
    K  = Si*transpose(m.H)*pinv(m.E+m.H*Si*transpose(m.H))
    return m.A*x - m.B*pinv(m.R)*transpose(Pi*m.B-m.S)*K*y
end
function optimal_control_func(m::OSC,Pi_vec,Si_vec,y)
    Pi = inv_vech(Pi_vec)
    Si = inv_vech(Si_vec)
    K  = Si*transpose(m.H)*pinv(m.E+m.H*Si*transpose(m.H))
    return -pinv(m.R)*transpose(Pi*m.B-m.S)*K*y
end

###------------------------------------------------------------------------------------------###
### Stochastic simulation
function stochastic_simulation(m::OSC,Pi_vec::Vector{Float64},Si_vec::Vector{Float64},rng,t_ini,t_end,dt)
    Si_ini = zeros(Float64,2,2)
    Si_ini[1,1] = m.D[1,1]/(-2.0*m.A[1,1])
    x = sqrt.(Si_ini)*randn(rng,Float64,2)
    y = m.H*x + sqrt.(m.E)*randn(rng,Float64,2)
    u = optimal_control_func(m,Pi_vec,Si_vec,y)
    t = t_ini
    time_course = vcat(x,y,u)
    while t < t_end-dt/2.0
        x = x + (m.A*x + m.B*u)*dt + sqrt.(m.D)*randn(rng,Float64,2)*sqrt(dt)
        y = m.H*x + sqrt.(m.E)*randn(rng,Float64,2)
        u = optimal_control_func(m,Pi_vec,Si_vec,y)
        t = t + dt
        time_course = hcat(time_course,vcat(x,y,u))
    end
    return time_course
end
function stochastic_simulation_repeat(m::OSC,Pi_vec::Vector{Float64},Si_vec::Vector{Float64},dirname,t_ini,t_end,dt)
    rng = MersenneTwister(1234)
    for i in 1:100
        x_seq = stochastic_simulation(m,Pi_vec,Si_vec,rng,t_ini,t_end,dt)
        output_data_file(x_seq,dirname*"/data_x_i="*string(i)*".txt")
    end
end
function do_stochastic_simulation_OMSF_fm(flag_memory=true)
    ### Set parameters
    A = -1.0
    D = 100.0
    F = 1.0
    E = 500.0
    Q = 10.0
    M = 1.0
    m = set_OMSF(A,D,F,E,Q,M)

    ### Set dirname and filename
    dirname   = "data_stochastic_simulation_OMSF"
    try mkdir(dirname) catch end
    dirname_x = dirname*"/data_A="*string(Int(round(-1000.0*A)))*"_D="*string(Int(round(1000.0*D)))*"_F="*string(Int(round(1000.0*F)))*"_E="*string(Int(round(1000.0*E)))*"_Q="*string(Int(round(1000.0*Q)))*"_M="*string(Int(round(1000.0*M)))*"_x_memory="*string(flag_memory)
    try mkdir(dirname_x) catch end
    
    ### Set Pi and La
    Pi_vec = zeros(Float64,3)
    La_vec = zeros(Float64,3)
    Si_vec = zeros(Float64,3)
    if flag_memory
        filename_newton_method = "data_newton_method_OMSF_positive_clustering/data_A="*string(Int(round(-1000.0*A)))*"_D="*string(Int(round(1000.0*D)))*"_F="*string(Int(round(1000.0*F)))*"_E="*string(Int(round(1000.0*E)))*"_Q="*string(Int(round(1000.0*Q)))*"_M="*string(Int(round(1000.0*M)))*".txt"
        x = input_data_file_matrix(filename_newton_method)
        if typeof(x) == Vector{Float64}
            Pi_vec[1] = -2.0*A*E*E*Q/((-2.0*A*E+D)*(-2.0*A*E+D))
            Si_vec[1] = D/(-2.0*A)
        else
            for i in 1:size(x,2)
                if x[2,i] < 0.0 && x[3,i] > 1.0
                    Pi_vec = x[1:3,i]
                    La_vec = x[4:6,i]
                end
            end
            Si_vec = vech(pinv(inv_vech(La_vec)))
        end
    else
        Pi_vec[1] = -2.0*A*E*E*Q/((-2.0*A*E+D)*(-2.0*A*E+D))
        Si_vec[1] = D/(-2.0*A)
    end

    ### Stochastic simulation
    t_ini = 0.0
    t_end = 10.0
    dt    = 0.01
    output_data_file(collect(t_ini:dt:t_end),dirname*"/data_time.txt")
    stochastic_simulation_repeat(m,Pi_vec,Si_vec,dirname_x,t_ini,t_end,dt)
end
function do_stochastic_simulation_OMSF()
    do_stochastic_simulation_OMSF_fm(true)  ### w/  memory
    do_stochastic_simulation_OMSF_fm(false) ### w/o memory
end
# do_stochastic_simulation_OMSF()
