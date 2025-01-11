using NLsolve
using Clustering
include("OSC.jl")
include("data_file_input.jl")
include("data_file_output.jl")

###------------------------------------------------------------------------------------------###
### nls can be applied to both scalar and vector equatons
function nls(func,params...;ini=[0.0])
    my_ftol=1e-8
    if typeof(ini) <: Number
        r = nlsolve((vout,vin)->vout[1]=func(vin[1],params...),[ini],ftol=my_ftol)
        v = r.zero[1]
    else
        r = nlsolve((vout,vin)->vout .= func(vin,params...),ini,ftol=my_ftol)
        v = r.zero
    end
    return v, r.f_converged
end

###------------------------------------------------------------------------------------------###
### generate random initial values
function generate_random_positive_semidefinite_matrix(d::Int,Smin=-1000.0,Smax=+1000.0)
    S = (Smax-Smin)*rand(d,d).+Smin
    return S*transpose(S)
end
function generate_random_x_ini(d)
    Pi_ini = vech(generate_random_positive_semidefinite_matrix(d))
    La_ini = vech(generate_random_positive_semidefinite_matrix(d))
    return vcat(Pi_ini,La_ini)
end

###------------------------------------------------------------------------------------------###
### solve steady-state equations of control gain matrix and precision matrix
function newton_method_OSC(m::OSC)
    flag_first = true
    x_sol_list = zeros(Float64,size(m.D,1)*(size(m.D,1)+1),1)

    N = 1000 
    for i in 1:N
        ### newton method
        x_ini = generate_random_x_ini(size(m.A,1))
        x_sol = nls(rate_func_OSC,m,ini=x_ini)
        
        ### success or failure
        if x_sol[2]
            ### first or second and subsequent
            if flag_first
                x_sol_list[:,1] = x_sol[1]
                flag_first = false
            else
                x_sol_list = hcat(x_sol_list,x_sol[1])
            end
        end
    end
    return x_sol_list
end

###------------------------------------------------------------------------------------------###
### remove inappropriate solutions
### 1: precision matrix should be positive semidefinite
### 2: objective function should be positive
function newton_method_positive(x_sol_list,m::OSC)
    flag_first = true
    y_sol_list = zeros(Float64,size(x_sol_list,1),1)
    for i in 1:size(x_sol_list,2)
        x_sol = x_sol_list[:,i]
        Pi = inv_vech(x_sol[1:Int(length(x_sol)/2)])
        La = inv_vech(x_sol[Int(length(x_sol)/2)+1:end])

        ### check if precision matrix is positive semidefinite
        lambda_La_min = eigen(La).values[1]
        if lambda_La_min > -1e-4 && La[1,1] > 1e-4
            ### check if objective function is positive
            Si = pinv(La)
            K  = Si*transpose(m.H)*pinv(m.E+m.H*Si*transpose(m.H))
            J  = tr(Pi*m.D + (Pi*m.B-m.S)*pinv(m.R)*transpose(Pi*m.B-m.S)*K*m.E*transpose(K))
            if J > -1e-4
                if flag_first
                    y_sol_list[:,1] = x_sol
                    flag_first = false
                else
                    y_sol_list = hcat(y_sol_list,x_sol)
                end
            end
        end
    end
    ### vector -> matrix
    if typeof(y_sol_list) == Vector{Float64}
        y_sol_list = reshape(y_sol_list,length(y_sol_list),1)
    end
    return y_sol_list
end

###------------------------------------------------------------------------------------------###
### unify duplicate solutions using clustering
function newton_method_clustering(x_sol_list)
    M = zeros(Float64,size(x_sol_list,1),1)
    for n in 1:10
        R = kmeans(x_sol_list,n)
        c = sum(R.costs)
        if c < 1e-4
            M = R.centers
            break
        end
    end
    return M
end
