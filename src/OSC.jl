using LinearAlgebra
include("vech.jl")

###------------------------------------------------------------------------------------------###
### model structure
struct OSC
    ### state
    A::Array{Float64}
    B::Array{Float64}
    D::Array{Float64}
    
    ### observation
    H::Array{Float64}
    E::Array{Float64}
    
    ### cost
    Q::Array{Float64}
    S::Array{Float64}
    R::Array{Float64}
end

###------------------------------------------------------------------------------------------###
### target estimation problem
function set_OMSF(A,D,F,E,Q,M)
    ### extended state
    tilde_A = zeros(Float64,2,2)
    tilde_A[1,1] = A
    tilde_B = Matrix{Float64}(I,2,2)
    tilde_B[1,1] = 0.0
    tilde_D = zeros(Float64,2,2)
    tilde_D[1,1] = D
    tilde_D[2,2] = F
    
    ### extended observation
    tilde_H = Matrix{Float64}(I,2,2)
    tilde_E = zeros(Float64,2,2)
    tilde_E[1,1] = E

    ### extended cost
    tilde_Q = zeros(Float64,2,2)
    tilde_Q[1,1] = Q
    tilde_S = zeros(Float64,2,2)
    tilde_S[1,1] = Q
    tilde_R = zeros(Float64,2,2)
    tilde_R[1,1] = Q
    tilde_R[2,2] = M
    return OSC(tilde_A,tilde_B,tilde_D,tilde_H,tilde_E,tilde_Q,tilde_S,tilde_R)
end
### target tracking problem
function set_OMSC(D,F,E,Q,R,M)
    ### extended state
    tilde_A = zeros(Float64,2,2)
    tilde_B = Matrix{Float64}(I,2,2)
    tilde_D = zeros(Float64,2,2)
    tilde_D[1,1] = D
    tilde_D[2,2] = F
    
    ### extended observation
    tilde_H = Matrix{Float64}(I,2,2)
    tilde_E = zeros(Float64,2,2)
    tilde_E[1,1] = E

    ### extended cost
    tilde_Q = zeros(Float64,2,2)
    tilde_Q[1,1] = Q
    tilde_S = zeros(Float64,2,2)
    tilde_R = zeros(Float64,2,2)
    tilde_R[1,1] = R
    tilde_R[2,2] = M
    return OSC(tilde_A,tilde_B,tilde_D,tilde_H,tilde_E,tilde_Q,tilde_S,tilde_R)
end
### target tracking problem (trap diffusion case)
function set_OMSC_bounded(D,F,E,Q,R,M)
    ### extended state
    tilde_A = zeros(Float64,2,2)
    tilde_A[1,1] = -1.0
    tilde_B = Matrix{Float64}(I,2,2)
    tilde_D = zeros(Float64,2,2)
    tilde_D[1,1] = D
    tilde_D[2,2] = F
    
    ### extended observation
    tilde_H = Matrix{Float64}(I,2,2)
    tilde_E = zeros(Float64,2,2)
    tilde_E[1,1] = E

    ### extended cost
    tilde_Q = zeros(Float64,2,2)
    tilde_Q[1,1] = Q
    tilde_S = zeros(Float64,2,2)
    tilde_R = zeros(Float64,2,2)
    tilde_R[1,1] = R
    tilde_R[2,2] = M
    return OSC(tilde_A,tilde_B,tilde_D,tilde_H,tilde_E,tilde_Q,tilde_S,tilde_R)
end

###------------------------------------------------------------------------------------------###
### steady-state equations of Pi and La
### Pi: control gain matrix
### La: precision matrix
### Si: covariance matrix
function rate_func_OSC(x,m::OSC)
    Pi = inv_vech(x[1:Int(length(x)/2)])
    La = inv_vech(x[Int(length(x)/2)+1:end])
    
    Si  = pinv(La)
    KH  = Si*transpose(m.H)*pinv(m.E+m.H*Si*transpose(m.H))*m.H
    IKH = Matrix{Float64}(I,size(KH,1),size(KH,2)) - KH
    
    PPP = (Pi*m.B-m.S)*pinv(m.R)*transpose(Pi*m.B-m.S)
    QQQ = transpose(IKH)*PPP*IKH
    AAA = m.A - m.B*pinv(m.R)*transpose(Pi*m.B-m.S)*KH

    Pi_func = vech(m.Q + transpose(m.A)*Pi + Pi*m.A - PPP + QQQ)
    La_func = vech(transpose(AAA)*La + La*(AAA) + La*m.D*La)
    return vcat(Pi_func,La_func)
end
### steady-state equations of Pi
function rate_func_OSC_Pi(Pi_vec,La_vec,m::OSC)
    Pi = inv_vech(Pi_vec)
    La = inv_vech(La_vec)
    Si = pinv(La)
    KH  = Si*transpose(m.H)*pinv(m.E+m.H*Si*transpose(m.H))*m.H
    IKH = Matrix{Float64}(I,size(KH,1),size(KH,2)) - KH
    PPP = (Pi*m.B-m.S)*pinv(m.R)*transpose(Pi*m.B-m.S)
    QQQ = transpose(IKH)*PPP*IKH
    return vech(m.Q + transpose(m.A)*Pi + Pi*m.A - PPP + QQQ)
end
### steady-state equations of La
function rate_func_OSC_La(La_vec,Pi_vec,m::OSC)
    La = inv_vech(La_vec)    
    Pi = inv_vech(Pi_vec)
    Si = pinv(La)
    KH  = Si*transpose(m.H)*pinv(m.E+m.H*Si*transpose(m.H))*m.H
    AAA = m.A - m.B*pinv(m.R)*transpose(Pi*m.B-m.S)*KH
    return vech(- transpose(AAA)*La - La*(AAA) - La*m.D*La)
end
