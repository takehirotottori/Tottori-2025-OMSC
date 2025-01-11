### Functions for loading data from files
### For vector data
function input_data_file_vector(filename::String)
    f = open(filename,"r")
    x::Vector{Float64} = []
    for line in eachline(f)
        push!(x,parse(Float64,line))
    end
    close(f)
    return x
end

### For matrix data
function input_data_file_matrix(filename::String)
    f = open(filename,"r")
    A = parse.(Float64,split(readline(f)))
    for line in eachline(f)
        A = hcat(A,parse.(Float64,split(line)))
    end
    close(f)
    return A
end