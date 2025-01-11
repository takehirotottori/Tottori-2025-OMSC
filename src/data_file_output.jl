### Functions for writing data to files
### For vector data
function output_data_file(A::Matrix{Float64},filename::String)
    f = open(filename,"w")
    for j in 1:size(A,2)
        for i in 1:size(A,1)
            print(f,string(A[i,j])*"   ")
        end
        println(f,"")
    end
    close(f)
    return A
end

### For matrix data
function output_data_file(x::Vector{Float64},filename::String)
    f = open(filename,"w")
    for i in 1:length(x)
        println(f,x[i])
    end
    close(f)
    return x
end