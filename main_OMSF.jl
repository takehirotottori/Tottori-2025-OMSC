include("src/do_newton_method_OMSF.jl")
include("src/do_stochastic_simulation_OMSF.jl")
include("src/plot_newton_method_OMSF.jl")
include("src/plot_stochastic_simulation_OMSF.jl")

function main()
    do_newton_method_OMSF()
    do_stochastic_simulation_OMSF()
    plot_newton_method_OMSF()
    plot_stochastic_simulation_OMSF()
end
main()
