

import Plots
color_falves = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728","#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"];
f = Plots.font(11)
using LaTeXStrings
Plots.pyplot()

### test plots
    Plots.plot(linspace(0,35,50), Plots.fakedata(50,5),w=3, legend = (0.8,0.3), title = "Ubuntu")
    Plots.plot(linspace(0,35,50), Plots.fakedata(50,5),w=3, legend = (0.8,0.3), tickfont = f)
    Plots.hline!([0.0], line=(1.0,:solid, :red) , label = "" )
    Plots.savefig("test.pdf")
    # import PlotlyJS
    # PlotlyJS.savefig(sp::ElectronPlot, "output_filename.pdf")

### FILE
include("init.jl")
fd  = FiniteDiff(λ2 = 0.4); sol = Solution(fd);
rr = 0.02; solve_hjb!(rr, sol, fd);

### Consumption and zero asset motion lines
Plots.plot(fd.bgrid, [sol.cons[:,1] (rr*fd.bgrid+fd.zgrid[1]) sol.cons[:,2] (rr*fd.bgrid+fd.zgrid[2]) ], w = 1.5,
        color=[color_falves[1] color_falves[1] color_falves[2] color_falves[2]], style = [:solid :dash :solid :dash],
        label=[L"$c_u(a)$" L"$\dot{a}_u = 0$" L"$c_e(a)$" L"$\dot{a}_e = 0$"], legendfont = f,
        xaxis = ("bonds",(-0.15, 0.75)), yaxis = ("consumption", (0.09,0.25)))
Plots.savefig("cons_policy_case01.pdf")

Plots.plot(fd.bgrid, sol.lsavings, w = 1.5,
            label = [L"$s_u(a)$" L"$s_e(a)$"], legendfont = f,
            xaxis = ("bonds",(-0.15, 1.25)),
            yaxis = ("savings", (-0.125,0.075)))
Plots.savefig("lsavings_policy_case01.pdf")

### GRAPH for stationary equilibrium
rr_vals = linspace(0.0, __pars.ρ-1e-3)
BOND_mkt = [steady_state_resid(rr, sol, fd) for rr in rr_vals]
rr_equil = solve_steady_state(sol, fd)
    Plots.plot(BOND_mkt, rr_vals, w = 2.5,
                label =L"$S(r)$", legendfont = f,  # L for LaTeX string
                alpha = 0.6,
                xaxis = ("savings",(-0.04,0.5)),
                yaxis =L"$r$")
    Plots.vline!([0.0], line=(2.0,:solid, :red) , label = "zero net supply", legendfont = f)
    Plots.savefig("stationary_equil.pdf")

### Preparing transition
fd_init  = FiniteDiff(λ2 = 0.4); sol_init = Solution(fd_init);
fd_end  = FiniteDiff(λ2 = 0.8); sol_end = Solution(fd_end);

rr_init = solve_steady_state(sol_init, fd_init)
rr_end  = solve_steady_state(sol_end, fd_end)

### Consumption
Plots.plot(fd.bgrid, [sol_init.cons sol_end.cons], w = 1.5,
        color=[:red :blue :red :blue], style = [:dash :dash :solid :solid], label=["c_unemp (low)" "c_emp (low)" "c_unemp (high)" "c_emp (high)"],
        xaxis = ("bonds",(-0.15, 1.5)), yaxis = ("c", (0.05,0.30)))
Plots.savefig("cons_policy.pdf")

### Savings
Plots.plot(fd.bgrid, [sol_init.lsavings sol_end.lsavings], w = 1.5,
            color=[color_falves[1] color_falves[2] color_falves[1] color_falves[2]],
            style = [:dash :dash :solid :solid],
            label=[L"$s_u(a),\ \lambda_2 = 0.4$" L"$s_e(a),\ \lambda_2 = 0.4$" L"$s_u(a),\ \lambda_2 = 0.6$" L"$s_e(a),\ \lambda_2 = 0.6$"], legendfont = f,
            # label=["sav_unemp (low)" "sav_emp (low)" "sav_unemp (high)" "sav_emp (high)"],
            xaxis = ("bonds",(-0.15, 1.5)),
            yaxis = ("savings", (-0.125,0.075)))
Plots.savefig("lsav_pol_init_x_end.pdf")

### Compute transition
BOND_mkt = zeros(401, 201);         # line is time, col is it
BOND_mkt_clear_dist = zeros(201);
R = zeros(401, 201);                # rate over it
# r_guess = squeeze(readdlm("output/r_guess.txt"),2)
sol_transition = transition_dynamics(R, BOND_mkt, BOND_mkt_clear_dist, sol_init, sol_end, rr_end, fd_end, r_guess)
sol_transition = transition_dynamics(R, BOND_mkt, BOND_mkt_clear_dist, sol_init, sol_end, rr_end, fd_end)

# **************************************************************************************
# PLOTS
#
# ======================================================================================
# writedlm("output/r_guess.txt", r_guess)
#

### Distributions
Plots.plot(fd.bgrid, [reshape(sol_init.g,1000,2) reshape(sol_end.g,1000,2)], w = 1.5,
            color=[:red :blue :red :blue], style = [:dash :dash :solid :solid], label=["dist_unemp (low)" "dist_emp (low)" "dist_unemp (high)" "dist_emp (high)"],
            xaxis = ("bonds",(-0.15, 0.50), -0.10:0.10:0.50), yaxis= ("density",(0.0,3.0) ) , legend = :right)
Plots.savefig("distribution.pdf")

Plots.plot(fd.bgrid, [reshape(sol_init.g,1000,2)[:,2]*(0.6/(0.6+0.4))^(-1) reshape(sol_end.g,1000,2)[:,2]*(0.6/(0.6+0.8))^(-1) ], w = 1.5,
            color=[:red :blue :red :blue], style = [:dash :dash :solid :solid], label=["dist_emp (low)" "dist_emp (high)"],
            xaxis = ("bonds",(-0.15, 0.50)), yaxis= ("density",(0.0,5.0) ) )

Plots.plot(fd.bgrid, [reshape(sol_init.g,1000,2) reshape(sol_transition[2].g,1000,2)], w = 1.5,
            color=[:red :blue :red :blue], style = [:dash :dash :solid :solid], label=["dist_unemp (init)" "dist_emp (init)" "dist_unemp" "dist_emp"],
            xaxis = ("bonds",(-0.15, 0.50)), yaxis= ("density",(0.0,3.0) ) )
Plots.savefig("distribution_t2.pdf")

# t = 2.0
Plots.plot(fd.bgrid, [reshape(sol_init.g,1000,2) reshape(sol_transition[17].g,1000,2)], w = 1.5,
            color=[:red :blue :red :blue], style = [:dash :dash :solid :solid], label=["dist_unemp (init)" "dist_emp (init)" "dist_unemp" "dist_emp"],
            xaxis = ("bonds",(-0.15, 0.50)), yaxis= ("density",(0.0,3.0) ) )
Plots.savefig("distribution_t17.pdf")

# t = 5.0
Plots.plot(fd.bgrid, [reshape(sol_init.g,1000,2) reshape(sol_transition[41].g,1000,2)], w = 1.5,
            color=[:red :blue :red :blue], style = [:dash :dash :solid :solid], label=["dist_unemp (init)" "dist_emp (init)" "dist_unemp (41)" "dist_emp (41)"],
            xaxis = ("bonds",(-0.15, 0.50)), yaxis= ("density",(0.0,3.0) ) )
Plots.savefig("distribution_t41.pdf")


Plots.plot(fd.bgrid, [reshape(sol_end.g,1000,2) reshape(sol_transition[50].g,1000,2)], w = 1.5,
            color=[:red :blue :red :blue], style = [:dash :dash :solid :solid], label=["dist_unemp (end)" "dist_emp (end)" "dist_unemp (50)" "dist_emp (50)"],
            xaxis = ("bonds",(-0.15, 0.50)), yaxis= ("density",(0.0,3.0) ) )

Plots.plot(fd.bgrid, [reshape(sol_init.g,1000,2) reshape(sol_transition[100].g,1000,2)], w = 1.5,
            color=[:red :blue :red :blue], style = [:dash :dash :solid :solid], label=["dist_unemp (init)" "dist_emp (init)" "dist_unemp (25)" "dist_emp (25)"],
            xaxis = ("bonds",(-0.15, 0.50)), yaxis= ("density",(0.0,3.0) ) )

### Transition
time = linspace(-0.25, 50., 403)
Plots.plot(time, [0.0,0.0, R[:,1]...], w = 2.5, xaxis = ("time",(-0.25,10.), 0:1:10), yaxis=("real rate",(-0.02, 0.0325)), label = L"$r(t)$", legendfont = f )
Plots.hline!([rr_end], line=(1.5,:dash, :red) , label = L"$r^{end}$", legendfont = f)
Plots.savefig("transition_realrate.pdf")

# comparing the error on the BOND_mkt
time = linspace(0.0, 50., 401)
Plots.plot(time, BOND_mkt[:,[1,200]], ylabel = "mkt clear res", xaxis = ("time",(0.0,50.)), w = 2., label = [L"$m=1$" L"$m=200$"], legendfont = f)
Plots.savefig("transition_res.pdf")
