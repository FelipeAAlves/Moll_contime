fd  = FiniteDiff(λ2 = 0.4); sol = Solution(fd);
fd_init  = FiniteDiff(λ2 = 0.4); sol_init = Solution(fd_init);
fd_end  = FiniteDiff(λ2 = 0.8); sol_end = Solution(fd_end);

r_init = solve_steady_state(sol_init, fd_init)
rr_end = solve_steady_state(sol_end, fd_end)
BOND_mkt          = zeros(401, 201);
BOND_mkt_clear_dist = zeros(201);
R = zeros(401, 201);
Rnew = zeros(401, 201);


**************************************************************************************
  PLOTS

======================================================================================
writedlm("output/r_guess.txt", r_guess)
r_guess = squeeze(readdlm("output/r_guess.txt"),2)
sol_transition = transition_dynamics(R, BOND_mkt, BOND_mkt_clear_dist, sol_init, sol_end, rr_end, fd_end, r_guess)
#
### Consumption
Plots.plot(fd.bgrid, [sol_init.cons sol_end.cons], w = 1.5,
            color=[:red :blue :red :blue], style = [:dash :dash :solid :solid], label=["c_unemp (low)" "c_emp (low)" "c_unemp (high)" "c_emp (high)"],
            xaxis = ("bonds",(-0.15, 1.5)), yaxis = ("c", (0.05,0.30)))
Plots.savefig("cons_policy.pdf")

## Distributions

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


time = linspace(-0.25, 50., 403)
Plots.plot(time, [0.0,0.0,R[:,200]...], w = 2, xaxis = ("time",(-0.25,10.), 0:1:10), yaxis=("real rate",(-0.02, 0.0325)) )
Plots.hline!([rr_end], line=(1.0,:dash, :red) , label = "")

# comparing the error on the BOND_mkt
time = linspace(0.0, 50., 401)
Plots.plot(time, BOND_mkt[:,[1,200]], xaxis = ("time",(0.0,50.)), w = 2.)
