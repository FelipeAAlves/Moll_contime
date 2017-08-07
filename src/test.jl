

import Plots
Plots.pyplot()

include("pars.jl")
include("labor.jl")
include("sparse.jl")
include("kfe.jl")

# x = linspace(0,15.0)
# Plots.plot(x, p.(x))

fd = FiniteDiff(rr = 0.005)
solve_equil(0.005, 0.25, fd)

reshape(fd.search_effort[2],150,61)
At = kfe_stationarity!(fd.g, fd.lsavings, fd.wage_search_pol, fd.tightness, fd.prob_job, fd);
At_long = kfe_stationarity_long!(fd.glong, fd.lsavings, fd.wage_search_pol, fd.tightness, fd.prob_job, fd);


## density over wages
@unpack na, agrid, nw, wgrid, wgrid_long, wmin, wmax, wmaximum, Δw, wgrid_, Δtildea = _grid_
dens_wage = sum(fd.g_mat .* Δtildea,1)[:];
dens_wage_long = sum(fd.glong_mat .* Δtildea,1)[:];

Plots.plot(wgrid, dens_wage[2:end], line = 2, xaxis = ("wage", wgrid) )
Plots.plot!(wgrid_long, dens_wage_long[2:end], line = 2, xaxis = ("wage", wgrid) )
Plots.savefig("labor_output/wage_dens_with_no_search.pdf")

fd.g_mat .* (Δtildea .* [1.0,Δw*nes(61)...]')

## conditional density of assets
dens_asset_condwage = zeros(fd.g_mat)
for iw = 1:nw+1
     dens_wage[iw] > 1e-5 && (dens_asset_condwage[:,iw] = fd.g_mat[:,iw] / dens_wage[iw])
end
Plots.plot(agrid, dens_asset_condwage[:,[1, 20:3:59...]], line = 2, ylim = (0.0,5.0))


Plots.plot(agrid, sum(fd.g_mat .* [1.0, Δw*ones(61)...]',2), line = 2, ylim = (0.0,5.0))


@unpack na, agrid, nw, wgrid, wmin, wmax, wmaximum, Δw, wgrid_ = _grid_
dens = zeros(na*(nw+1), 250)
dens[1] = 1.0 ###  WARN:  mass at first node     ###

At = kfeupdate!(dens, fd.lsavings, fd.wage_search_pol, fd.tightness, fd.prob_job, fd)
dens_wage = zeros(nw+1, 250)
for t = 1 :250
    dens_wage[:,t] = sum(reshape(dens[:,t], na, nw+1) .* _grid_.Δtildea,1)[:]
end
ind = (0:20:240).+1
Plots.plot(wgrid, dens_wage[2:end,ind[1]], line = 2, xaxis = ("wage", wgrid_) )
Plots.plot!(wgrid, dens_wage[2:end,ind[3]], line = 2, xaxis = ("wage", wgrid_) )



fd.g_mat .* (_grid_.Δtildea .* [1.0,_grid_.Δw*ones(61)...]')
dens_wage = sum(fd.g_mat .* _grid_.Δtildea,1)[:]
dens_asset_condwage = (fd.g_mat ./ dens_wage')

Plots.plot(wgrid, dens_wage[2:end], line = 2, xaxis = ("wage", wgrid_) )
Plots.plot(agrid, fd.g_mat[:,[1, 20, 25, 30, 40]], line = 2, ylim = (0.0,0.45))


wvals = linspace(wmin, wmax, 1601)
Plots.plot(wgrid, p.(fd.Θ[1,:]), line = :scatter)
Plots.plot!(wvals,  [p.(fd.tightness.(1:5,wvals)) fd.prob_job.(1:5,wvals)] , line = 2, xaxis = ("wage", wgrid_))
Plots.savefig("labor_output/firm_value_l0.23_81.pdf")
Plots.plot(wvals,  p.(fd.tightness.(0:25:125+1,wvals'))' , line = 2, xaxis = ("wage", wgrid_))

### Wage search
@unpack na, agrid, nw, wgrid, wmin, wmax, wmaximum, Δw = _grid_
Plots.plot(agrid, [fd.wage_search_pol[1] reshape(fd.wage_search_pol[2],na,nw)], yaxis = ("wage", linspace(0.70,0.90,41)), ylim = (0.75,0.905) )
Plots.savefig("labor_output/wage_search_l0.23_61.pdf")
Plots.plot(agrid, reshape(fd.wage_search_pol[2],na,nw)[:,1:10], yaxis = ("wage",wgrid) )


fd.prob_job.(1,wgrid') .* (fd.Ve[1,:]-fd.Ve[1,1])
# --------------------------------------------------------------------------------------
# % %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%%
# --------------------------------------------------------------------------------------

# fd.g .= fd.g/sum(fd.g
# fd.g.*repmat(Δtildea,52)
# reshape(fd.g,na,nw)
#
# ind = vcat(collect(1:3:31),collect(32:40))
# spl = Spline2D(1:na, wgrid[ind], fd.Θ[:, ind], ky =1)
# Plots.plot!(wvals,  [p.(fd.tightness.(10,wvals)) fd.prob_job.(10,wvals) spl.(10,wvals)] , line = 2)

@unpack na, agrid, nw, wgrid, wmin, wmax, wmaximum, Δw, wgrid_ = _grid_
wvals = linspace(wmin, wmax, 1601)
Plots.plot(wgrid, p.(fd.Θ[1,:]), line = :scatter)
Plots.plot!(wvals,  [p.(fd.tightness.(1,wvals)) fd.prob_job.(1,wvals)] , line = 2, xaxis = ("wage", wgrid_))
# Plots.plot!(wvals,  [p.(tightness.(1,wvals)) prob_job.(1,wvals)] , line = 2, xaxis = ("wage", wgrid_))

# value of the firm
Plots.plot(wgrid, [(0.90 - wgrid')/(0.005+__pars.δ); reshape(fd.Jfirm,na,nw)[1:25:na,:]]', yaxis=("val",(0.0,3.0)), xaxis = ("wage",wgrid) )

# Plots.plot(wgrid, p.(fd.Θ[1,:]), line = :scatter)
# Plots.plot!(wvals,  [p.(Θ_slp.(1,wvals)) p.(fd.tightness.(1,wvals))] , line = 2)
#
# Plots.plot(wgrid, fd.Θ[10,:], line = :scatter)
# Plots.plot!(wvals,  [fd.tightness.(10,wvals) spl.(10,wvals)], line = 2)

# Plots.plot!(w, [prob_job.(10,w) prob_job_quad.(10,w)] , line=([:path :path],3) )

Plots.plot(linspace(0.60,0.9), fd.prob_job.(10, linspace(0.6,0.9)))

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

@unpack na, agrid, nw, wgrid, wmin, wmax, wmaximum, Δw = _grid_
Plots.plot(agrid, [fd.wage_search_pol[1] reshape(fd.wage_search_pol[2],na,nw)[:,1:10]], yaxis = ("wage",wgrid) )

Plots.plot(agrid, reshape(fd.wage_search_pol[2],na,nw)[:,1:10], yaxis = ("wage",wgrid) )

max_wage1!(fd.wage_search_pol2, fd.R, fd.V[1], fd.Ve, fd.tightness, fd.prob_job)
Plots.plot(agrid, [fd.wage_search_pol[1] fd.wage_search_pol2[1]], yaxis = ("wage",wgrid))
Plots.plot(agrid, [reshape(fd.wage_search_pol[2],na,nw)[:,1] reshape(fd.wage_search_pol2[2],na,nw)[:,1]], yaxis = ("wage",wgrid) )

# Plots.plot(agrid, [fd.wage_search_pol[1] reshape(fd.wage_search_pol[2],na,nw)], yaxis = ("wage",wgrid) )
# Plots.plot(agrid, fd.wage_search_pol[1]) #, yaxis=("wage",(0.75, 0.875)) )
# Plots.plot(agrid, [fd.wage_search_pol[1] reshape(fd.wage_search_pol[2],na,nw)[:,1:10]], yaxis = ("wage",wgrid) )

# @unpack na, agrid, nw, wgrid, wmin, wmax, wmaximum, Δw = _grid_
# Plots.plot(agrid, [fd.cons[1] reshape(fd.cons[2][1:25*na],na,25)] ,yaxis = ("c",(0.25,1.00)))
# Plots.plot(agrid, [fd.cons0[1] reshape(fd.cons0[2][1:25*na],na,25)],yaxis = ("c",(0.25,1.00)))
# Plots.plot(agrid, [fd.lsavings[1] reshape(fd.lsavings[2],na,nw)], yaxis = ("sav",(-0.30,0.25)))

# Plots.plot(agrid, reshape(fd.lsavings0[2],na,nw), yaxis = ("sav",(-0.30,0.25)))
# Plots.plot(agrid, reshape(fd.lsavings[2][1:25*na],na,25))
# Plots.plot(agrid, [fd.V[1] reshape(fd.V[2][1:10*na],na,10)] )
# Plots.plot(agrid, [fd.V0[1] reshape(fd.V0[2][1:10*na],na,10)] )


# Plots.plot(agrid, reshape(fd.V[2],na,nw) )
# Plots.plot(agrid, fd.cons[1] )
# Plots.plot(agrid, fd.lsavings[1] )
#
# reshape(fd.lsavings[2],na,nw)
# Plots.plot(agrid, fd.wage_search_pol[1], yaxis=("wage",(0.80, 0.875)) )

#
#
# Jfirm_itp = interpolate( reshape(fd.Jfirm,na,nw), (NoInterp(),BSpline(Cubic(Natural()))), OnGrid() );
# Jfirm_stp = scale(Jfirm_itp, 1:na, StepRangeLen(wgrid))
#

# w_search = fd.wage_search_pol
# Vu = fd.V[1]
# Ve = fd.Ve
# tightness = fd.tightness
#
# wmaximum_ = wmaximum[1]
# #== interpolate ==#
# Ve_stp  = interpolate((1:na, wgrid), Ve, (NoInterp(), Gridded(Linear())) )
# Vef(ia::Int64, wtilde) = Ve_stp[ia, wtilde]
#
# function f_objective(w_search, ia, val_cont)
#     # w_search = wvec[1]
#     Θ = tightness(ia, w_search)
#
#     return - p(Θ) * ( Vef(ia, w_search) - val_cont )
# end
#
# y_vals = Array{Vector{Float64}}(0)
# iw = 4; wage = wgrid[iw]
# wvals = linspace(wage, wmax,150);
# for ia = 75:100
#     asset = agrid[ia]
#     # y = [f_objective(w, ia, Ve[ia,iw]) for w in wvals]; push!(y_vals,y)
#     y = [f_objective(w, ia, Vu[ia]) for w in wvals]; push!(y_vals,y)
# end
# #
# Plots.plot(wvals, y_vals, xaxis=("wage",(0.60,0.85)) )#,yaxis=("wage",(-3,-1.5)))
