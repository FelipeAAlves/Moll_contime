#=
Some continuous time code

@author : Felipe Alves <>

@date: 2017-

References
----------

1)
=#
include("pars.jl")
include("sparse.jl")

# include("init.jl")
# fd = FiniteDiff()
# sol = Solution(fd)

const __pars = Param()
const _tol_hjb = 1e-8

using Roots
using QuantEcon


struct Solution

    ### Values
    g::Vector{Float64}          # Density over state space
    Vv::Vector{Float64}
    V::Array{Float64,2}             # Value function
    cons::Array{Float64,2}          # Optimal consumption
    lsavings::Array{Float64,2}      # Savings
end

struct FiniteDiff

    ### Asset grid
    bgrid::LinSpace{Float64}
    Δb::Float64
    nb::Int64

    ### Stochastic Information
    zgrid::Vector{Float64}
    zmarkov::Matrix{Float64}
    nz::Int64

    ### Storage
    b̃ ::Vector{Float64}                                  # RHS of hjb
    # Λ ::Base.SparseMatrix.SparseMatrixCSC{Float64, Int}  # exogenous movements
    # P ::Base.SparseMatrix.SparseMatrixCSC{Float64, Int}  # storage matrix

    ### Sparse Representations _storage_
    I ::Vector{Int64}
    J ::Vector{Int64}
    Values ::Vector{Float64}

    Iz ::Vector{Int64}
    Jz ::Vector{Int64}
    Valz ::Vector{Float64}
    Nz::Int64

    Ip ::Vector{Int64}
    Jp ::Vector{Int64}
    Valp ::Vector{Float64}
    Np::Int64

    ### Sparse construction _storage_
    csrrowptr ::Vector{Int64}
    # csrcolval ::Vector{Int64}
    # csrnzval ::Vector{Float64}

    csccolptr ::Vector{Int64}
    klasttouch ::Vector{Int64}
    cscrowval ::Vector{Int64}
    cscnzval ::Vector{Float64}

end

function Solution(fd::FiniteDiff)
    @unpack ρ = __pars

    nb = fd.nb; nz = fd.nz;
    Nn = nb * nz;
    bgrid = fd.bgrid; zgrid = fd.zgrid;

    ### Policies & value function
    g = zeros(Nn); Vv = zeros(Nn)
    for (iz,zval) in enumerate(zgrid), (ib,bonds) in enumerate(bgrid)
        ibz = ib + (iz-1)*fd.nb
        Vv[ibz] = 1/ρ * utilfn( 0.03*bonds + zval )
    end
    V = reshape(Vv,nb,nz) ###  WARN WARN
    cons = similar(V); lsavings = similar(V);

    return Solution(g, Vv, V, cons, lsavings)

end

function FiniteDiff(;λ2 = 0.4)

    @unpack ρ, Δhjb = __pars

    ### GRIDs
    zgrid = [.1, .2]; nz = 2
    λ1 = 0.6;
    # λ2 = 0.4;
    zmarkov = [ -λ1 λ1; λ2 -λ2] # generator matrix

    nb = 1000; bmin = -0.15; bmax = 4.;
    bgrid = linspace(bmin,bmax,nb)
    Δb = bgrid[2]-bgrid[1]

    Nn = nz*nb; m = Nn; n = Nn

    ### Sparse matrices

    # for the stochastic transition on RHS,RHS
    lambdas = ( kron(diag(zmarkov),ones(nb)), ones(nb)*λ1, ones(nb)*λ2 )
    Λ = spdiagm(lambdas, (0,nb,-nb))
    Iz, Jz, Valz = findnz(Λ); Nz = length(Iz)

    # for the HJB eq on the LHS,LHS
    ps = ( ρ+1./Δhjb + kron(-diag(zmarkov),ones(nb)) , -ones(nb)*λ1, -ones(nb)*λ2 )
    P = spdiagm(ps, (0,nb,-nb))
    Ip, Jp, Valp = findnz(P); Np = length(Ip)

    I = zeros(Int64,4*Nn-2+2*nb); J = similar(I); #count the diagonal twice
    Values = zeros(Float64,4*Nn-2+2*nb);

    ###  for sparse!

    # Allocate storage for CSR form
    csrrowptr = Vector{Int64}(m+1)
    # csrcolval = Vector{Int64}(coolen)
    # csrnzval = Vector{Tv}(coolen)

    # Allocate storage for the CSC form's column pointers and a necessary workspace
    csccolptr = Vector{Int64}(n+1)
    klasttouch = Vector{Int64}(n)

    # Allocate empty arrays for the CSC form's row and nonzero value arrays
    # The parent method called below automagically resizes these arrays
    cscrowval = Vector{Int64}()
    cscnzval = Vector{Float64}()

    FiniteDiff(bgrid, Δb, nb, zgrid, zmarkov, nz, zeros(Nn),
                # g, Vv, V, c, ls,
                I, J, Values, Iz, Jz, Valz, Nz, Ip, Jp, Valp, Np,
                csrrowptr, csccolptr, klasttouch, cscrowval, cscnzval)

end

#==================================#
##             Utility             ##
#==================================#
utilfn(c) = 1/(1.-__pars.γ) * c^(1.-__pars.γ)
inv_mu(dV) = dV^(-__pars.invγ)

function inv_mu( dV, bond, zval, r)

    c = inv_mu(dV)
    ls   = r * bond + zval - c
    utilV = utilfn(c) + dV*ls

    return c, ls, utilV
end

function inv_mu!(MU, dV, bond, zval, r)

    MU[1] = inv_mu(dV)
    MU[2]   = r * bond + zval - MU[1]
    MU[3] = utilfn(MU[1]) + dV*MU[2]

    return nothing
end

function init_V!(r, sol, fd)

    @unpack ρ = __pars
    zgrid, bgrid = fd.zgrid, fd.bgrid; nb = fd.nb
    rr = max(r,0.001)

    for (iz,zval) in enumerate(zgrid), (ib,bonds) in enumerate(bgrid)
        ibz = ib + (iz-1)*nb
        sol.V[ibz] = 1/ρ * utilfn( rr*bonds + zval )
    end

end

function init_V!(r, V::Matrix{Float64},fd)

    @unpack ρ = __pars
    zgrid, bgrid = fd.zgrid, fd.bgrid; nb = fd.nb

    for (iz,zval) in enumerate(zgrid), (ib,bonds) in enumerate(bgrid)
        V[ib,iz] = 1/ρ * utilfn( r*bonds + zval )
    end

end

## *************************************************************************************
##   Functions
##
##
## =====================================================================================

"""
Find the steady state of the model

#### Arguments

- `arg1`

"""
function solve_steady_state(sol, fd)

    f(r) = steady_state_resid(r, sol, fd)
    return fzero(f, 0.01, 0.04)
end


"""
    steady state bond market clearing residual
"""
function steady_state_resid(rr::Float64, sol::Solution, fd::FiniteDiff)

    #== Solve the HJB ==#
    solve_hjb!(rr, sol, fd)

    #== Solve the Stationary Distribution ==#
    kfe_equation(sol.g, sol.lsavings, fd)

    int_bonds = fd.Δb * repmat(fd.bgrid,fd.nz) .* sol.g
    return sum(int_bonds)
end


"""
    Transition dynamics
    *   Solve the HJB backwards
    *   Solve the KFE forward
"""
function transition_dynamics(R, BOND_mkt, BOND_mkt_clear_dist, sol_init, sol_end, rr_end, fd_end, r_guess = ones(401)*rr_end)

    nb = fd_end.nb; nz = fd_end.nz; Nn = nb*nz
    Δb = fd_end.Δb
    In = speye(Nn)
    Bgrid = repmat(fd_end.bgrid,nz)

    #== Time index ==#
    T = 50; ntime = 401;
    Δt = T/(ntime-1);
    ξ = 10.*ones(ntime);

    #== Create sol_transition ==#
    sol_transition = Vector{Solution}(ntime);
    for ii = 1:ntime
        sol_transition[ii] = Solution(fd_end)
    end
    sol_transition[1].g .= sol_init.g;
    sol_transition[end].Vv .= sol_end.Vv;

    r_transition = zeros( ntime ); r_transition .= r_guess;
    rnew_transition = zeros( ntime ); rnew_transition2 = zeros( ntime );
    bond_mkt = zeros( ntime ); dbond_mkt = zeros( ntime);
    V_out = zeros(Nn); dens = zeros(Nn);

    bond_mkt_clear_dist = 1.0
    it = 1
    while bond_mkt_clear_dist > 1e-5 && it <= 1 # 200

        #== Iterate backwards on the HJB ==#
        R[:,it] = r_transition;
        for tt in ntime-1:-1:1

            V = sol_transition[tt+1].V;
            Vv = sol_transition[tt+1].Vv;

            #== update value function ==#
            updateV!(V_out, V, sol_transition[tt], fd_end, r_transition[tt], Δt)

            distance = maximum(abs, Vv - V_out)
            # (tt%25 == 0) && @printf(     "        value function iteration %d, distance %.4f \n", tt, distance)
            # (tt < 25)    && @printf(     "        value function iteration %d, distance %.4f \n", tt, distance)

            sol_transition[tt].Vv .= V_out ###  WARN WARN update V as well
        end

        #== Iterate forward on the KFE ==#
        for tt in 1:ntime-1

            dens .= sol_transition[tt].g
            g_np1 = sol_transition[tt+1].g

            A_t = kfe_equation(sol_transition[tt].lsavings, fd_end);

            # A_ldiv_B!(lufact(In - Δt*A_t), dens) ###  CAREFUL: updates V_out in-place     ###
            g_np1 .= (In - Δt*A_t)\dens;

            bond_mkt[tt+1] = sum(Δb * Bgrid .* g_np1);
            # (tt%25 == 0) && @printf(     "        kfe period %d, market clear %.4f \n", tt, bond_mkt[tt])k
        end
        bond_mkt_clear_dist = maximum(abs,bond_mkt)

        @printf(     " ********** ITERATION %d: market clear (max) %.2e \n", it, bond_mkt_clear_dist)
        BOND_mkt[:,it] = bond_mkt;
        BOND_mkt_clear_dist[it] = bond_mkt_clear_dist;
        #Forward and backword approximations are used alternately because once
        #the excess savings level becomes small enough, the "rounding error" can
        #lead the update in the interest rate, and repeatedly adding the same
        #rounding error can lead to a propagating error.
        if mod(it,2)==0
            dbond_mkt[1:ntime-1] = bond_mkt[2:ntime]-bond_mkt[1:ntime-1];
            dbond_mkt[ntime] = bond_mkt[ntime]-bond_mkt[ntime-1];
        else
            dbond_mkt[2:ntime] = bond_mkt[2:ntime]-bond_mkt[1:ntime-1];
            dbond_mkt[1] = dbond_mkt[2];
        end
        #Update the interest rate to reduce aggregate savings amount.
        rnew_transition .= r_transition - ξ.*dbond_mkt;
        rnew_transition[ntime-11:ntime] = rnew_transition[ntime-11]; #This was done to minimize the "rounding error" at the end. Otherwise, the end point keeps decreasing
        # Rnew[:,it+1] = rnew_transition;

        #To improve speed, for the first few updates, the update will be done
        #fast, but to reduce wave pattern that gets created, updated interest
        #rate will be smoothed (since the wave patern is due to how things are
        #being updated). After getting "good" initial starting r_t, standard update
        #with declining update weight will be used for convergence. Note that the
        #smoothing will be stopped before the convergence.
        rnew_transition2 .= rnew_transition;
        if it<20
            rnew_transition2[ntime-101:ntime] = rr_end;                 #We know that r will approach the stationary value, so we will force smoothing around the stationary value.
            rnew_transition2[101:ntime] .= QuantEcon.smooth(rnew_transition2,51)[101:ntime];
        elseif it<40
            rnew_transition2[ntime-51:ntime] = rr_end;
            rnew_transition2[101:ntime]  = QuantEcon.smooth(rnew_transition2,15)[101:ntime]; #Amount of smooth will be reduced over time since not as strong of an update is necessary
        elseif it<80
            rnew_transition2[ntime-51:ntime] = rr_end;
            rnew_transition2[101:ntime]  = QuantEcon.smooth(rnew_transition2,7)[101:ntime]; #Amount of smooth will be reduced over time since not as strong of an update is necessary
        elseif it==80
            ξ .= 10.*exp.(-0.0005*collect(1:ntime));
        elseif it==125
            ξ .= 10.*exp.(-0.00005*collect(1:ntime)); #Give bigger update weight to later times.
        end

        #== New iteration ==#
        r_transition .= rnew_transition2;
        it += 1
    end

    return sol_transition
end

"""
    Solve for the HJB equation
    #### Arguments
    - `r`: interest rate
    -
"""
function solve_hjb!(rr::Float64, sol::Solution, fd::FiniteDiff)

    Vv = sol.Vv; V = sol.V
    V_out = zeros(length(V))

    distance =  1.0
    it = 1
    while distance > _tol_hjb && it <= 5000

        #== update value function ==#
        updateV!(V_out, V, sol, fd, rr)

        distance = maximum(abs, V_out - Vv)
        # @printf(     "        value function iteration %d, distance %.4f \n", it, distance)

        Vv .= V_out ###  WARN WARN update V as well
        it += 1
    end
    # (distance < _tol_hjb) && @printf("  Value Function Converged, Iteration %d \n ", it-1)

    return nothing
end

# """
#     vector version
# """
# function updateV!(fd::)
#
# end


"""
    UpdateV! loop version
"""
function updateV!(V_out::Vector{Float64}, V::Matrix{Float64}, sol::Solution, fd::FiniteDiff, rr::Float64, Δ::Float64 = 0.0, LHS::Bool =true)

    @unpack Δhjb,ρ = __pars

    #== Info ==#
    nz = fd.nz; nb = fd.nb; Nn = nb*nz;
    zgrid = fd.zgrid; bgrid = fd.bgrid; Δb = fd.Δb

    #== Sparse matrix values ==#
    I, J, Values = fd.I, fd.J, fd.Values
    # MUb = zeros(3); MUf = zeros(3)

    count = 0
    for (iz,zval) in enumerate(zgrid), (ib,bonds) in enumerate(bgrid)
        ibz = ib + (iz-1)*nb

        #== BACKWARD derivative ==#
        if ib > 1
            dVb = (V[ib, iz] - V[ib-1, iz])/Δb
            cb, lsb, util_b = inv_mu(dVb, bonds, zval, rr)
            # inv_mu!(MUb, dVb, bonds, zval, rr )
            # cb, lsb, util_b = MUb
        else
            lsb = 0.0; util_b = -1e12
        end
        (lsb < 0.0) ? (valid_b = true) : (valid_b = false)

        #== FORWARD derivative ==#
        if ib < nb
            dVf = (V[ib+1, iz] - V[ib, iz])/Δb
            cf, lsf, util_f = inv_mu( dVf, bonds, zval, rr )
            # inv_mu!(MUf, dVf, bonds, zval, rr )
            # cf, lsf, util_f = MUf
        else
            lsf = 0.0; util_f = -1e12
        end
        (lsf > 0.0) ? (valid_f = true) : (valid_f = false)

        #== No savings ==#
        c0     = rr * bonds + zval
        util_0 = utilfn( max(c0,0.0) )  ###  WARN WARN important for very low intereste rates

        #== DECIDE which one to use ==#
        if valid_b & (~valid_f | (util_b>util_f)) & (util_b > util_0)
            sol.cons[ib,iz] = cb
            sol.lsavings[ib,iz] = lsb
        elseif valid_f & (~valid_b | (util_f>util_b) ) & (util_f > util_0)
            sol.cons[ib,iz] = cf
            sol.lsavings[ib,iz] = lsf
        else
            sol.cons[ib,iz] = c0
            sol.lsavings[ib,iz] = 0.0
        end

        lbdrift_b = min(sol.lsavings[ib,iz], 0.0); lbdrift_f = max(sol.lsavings[ib,iz], 0.0)

        if lbdrift_b < 0.
            count += 1
            I[count] = ibz
            J[count] = ibz -1
            Values[count] = -lbdrift_b/Δb
        end

        count += 1
        lval = (lbdrift_b - lbdrift_f)/Δb
        I[count] = ibz
        J[count] = ibz
        Values[count] = lval

        if lbdrift_f > 0.
            count += 1
            I[count] = ibz
            J[count] = ibz+1
            Values[count] = lbdrift_f/Δb
        end

    end

    if LHS

        Values .*= (-1.) # pass it to the LHS
        # println(count)

        if Δ == 0.0

            #== Add the diagonals ==#
            Np = fd.Np;
            Ip, Jp, Valp = fd.Ip, fd.Jp, fd.Valp

            I[count+1:count+Np] .= Ip
            J[count+1:count+Np] .= Jp
            Values[count+1:count+Np] .= Valp

            B_mat = falves_sparse!(I[1:count+Np], J[1:count+Np], Values[1:count+Np], Nn, Nn, +, fd.klasttouch,
                    fd.csrrowptr, #csrcolval, csrnzval,
                    fd.csccolptr, fd.cscrowval, fd.cscnzval)

            ###  SOLVE the linear system
            @. V_out = utilfn( sol.cons[:] ) + 1./Δhjb * V[:]
            A_ldiv_B!(lufact(B_mat), V_out)                         ###  CAREFUL: updates V_out in-place     ###
        else

            #== Add the diagonal ==#
            Nz = fd.Nz;
            Iz, Jz, Valz = fd.Iz, fd.Jz, fd.Valz
            diag_ind = Iz .== Jz ## NOTE indices on the diagonal

            I[count+1:count+Nz] .= Iz
            J[count+1:count+Nz] .= Jz
            Values[count+1:count+Nz] .= -Valz .+ diag_ind.*(ρ + 1.0/Δ)

            B_mat = falves_sparse!(I[1:count+Nz], J[1:count+Nz], Values[1:count+Nz], Nn, Nn, +, fd.klasttouch,
                    fd.csrrowptr, #csrcolval, csrnzval,
                    fd.csccolptr, fd.cscrowval, fd.cscnzval)

            ###  SOLVE the linear system
            @. V_out = utilfn( sol.cons[:] ) + 1./Δ * V[:]
            A_ldiv_B!(lufact(B_mat), V_out)                         ###  CAREFUL: updates V_out in-place     ###
        end
        fill!(I,0); fill!(J,0); fill!(Values,0.0)
    end

    return nothing
end


"""
Solves for the invariant distribution

#### Arguments

- `arg1`

"""
function kfe_equation(lsavings, fd)

    #== Info ==#
    nz = fd.nz; nb = fd.nb; Nn = nb*nz;
    zgrid = fd.zgrid;
    bgrid = fd.bgrid; Δb = fd.Δb

    #== Sparse matrix values ==#
    I, J, Values = fd.I, fd.J, fd.Values
    fill!(I,0); fill!(J,0); fill!(Values,0.0)

    count = 0
    for (iz,zval) in enumerate(zgrid), (ib,bonds) in enumerate(bgrid)
        ibz = ib + (iz-1)*nb

        lbdrift_b = min(lsavings[ib,iz], 0.0); lbdrift_f = max(lsavings[ib,iz], 0.0)

        # i-1
        if lbdrift_b < 0.
            count += 1
            I[count] = ibz
            J[count] = ibz -1
            Values[count] = -lbdrift_b/Δb
        end

        # i
        count += 1
        lval = (lbdrift_b - lbdrift_f)/Δb
        I[count] = ibz
        J[count] = ibz
        Values[count] = lval

        # i+1
        if lbdrift_f > 0.
            count += 1
            I[count] = ibz
            J[count] = ibz+1
            Values[count] = lbdrift_f/Δb
        end
    end

    #== Add the diagonals ONLY for the stochastic part ==#
    Nz = fd.Nz;
    Iz, Jz, Valz = fd.Iz, fd.Jz, fd.Valz

    I[count+1:count+Nz] .= Iz
    J[count+1:count+Nz] .= Jz
    Values[count+1:count+Nz] .= Valz

    ###  WARN invert J,I indexes
    A_transpose = falves_sparse!(J[1:count+Nz], I[1:count+Nz] , Values[1:count+Nz], Nn, Nn, +, fd.klasttouch,
            fd.csrrowptr, #csrcolval, csrnzval,
            fd.csccolptr, fd.cscrowval, fd.cscnzval)

    return A_transpose

end
function kfe_equation(dens, lsavings, fd)

    #== Info ==#
    nz = fd.nz; nb = fd.nb; Nn = nb*nz;
    zgrid = fd.zgrid;
    bgrid = fd.bgrid; Δb = fd.Δb

    #== Sparse matrix values ==#
    I, J, Values = fd.I, fd.J, fd.Values
    fill!(I,0); fill!(J,0); fill!(Values,0.0)

    count = 0
    for (iz,zval) in enumerate(zgrid), (ib,bonds) in enumerate(bgrid)
        ibz = ib + (iz-1)*nb

        lbdrift_b = min(lsavings[ib,iz], 0.0); lbdrift_f = max(lsavings[ib,iz], 0.0)

        # i-1
        if lbdrift_b < 0.
            count += 1
            I[count] = ibz
            J[count] = ibz -1
            Values[count] = -lbdrift_b/Δb
        end

        # i
        count += 1
        lval = (lbdrift_b - lbdrift_f)/Δb
        I[count] = ibz
        J[count] = ibz
        Values[count] = lval

        # i+1
        if lbdrift_f > 0.
            count += 1
            I[count] = ibz
            J[count] = ibz+1
            Values[count] = lbdrift_f/Δb
        end
    end

    #== Add the diagonals ONLY for the stochastic part ==#
    Nz = fd.Nz;
    Iz, Jz, Valz = fd.Iz, fd.Jz, fd.Valz

    I[count+1:count+Nz] .= Iz
    J[count+1:count+Nz] .= Jz
    Values[count+1:count+Nz] .= Valz

    ###  WARN invert J,I indexes
    A_transpose = falves_sparse!(J[1:count+Nz], I[1:count+Nz] , Values[1:count+Nz], Nn, Nn, +, fd.klasttouch,
            fd.csrrowptr, #csrcolval, csrnzval,
            fd.csccolptr, fd.cscrowval, fd.cscnzval)


        #== SOLVE for the stationary distribution ==#
        dens[:] = 0.0; dens[1] = 1.0
        A_transpose[1,:] .= 0.0; A_transpose[1,1] = 1.0
        A_ldiv_B!(lufact(A_transpose), dens) ###  CAREFUL: updates V_out in-place     ###

        ###  NOTE:  distribution is solved in one-step onlu     ###

        #== normalization ==#
        mass = Δb*sum(dens)
        dens .= dens/mass
end
