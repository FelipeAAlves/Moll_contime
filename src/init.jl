
#=
Some continuous time code

@author : Felipe Alves <>

@date: 2017-

References
----------

    1)
=#
const __pars = Param()


immutable FiniteDiff

    ### asset grid
    bgrid::LinSpace{Float64}
    Δb::Float64
    nb::Int64

    ### Stochastic Information
    zgrid::Vector{Float64}
    zmarkov::Matrix{Float64}
    nz::Int64

    # zdist::Vector{Float64}
    # zmarkov_diag::Vector{Float64}
    # zmarkov_off::Matrix{Float64}

    ### Storage
    b̃ ::Vector{Float64}                                  # RHS of hjb
    # Λ ::Base.SparseMatrix.SparseMatrixCSC{Float64, Int}  # exogenous movements
    # P ::Base.SparseMatrix.SparseMatrixCSC{Float64, Int}  # storage matrix

    g::Vector{Float64}          # Density over state space

    Vv::Vector{Float64}
    V::Array{Float64,2}             # Value function
    cons::Array{Float64,2}          # Optimal consumption
    lsavings::Array{Float64,2}      # Savings

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

function FiniteDiff()

    @unpack ρ, r, Δhjb = __pars

    ### GRIDs
    zgrid = [.1, .2]; nz = 2
    λ1 = 0.02;  λ2 = 0.03;
    zmarkov = [ -λ1 λ1; λ2 -λ2]

    nb = 500; bmin = -0.02; bmax = 2;
    bgrid = linspace(bmin,bmax,nb)
    Δb = bgrid[2]-bgrid[1]

    Nn = nz*nb; m = Nn; n = Nn

    ### Policies & value function
    g = zeros(Nn); b̃ = zeros(Nn); Vv = zeros(Nn)
    for (iz,zval) in enumerate(zgrid), (ib,bonds) in enumerate(bgrid)
        ibz = ib + (iz-1)*nb
        Vv[ibz] = 1/ρ * utilfn( r*bonds + zval )
    end
    V = reshape(Vv,na,nz)
    c = similar(V); ls = similar(V);


    ### Sparse matrices
    lambdas = ( kron(diag(zmarkov),ones(nb)), ones(nb)*λ1, ones(nb)*λ2 )
    Λ = spdiagm(lambdas, (0,nb,-nb))
    Iz, Jz, Valz = findnz(Λ); Nz = length(Iz)

    ps = ( kron(-diag(zmarkov)+ρ+1/Δhjb,ones(nb)) , ones(nb)*-λ1, ones(nb)*-λ2 )
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

    FiniteDiff(bgrid, Δb, nb, zgrid, zmarkov, nz, b̃, g, Vv, V, c, ls, I, J, Values, Iz, Jz, Valz, Nz, Ip, Jp, Valp, Np,
               csrrowptr, csccolptr, klasttouch, cscrowval, cscnzval)

end

function init_V!(V::Vector{Float64},fd)

    @unpack ρ,r = __pars
    zgrid, bgrid = fd.zgrid, fd.bgrid; nb = fd.nb

    for (iz,zval) in enumerate(zgrid), (ib,bonds) in enumerate(bgrid)
        ibz = ib + (iz-1)*nb
        V[ibz] = 1/ρ * utilfn( r*bonds + zval )
    end

end

function init_V!(V::Matrix{Float64},fd)

    @unpack ρ,r = __pars
    zgrid, bgrid = fd.zgrid, fd.bgrid; nb = fd.nb

    for (iz,zval) in enumerate(zgrid), (ib,bonds) in enumerate(bgrid)
        V[ib,iz] = 1/ρ * utilfn( r*bonds + zval )
    end

end


utilfn(c) = 1/(1.-__pars.γ) * c^(1.-__pars.γ)
inv_mu(dV) = dV^(-__pars.invγ)

function inv_mu( dV, bond, zval )
    @unpack r = __pars

    c = inv_mu(dV)
    ls   = r * bond + zval - c
    utilV = utilfn(c) + dV*ls

    return c, ls, utilV
end

function inv_mu!(MU, dV, bond, zval )
    @unpack r = __pars

    MU[1] = inv_mu(dV)
    MU[2]   = r * bond + zval - MU[1]
    MU[3] = utilfn(MU[1]) + dV*MU[2]

    return nothing
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
function solve_steady_state()

end


"""
    Solve for the HJB equation
"""
function solve_hjb(fd::FiniteDiff)

    Vv = fd.Vv; V = fd.V
    V_out = zeros(length(V))

    dist =  1.0
    it = 1
    while distance < 1e-6 && it <= 1000

        updateV!(V_out, V, fd)

        distance = maxabs( V_out - Vv)
        @printf(     "        value function iteration %d, distance %.4f \n", it, distance)
        Vv .= V_out

        (distance < 1e-6) && @printf("Value Function Converged, Iteration %d ", it)
        it += 1
    end

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
function updateV!(V_out::Vector{Float64}, V::Matrix{Float64}, fd::FiniteDiff, LHS::Bool =true)

    @unpack r, Δhjb = __pars

    #== Info ==#
    nz = fd.nz; nb = fd.nb; Nn = nb*nz;
    zgrid = fd.zgrid;
    bgrid = fd.bgrid; Δb = fd.Δb

    #== Sparse matrix values ==#
    I, J, Values = fd.I, fd.J, fd.Values
    # MUb = zeros(3); MUf = zeros(3)

    count = 0
    for (iz,zval) in enumerate(zgrid), (ib,bonds) in enumerate(bgrid)
        ibz = ib + (iz-1)*nb

        #== Backward derivative ==#
        if ib > 1
            dVb = (V[ib, iz] - V[ib-1, iz])/Δb
            cb, lsb, util_b = inv_mu(dVb, bonds, zval)
            # inv_mu!(MUb, dVb, bonds, zval )
            # cb, lsb, util_b = MUb
        else
            lsb = 0.0; util_b = -1e12
        end
        (lsb < 0.0) ? (valid_b = true) : (valid_b = false)

        #== Forward derivative ==#
        if ib < nb
            dVf = (V[ib+1, iz] - V[ib, iz])/Δb
            cf, lsf, util_f = inv_mu( dVf, bonds, zval )
            # inv_mu!(MUf, dVf, bonds, zval )
            # cf, lsf, util_f = MUf
        else
            lsf = 0.0; util_f = -1e12
        end
        (lsf > 0.0) ? (valid_f = true) : (valid_f = false)

        #== No savings ==#
        c0     = r * bonds + zval
        util_0 = utilfn(c0)

        #== DECIDE which one to use ==#
        if valid_b & (~valid_f | (util_b>util_f)) & (util_b > util_0)
            fd.cons[ib,iz] = cb
            fd.lsavings[ib,iz] = lsb
        elseif valid_f & (~valid_b | (util_f>util_b) ) & (util_f > util_0)
            fd.cons[ib,iz] = cf
            fd.lsavings[ib,iz] = lsf
        else
            fd.cons[ib,iz] = c0
            fd.lsavings[ib,iz] = 0.0
        end

        lbdrift_b = min(fd.lsavings[ib,iz], 0.0); lbdrift_f = max(fd.lsavings[ib,iz], 0.0)

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

        #== Add the diagonals ==#
        Np = fd.Np;
        Ip, Jp, Valp = fd.Ip, fd.Jp, fd.Valp

        I[count+1:count+Np] .= Ip
        J[count+1:count+Np] .= Jp
        Values[count+1:count+Np] .= Valp

        B_mat = falves_sparse!(I[1:count+Np], J[1:count+Np], Values[1:count+Np], Nn, Nn, +, fd.klasttouch,
                fd.csrrowptr, #csrcolval, csrnzval,
                fd.csccolptr, fd.cscrowval, fd.cscnzval)

        # nn = length(fd.cscrowval)
        # deleteat!(fd.cscrowval,1:nn); deleteat!(fd.cscnzval,1:nn);

        ###  SOLVE the linear system
            ###  CASE 01
            @. V_out = utilfn( vec(fd.cons) ) + 1/Δhjb * fd.Vv
            A_ldiv_B!(lufact(B_mat), V_out)                         ###  CAREFUL: updates V_out in-place     ###

        fill!(I,0); fill!(J,0); fill!(Values,0.0)
    end

    return nothing
end


"""
Solves for the invariant distribution

#### Arguments

- `arg1`

"""
function solve_kfe(dens, lsavings, fd)

    @unpack r, Δhjb = __pars

    #== Info ==#
    nz = fd.nz; nb = fd.nb; Nn = nb*nz;
    zgrid = fd.zgrid;
    bgrid = fd.bgrid; Δb = fd.Δb

    count = 0
    for (iz,zval) in enumerate(zgrid), (ib,bonds) in enumerate(bgrid)
        ibz = ib + (iz-1)*nb

        lbdrift_b = min(fd.lsavings[ib,iz], 0.0); lbdrift_f = max(fd.lsavings[ib,iz], 0.0)

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

    #== Add the diagonals ==#
    Nz = fd.Nz;
    Iz, Jz, Valz = fd.Iz, fd.Jz, fd.Valz

    I[count+1:count+Nz] .= Iz
    J[count+1:count+Nz] .= Jz
    Values[count+1:count+Nz] .= Valz

    A_transpose = falves_sparse!(J[1:count+Nz], I[1:count+Nz] , Values[1:count+Nz], Nn, Nn, +, fd.klasttouch,
            fd.csrrowptr, #csrcolval, csrnzval,
            fd.csccolptr, fd.cscrowval, fd.


    #== SOLVE for the stationary distribution ==#
    dens[:] = 0.0; dens[1] = 1.0
    A_transpose[1,:] .= 0.0; A[1,1] = 1.0
    A_ldiv_B!(lufact(A_transpose), dens)                         ###  CAREFUL: updates V_out in-place     ###

    return A_transpose

end
