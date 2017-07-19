"""
    look at the Julia sparse function
"""
function falves_sparse!{Tv,Ti<:Integer}(I::AbstractVector{Ti}, J::AbstractVector{Ti}, V::AbstractVector{Tv},
        m::Integer, n::Integer, combine, klasttouch::Vector{Ti},
        csrrowptr::Vector{Ti}, #csrcolval::Vector{Ti}, csrnzval::Vector{Tv},
        csccolptr::Vector{Ti}, cscrowval::Vector{Ti}, cscnzval::Vector{Tv})

    nI = length(I)
    csrcolval = Vector{Int64}(nI)
    csrnzval  = Vector{Float64}(nI)

    # Compute the CSR form's row counts and store them shifted forward by one in csrrowptr
    fill!(csrrowptr, 0)
    coolen = length(I)
    @inbounds for k in 1:coolen
        Ik = I[k]
        if 1 > Ik || m < Ik
            throw(ArgumentError("row indices I[k] must satisfy 1 <= I[k] <= m"))
        end
        csrrowptr[Ik+1] += 1
    end

    # Compute the CSR form's rowptrs and store them shifted forward by one in csrrowptr
    countsum = 1
    csrrowptr[1] = 1
    @inbounds for i in 2:(m+1)
        overwritten = csrrowptr[i]
        csrrowptr[i] = countsum
        countsum += overwritten
    end

    # Counting-sort the column and nonzero values from J and V into csrcolval and csrnzval
    # Tracking write positions in csrrowptr corrects the row pointers
    @inbounds for k in 1:coolen
        Ik, Jk = I[k], J[k]
        if 1 > Jk || n < Jk
            throw(ArgumentError("column indices J[k] must satisfy 1 <= J[k] <= n"))
        end
        csrk = csrrowptr[Ik+1]
        csrrowptr[Ik+1] = csrk+1
        csrcolval[csrk] = Jk
        csrnzval[csrk] = V[k]
    end
    # This completes the unsorted-row, has-repeats CSR form's construction

    # Sweep through the CSR form, simultaneously (1) caculating the CSC form's column
    # counts and storing them shifted forward by one in csccolptr; (2) detecting repeated
    # entries; and (3) repacking the CSR form with the repeated entries combined.
    #
    # Minimizing extraneous communication and nonlocality of reference, primarily by using
    # only a single auxiliary array in this step, is the key to this method's performance.
    fill!(csccolptr, 0)
    fill!(klasttouch, 0)
    writek = 1
    newcsrrowptri = 1
    origcsrrowptri = 1
    origcsrrowptrip1 = csrrowptr[2]
    @inbounds for i in 1:m
        for readk in origcsrrowptri:(origcsrrowptrip1-1)
            j = csrcolval[readk]
            if klasttouch[j] < newcsrrowptri
                klasttouch[j] = writek
                if writek != readk
                    csrcolval[writek] = j
                    csrnzval[writek] = csrnzval[readk]
                end
                writek += 1
                csccolptr[j+1] += 1
            else
                klt = klasttouch[j]
                csrnzval[klt] = combine(csrnzval[klt], csrnzval[readk])
            end
        end
        newcsrrowptri = writek
        origcsrrowptri = origcsrrowptrip1
        origcsrrowptrip1 != writek && (csrrowptr[i+1] = writek)
        i < m && (origcsrrowptrip1 = csrrowptr[i+2])
    end

    # Compute the CSC form's colptrs and store them shifted forward by one in csccolptr
    countsum = 1
    csccolptr[1] = 1
    @inbounds for j in 2:(n+1)
        overwritten = csccolptr[j]
        csccolptr[j] = countsum
        countsum += overwritten
    end

    # Now knowing the CSC form's entry count, resize cscrowval and cscnzval if necessary
    cscnnz = countsum - 1
    length(cscrowval) < cscnnz && resize!(cscrowval, cscnnz)
    length(cscnzval) < cscnnz && resize!(cscnzval, cscnnz)

    # Finally counting-sort the row and nonzero values from the CSR form into cscrowval and
    # cscnzval. Tracking write positions in csccolptr corrects the column pointers.
    @inbounds for i in 1:m
        for csrk in csrrowptr[i]:(csrrowptr[i+1]-1)
            j = csrcolval[csrk]
            x = csrnzval[csrk]
            csck = csccolptr[j+1]
            csccolptr[j+1] = csck+1
            cscrowval[csck] = i
            cscnzval[csck] = x
        end
    end

    SparseMatrixCSC(n, m, csrrowptr, csrcolval, csrnzval)
end

function falves_sparse!{Tv,Ti<:Integer}(I::AbstractVector{Ti}, J::AbstractVector{Ti}, V::AbstractVector{Tv},
        m::Integer, n::Integer, combine, klasttouch::Vector{Ti},
        csrrowptr::Vector{Ti}, csrcolval::Vector{Ti}, csrnzval::Vector{Tv},
        csccolptr::Vector{Ti}, cscrowval::Vector{Ti}, cscnzval::Vector{Tv})

    # nI = length(I)
    # csrcolval = Vector{Int64}(nI)
    # csrnzval  = Vector{Float64}(nI)

    # Compute the CSR form's row counts and store them shifted forward by one in csrrowptr
    fill!(csrrowptr, 0)
    coolen = length(I)
    @inbounds for k in 1:coolen
        Ik = I[k]
        if 1 > Ik || m < Ik
            throw(ArgumentError("row indices I[k] must satisfy 1 <= I[k] <= m"))
        end
        csrrowptr[Ik+1] += 1
    end

    # Compute the CSR form's rowptrs and store them shifted forward by one in csrrowptr
    countsum = 1
    csrrowptr[1] = 1
    @inbounds for i in 2:(m+1)
        overwritten = csrrowptr[i]
        csrrowptr[i] = countsum
        countsum += overwritten
    end

    # Counting-sort the column and nonzero values from J and V into csrcolval and csrnzval
    # Tracking write positions in csrrowptr corrects the row pointers
    @inbounds for k in 1:coolen
        Ik, Jk = I[k], J[k]
        if 1 > Jk || n < Jk
            throw(ArgumentError("column indices J[k] must satisfy 1 <= J[k] <= n"))
        end
        csrk = csrrowptr[Ik+1]
        csrrowptr[Ik+1] = csrk+1
        csrcolval[csrk] = Jk
        csrnzval[csrk] = V[k]
    end
    # This completes the unsorted-row, has-repeats CSR form's construction

    # Sweep through the CSR form, simultaneously (1) caculating the CSC form's column
    # counts and storing them shifted forward by one in csccolptr; (2) detecting repeated
    # entries; and (3) repacking the CSR form with the repeated entries combined.
    #
    # Minimizing extraneous communication and nonlocality of reference, primarily by using
    # only a single auxiliary array in this step, is the key to this method's performance.
    fill!(csccolptr, 0)
    fill!(klasttouch, 0)
    writek = 1
    newcsrrowptri = 1
    origcsrrowptri = 1
    origcsrrowptrip1 = csrrowptr[2]
    @inbounds for i in 1:m
        for readk in origcsrrowptri:(origcsrrowptrip1-1)
            j = csrcolval[readk]
            if klasttouch[j] < newcsrrowptri
                klasttouch[j] = writek
                if writek != readk
                    csrcolval[writek] = j
                    csrnzval[writek] = csrnzval[readk]
                end
                writek += 1
                csccolptr[j+1] += 1
            else
                klt = klasttouch[j]
                csrnzval[klt] = combine(csrnzval[klt], csrnzval[readk])
            end
        end
        newcsrrowptri = writek
        origcsrrowptri = origcsrrowptrip1
        origcsrrowptrip1 != writek && (csrrowptr[i+1] = writek)
        i < m && (origcsrrowptrip1 = csrrowptr[i+2])
    end

    # Compute the CSC form's colptrs and store them shifted forward by one in csccolptr
    countsum = 1
    csccolptr[1] = 1
    @inbounds for j in 2:(n+1)
        overwritten = csccolptr[j]
        csccolptr[j] = countsum
        countsum += overwritten
    end

    # Now knowing the CSC form's entry count, resize cscrowval and cscnzval if necessary
    cscnnz = countsum - 1
    length(cscrowval) < cscnnz && resize!(cscrowval, cscnnz)
    length(cscnzval) < cscnnz && resize!(cscnzval, cscnnz)

    # Finally counting-sort the row and nonzero values from the CSR form into cscrowval and
    # cscnzval. Tracking write positions in csccolptr corrects the column pointers.
    @inbounds for i in 1:m
        for csrk in csrrowptr[i]:(csrrowptr[i+1]-1)
            j = csrcolval[csrk]
            x = csrnzval[csrk]
            csck = csccolptr[j+1]
            csccolptr[j+1] = csck+1
            cscrowval[csck] = i
            cscnzval[csck] = x
        end
    end

    # SparseMatrixCSC(m, n, csccolptr, cscrowval, cscnzval)
end
