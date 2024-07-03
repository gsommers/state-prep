#!/usr/licensed/julia/1.10/bin/julia

#= simulate bipartite entropy on tree with Clifford gates =#
module CircuitTools

using QuantumClifford

export transversal_cnot!, trace_layer!, remove_qubits!, track_substabilizers

function transversal_cnot!(s, sites, measure_type)
    for site_pair in sites
        if measure_type==1 # Z checks: ancilla is target
	    apply!(s, sCNOT(site_pair[1], site_pair[2]))
	else # X checks: system is target
	    apply!(s, sCNOT(site_pair[2], site_pair[1]))
	    apply!(s, sHadamard(site_pair[2]))
	end
    end
    s
end

# trace out the qubits at idxs[sites]
function trace_layer!(s, idxs, sites)
    if !isempty(sites)
        traceout!(s, idxs[sites])
    end
    s
end

# erase qubits at site idxs, and remove from system
# there's probably a faster way, but this is foolproof
function remove_qubits!(md::MixedDestabilizer, idxs)
    traceout!(md, idxs)
    if md.rank == 0
        return one(MixedDestabilizer,0, nqubits(md)-length(idxs))
    end
    md = MixedDestabilizer(stabilizerview(md)[:,setdiff(1:nqubits(md), idxs)])
end

#= Tracking operator spreading from root or fresh leaf =#
"""
get the stabilizer tableaus up to depth t, for t=1:tmax
Note: this function and its dependencies are taken from a larger project so are overly general.
For this experiment, it is always (H otimes H)CNOT, so the set of stabilizers to be measured could be precomputed and just stored as a global variable.
"""
function track_substabilizers(gate, tmax; init_pauli = P"Z")
    stab = Stabilizer(track_stabilizer_generators(gate, tmax; init_pauli = init_pauli))
    starts = vcat([0], cumsum([2^(tmax-t) for t=1:tmax-1]))
    vcat([stab[vcat([(1:2^(t1-t)) .+ starts[t] for t=1:t1]...),1:2^t1] for t1=1:tmax-1], [stab])
end

function track_operator_spreading(cliff, pauli, tmax; fresh::Bool = true)
    tabs = Array{Array}(undef, tmax)
    if fresh
        # feed in fresh, on right
	tabs[1] = Bool.(mod.(cliff * [0,pauli[1],0,pauli[2]], 2))
    else
        # feed in on "logical leg"
	tabs[1] = Bool.(mod.(cliff * [pauli[1],0,pauli[2],0], 2))
    end

    image = [tabs[1][[1,3]], tabs[1][[2,4]]]
    for t=2:tmax
        image = vcat([inflate_pauli(cliff, pauli) for pauli in image]...)
	bits = vcat(image...)
	tabs[t] = vcat(bits[1:2:end], bits[2:2:end])
    end

    tabs
end

# track the operators fed in on the right, on fresh branches, in the given levels
function track_fresh_stabilizers(gate, tmax, levels; init_pauli = P"Z")
    paulis = track_operator_spreading(stab_to_gf2(gate.tab)', Bool.(init_pauli.xz), tmax; fresh = true)

    # only get the stabilizers associated with the levels t
    stabs = [zeros(Bool, 2^(t-1), 2^(tmax+1)) for t=levels]
    for (t_i, t)=enumerate(levels)
        tp = tmax -t + 1
	for i=1:2^(t-1)
	    stabs[t_i][i, ((i-1)*2^tp+1):i*2^tp] = paulis[tp][1:end÷2]
	    stabs[t_i][i, (end÷2+(i-1)*2^tp+1):(end÷2+i*2^tp)] = paulis[tp][end÷2+1:end]
	end
    end
    stabs
end

# get the parity check matrix for the tree code with one logical
function track_stabilizer_generators(gate, tmax; init_pauli = P"Z")
    stabs = track_fresh_stabilizers(gate, tmax, 1:tmax; init_pauli = init_pauli)
    vcat(stabs[end:-1:1]...)
end

# apply inflation rule to input pauli, with I on right qubit
function inflate_pauli(cliff, pauli)
    inflated = Bool.(mod.(cliff * [pauli[1],0,pauli[2],0], 2))
    [inflated[[1,3]], inflated[[2,4]]]
end


end # module