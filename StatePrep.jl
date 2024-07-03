#!/usr/licensed/julia/1.10/bin/julia
module StatePrep

using QuantumClifford
using CircuitTools

export steane_erasures_blocks, prep_state_blocks, block_layer, pair_blocks

function sample_sites(r, n)
    return findall(el->el<=r, rand(n))
end

"""
rates[1]: rate of erasures on fresh qubits
rates[2]: rate of erasures after two-qubit gates (encoding gate and check gates)
rates[3]: rate of erasures before measurements

input_sites: locations of erasures on fresh qubits
encoding_sites: locations of erasures after encoding gate
check_sites: (1) locations of erasures that will be used to simulate possible wirings
	     (2) actual locations of erasures after check gates
measure_sites: (1) locations of erasures that will be used to simulate possible wirings
	       (2) actual locations of erasures after before measurements
"""
function get_erasure_blocks(rates, tmax)
    encoding_sites = [[sample_sites(rates[2], 2^t) for i=1:2^(tmax-t)] for t=1:tmax]

    check_sites = [[[sample_sites(rates[2], 2^t) for i=1:2^(tmax-t)] for t=2:tmax] for j=1:2]
    measure_sites = [[[sample_sites(rates[3], 2^t) for i=1:2^(tmax-t-1)] for t=1:tmax-1] for j=1:2]
    input_sites = vcat([[sample_sites(rates[1], 2) for i=1:2^(tmax-1)]], [[sample_sites(rates[1], 2^t) for i=1:2^(tmax-t-1)] for t=1:tmax-1])
    [input_sites, encoding_sites, check_sites, measure_sites]
end

function get_gate_idxs(t, tmax)
    return reshape([(j+i) .+ [0, 2^(t-1)] for i=1:2^(t-1), j=0:2^t:2^tmax-1],2^(tmax-1))
end

function get_system_idxs(tmax)
    system_idxs = [1]
    for t=1:tmax
        system_idxs = vcat(2 .* system_idxs .- 1, 2 .* system_idxs)
    end
    system_idxs
end

"""
Simulate possible wiring between "system" (control) and "ancilla" (control%2+1)
If check gate is performed, simulation includes sampled erasures (which may not be the actual erasures)

Return: 
    - `s`: updated stabilizer tableau of system (before perfect stabilizer measurements
    - `ent`: entropy of system following perfect stabilizer measurements
"""
function pair_blocks(stabs, erasure_sites, to_measure; check_gate::Bool = true, input =:Z, control::Int = 1, measure_type::Int = 1)

    sites = [[i,nqubits(stabs[1])+i] for i=1:nqubits(stabs[1])]
    if !check_gate # just discard one of the two states
        s = stabs[control]⊗MixedDestabilizer(one(Stabilizer,nqubits(stabs[control]), basis =input))
    else # use "control" state as system, measure target
    	s = stabs[control]⊗stabs[control%2+1]

    	transversal_cnot!(s, sites, measure_type)

    	# errors associated with check gates
    	trace_layer!(s, 1:nqubits(s), erasure_sites[1])

    	# errors associated with measurement
    	trace_layer!(s, nqubits(s)÷2+1:nqubits(s), erasure_sites[2])

    	# measure
    	for i=nqubits(s)÷2+1:nqubits(s)
            projectZ!(s, i; keep_result = false)
        end
    	if input !=:Z
            reset_qubits!(s, one(Stabilizer,nqubits(s)÷2, basis=input), nqubits(s)÷2+1:nqubits(s))
        end
    end
    
    # get entropy of state, following perfect stabilizer measurements
    ent = system_entropy!(copy(s), to_measure,nqubits(s))

    # reorder qubits
    # s.tab = s.tab[:,vcat(sites...)]
    s, ent
end

"""
Prepare logical state with encoding gate "gate", up to depth tmax
"""
function prep_state_blocks(gate, tmax, rates, to_measure; log_state = :Z, erasure_sites = nothing)
    my_s = MixedDestabilizer(one(Stabilizer,1,basis=log_state)⊗one(Stabilizer,1,basis=:Z))

    # erasure locations of each type
    if isnothing(erasure_sites)
        erasure_sites = get_erasure_blocks(rates, tmax)
    end

    if log_state==:Z # measure in X basis, Z basis, X basis...
        measure_first = 2
    else # log_state=:X
        measure_first = 1
    end
    
    big_gate = gate
    stabs = [deepcopy(my_s) for i=1:2^(tmax-1)]
    # up to depth tmax...
    for t=1:tmax-1
        stabs = block_layer(big_gate, t, to_measure[t], erasure_sites, stabs; parity = (t+measure_first)%2+1)
	big_gate = big_gate⊗big_gate
    end

    final_ent = system_entropy!(stabs[end], to_measure[end], 2^tmax)
end

# one layer of optimizing block pairs
function block_layer(gate, t, to_measure, erasure_sites, stabs; parity = 1)
    # up to depth tmax...
    gate_idxs = get_gate_idxs(t,t)
    if t==1
        input_sites = 1:2
    else
        input_sites = [idx[2] for idx in gate_idxs]
    end
    # Step 0: erasure errors on input legs (not actually used in Quantinuum experiment)
    for i=1:length(stabs)
	trace_layer!(stabs[i], input_sites, erasure_sites[1][t][i])
    end
	
    # Step 1: encoding gates
    for i=1:length(stabs)
        apply!(stabs[i], gate, vcat(gate_idxs...))
	trace_layer!(stabs[i], 1:2^t, erasure_sites[2][t][i])
    end

    # Step 2: CNOTs with ancillas
    gate_idxs = get_gate_idxs(t+1,t+1)
    input_sites = [idx[2] for idx in gate_idxs]

    # try four different ways of pairing up
    new_stabs = Array{MixedDestabilizer}(undef, length(stabs)÷2)
    entropies = zeros(Int,6, length(stabs)÷2)
    for i=1:length(stabs)÷2
    	ents = zeros(Int,4)
        j=1
	while j<=4
	    new_stabs[i], ents[j] = pair_blocks(stabs[2*i-1:2*i], [erasure_sites[3][1][t][i], erasure_sites[4][1][t][i]], to_measure; check_gate = (j<=2), measure_type = parity, control = (j+1)%2+1)
	    if entropies[j]==0
	        break
	    end
	    j+=1
	end

	if all(isone, ents) # if all wirings result in entropy 1, just default to first wiring
	    j=1
	end

	# if check gates are performed, evolve forward in time with the actual erasure pattern
	if j<=2
	    new_stabs[i],_ = pair_blocks(stabs[2*i-1:2*i], [erasure_sites[3][2][t][i], erasure_sites[4][2][t][i]], to_measure; check_gate = true, measure_type = parity, control = (j+1)%2+1)
	end
    end

    new_stabs
end


function system_entropy!(s, to_measure, L::Int)
    system_entropy!(s, to_measure, [1:L;])
end

# do a perfect measurement of all the stabilizers (not the logical)
function system_entropy!(s, to_measure, system_qubits::Array)
    for p in to_measure
        project!(s, p; keep_result = false)
    end
    # is my logical state pure?
    ent = entanglement_entropy(s, system_qubits, Val(:rref))
    @assert ent <= 1
    ent
end

"""
Main function to do classical simulation, at many different erasure rates
(When erasures are the only error source)
Parameters:
   - `gate`: 2-qubit gate applied at each node of the binary tree. For this experiment, it is always (H otimes H)CNOT, so the set of stabilizers to be measured could be precomputed.
   - `tmax`: total depth of logical state to be prepared
   - `rates`: rate of erasures, (1) on inputs, (2) after two-qubit gates, (3) before measurements
   - `log_state`: logical state to prepare. Either :Z or :X
   - `num_samples`: number of samples to take at each erasure rate
Return:
   - `entropies`: num_samples x length(rates) array of final entropies
"""
function steane_erasures_blocks(gate, tmax, rates; log_state =:Z, num_samples::Int = 100)

    # stabilizer generators that get measured in each layer, when determining entropy
    to_measure = [s[:,get_system_idxs(t)] for (t,s)=enumerate(track_substabilizers(gate, tmax,init_pauli=P"Z"))]
    for t=1:tmax-1
        id_string = PauliOperator(zeros(Bool, 2*nqubits(to_measure[t])))
	to_measure[t] = Stabilizer([p⊗id_string for p in to_measure[t]])
    end
    
    entropies = zeros(Int, num_samples, length(rates))
    for r_i=1:length(rates)
        for num_i=1:num_samples
	    entropies[num_i,r_i] = prep_state_blocks(gate, tmax, rates[r_i], to_measure; log_state = log_state)
	end
    end
    entropies
end

end # module