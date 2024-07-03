using StatePrep
using QuantumClifford

# Example of how I would do this experiment fully classically
function run_classical()
    gate = (tHadamard⊗tHadamard)*tCNOT
    tmax = 5
    rates = [[0,0.05,0.05],[0,0.1,0.1],[0,0.15,0.15]]
    steane_erasures_blocks(gate,tmax,rates; log_state=:X)
end