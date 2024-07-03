# state-prep

To use the modules included in this repo, execute the following from the command line:

```export JULIA_LOAD_PATH="~/path/to/repo:"```

or, inside a Julia session, do:

```push!(LOAD_PATH, "~/path/to/repo")```

Then you can access the functions in these modules through the command `using <module name>`.

The code we'd like to call using hybrid-compute is lines 146-157 in StatePrep.jl
Would then continue to simulate classically (the rest of the function `block_layer`), while running the chosen wiring on the quantum computer. 
Required packages
-----------------
  - `QuantumClifford` (for stabilizer formalism)