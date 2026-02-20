# Quantum Emulator Specification

## Individual Logic Tests

This file defines every capability a quantum emulator should have, expressed as individual testable properties. Each test is self-contained and verifiable.

---

# 1. FOUNDATIONAL STATE TESTS

## 1.1 Basis States
- [ ] `test_init_zero_state`: Initialize qubits to |0⟩ produces correct state vector
- [ ] `test_init_one_state`: Initialize qubits to |1⟩ produces correct state vector
- [ ] `test_init superposition`: Initialize to arbitrary superposition
- [ ] `test_init_custom_state`: Initialize to user-defined complex amplitudes

## 1.2 State Vector Properties
- [ ] `test_state_norm_one`: State vector is always normalized (norm = 1)
- [ ] `test_state_probabilities_sum_one`: All measurement probabilities sum to 1
- [ ] `test_state_dimensions`: 2^n dimensions for n qubits
- [ ] `test_partial_trace_dimensions`: Partial trace produces valid reduced density matrix

## 1.3 State Validity
- [ ] `test_state_positive`: Density matrix is positive semidefinite
- [ ] `test_state_trace_one`: Density matrix has trace = 1
- [ ] `test_state_hermitian`: Density matrix is Hermitian

---

# 2. SINGLE-QUBIT GATE TESTS

## 2.1 Pauli Gates
- [ ] `test_pauli_x_on_zero`: X|0⟩ = |1⟩
- [ ] `test_pauli_x_on_one`: X|1⟩ = |0⟩
- [ ] `test_pauli_x_on_superposition`: X applied to superposition
- [ ] `test_pauli_y_on_zero`: Y|0⟩ = i|1⟩
- [ ] `test_pauli_y_on_one`: Y|1⟩ = -i|0⟩
- [ ] `test_pauli_z_on_zero`: Z|0⟩ = |0⟩
- [ ] `test_pauli_z_on_one`: Z|1⟩ = -|1⟩
- [ ] `test_pauli_z_preserves_superposition_phase`: Z preserves relative phase
- [ ] `test_pauli_x_squared_identity`: X² = I
- [ ] `test_pauli_y_squared_identity`: Y² = I
- [ ] `test_pauli_z_squared_identity`: Z² = I
- [ ] `test_pauli_anticommutation`: {X, Y} = 0 (anticommute)
- [ ] `test_pauli_x_y_commutation`: [X, Y] = 2iZ

## 2.2 Hadamard Gate
- [ ] `test_hadamard_on_zero`: H|0⟩ = (|0⟩ + |1⟩)/√2
- [ ] `test_hadamard_on_one`: H|1⟩ = (|0⟩ - |1⟩)/√2
- [ ] `test_hadamard_superposition`: H on superposition maintains proper amplitudes
- [ ] `test_hadamard_squared_identity`: H² = I
- [ ] `test_hadamard_self_inverse`: H = H⁻¹
- [ ] `test_hadamard_preserves_norm`: H preserves vector norm

## 2.3 Phase Gates
- [ ] `test_phase_gate_S_on_zero`: S|0⟩ = |0⟩
- [ ] `test_phase_gate_S_on_one`: S|1⟩ = i|1⟩
- [ ] `test_phase_gate_T_on_zero`: T|0⟩ = |0⟩
- [ ] `test_phase_gate_T_on_one`: T|1⟩ = e^(iπ/4)|1⟩
- [ ] `test_phase_gate_S_squared`: S² = Z
- [ ] `test_phase_gate_T_fourth_power`: T⁴ = Z
- [ ] `test_phase_gate_commutation`: S and T commute (both diagonal)

## 2.4 Rotation Gates
- [ ] `test_rotation_X_gate`: Rx(θ) = cos(θ/2)I - i sin(θ/2)X
- [ ] `test_rotation_Y_gate`: Ry(θ) = cos(θ/2)I - i sin(θ/2)Y
- [ ] `test_rotation_Z_gate`: Rz(θ) = e^(-iθ/2)|0⟩⟨0| + e^(iθ/2)|1⟩⟨1|
- [ ] `test_rotation_X_zero`: Rx(0) = I
- [ ] `test_rotation_Y_pi_half`: Ry(π/2) creates expected superposition
- [ ] `test_rotation_Z_phases`: Rz correctly applies phase to |1⟩
- [ ] `test_rotation_X_pi`: Rx(π) = -iX
- [ ] `test_rotation_continuous`: Rotations form continuous group

## 2.5 General Single-Qubit Gates
- [ ] `test_u3_gate_parameterization`: U3(θ, φ, λ) covers all single-qubit unitaries
- [ ] `test_identity_preserves_state`: I|ψ⟩ = |ψ⟩
- [ ] `test_gate_unitarity`: All gates are unitary (U†U = I)
- [ ] `test_gate_determinant_one`: All gates have determinant 1

---

# 3. MULTI-QUBIT GATE TESTS

## 3.1 CNOT Gate
- [ ] `test_cnot_00_to_00`: CNOT|00⟩ = |00⟩
- [ ] `test_cnot_01_to_01`: CNOT|01⟩ = |01⟩
- [ ] `test_cnot_10_to_11`: CNOT|10⟩ = |11⟩
- [ ] `test_cnot_11_to_10`: CNOT|11⟩ = |10⟩
- [ ] `test_cnot_control_flips_target`: Control |1⟩ flips target
- [ ] `test_cnot_control_preserves_target`: Control |0⟩ preserves target

## 3.2 Controlled Gates (General)
- [ ] `test_controlled_gate_cz`: CZ|++⟩ produces expected sign flip
- [ ] `test_controlled_gate_cs`: CS applies phase to |11⟩
- [ ] `test_controlled_gate_cnot`: Controlled-X is CNOT
- [ ] `test_controlled_gate_cswap`: Fredkin gate (CSWAP)
- [ ] `test_controlled_gate_reversed`: Control on second qubit

## 3.3 Toffoli Gate (CCNOT)
- [ ] `test_toffoli_000_to_000`: CCNOT|000⟩ = |000⟩
- [ ] `test_toffoli_001_to_001`: CCNOT|001⟩ = |001⟩
- [ ] `test_toffoli_010_to_010`: CCNOT|010⟩ = |010⟩
- [ ] `test_toffoli_011_to_011`: CCNOT|011⟩ = |011⟩
- [ ] `test_toffoli_100_to_100`: CCNOT|100⟩ = |100⟩
- [ ] `test_toffoli_101_to_101`: CCNOT|101⟩ = |101⟩
- [ ] `test_toffoli_110_to_111`: CCNOT|110⟩ = |111⟩
- [ ] `test_toffoli_111_to_110`: CCNOT|111⟩ = |110⟩
- [ ] `test_toffoli_three_control`: Only flips when all three controls are |1⟩

## 3.4 SWAP and ISWAP
- [ ] `test_swap_gate`: SWAP|01⟩ = |10⟩, SWAP|10⟩ = |01⟩
- [ ] `test_swap_symmetric`: SWAP is its own inverse
- [ ] `test_iswap_gate`: ISWAP adds π/2 phase
- [ ] `test_iswap_plus_gate`: ISWAP† is ISWAP⁻¹
- [ ] `test_swap_three_qubits`: SWAP can reorder qubits

## 3.5 Fredkin Gate (CSWAP)
- [ ] `test_fredkin_000`: Fredkin|000⟩ = |000⟩
- [ ] `test_fredkin_001`: Fredkin|001⟩ = |001⟩
- [ ] `test_fredkin_010`: Fredkin|010⟩ = |010⟩
- [ ] `test_fredkin_011`: Fredkin|011⟩ = |110⟩ (swaps)
- [ ] `test_fredkin_100`: Fredkin|100⟩ = |100⟩
- [ ] `test_fredkin_101`: Fredkin|101⟩ = |101⟩
- [ ] `test_fredkin_110`: Fredkin|110⟩ = |011⟩ (swaps)
- [ ] `test_fredkin_111`: Fredkin|111⟩ = |111⟩

## 3.6 Gate Decompositions
- [ ] `test_cnot_from_hadamard_cz_hadamard`: CNOT = (I⊗H)CZ(I⊗H)
- [ ] `test_toffoli_decomposition`: Toffoli can be decomposed into single-qubit and CNOT
- [ ] `test_swap_from_three_cnots`: SWAP = three CNOTs
- [ ] `test_hadamard_from_rotations`: H can be decomposed into rotations

---

# 4. ENTANGLEMENT TESTS

## 4.1 Bell States
- [ ] `test_bell_state_phi_plus`: Create |Φ⁺⟩ = (|00⟩ + |11⟩)/√2
- [ ] `test_bell_state_phi_minus`: Create |Φ⁻⟩ = (|00⟩ - |11⟩)/√2
- [ ] `test_bell_state_psi_plus`: Create |Ψ⁺⟩ = (|01⟩ + |10⟩)/√2
- [ ] `test_bell_state_psi_minus`: Create |Ψ⁻⟩ = (|01⟩ - |10⟩)/√2
- [ ] `test_bell_states_orthogonal`: Bell states are mutually orthogonal
- [ ] `test_bell_measurement_correlation`: Measuring Bell state gives correlated results

## 4.2 GHZ and W States
- [ ] `test_ghz_state_n3`: Create |GHZ₃⟩ = (|000⟩ + |111⟩)/√2
- [ ] `test_ghz_state_n4`: Create |GHZ₄⟩ = (|0000⟩ + |1111⟩)/√2
- [ ] `test_w_state_n3`: Create |W₃⟩ = (|001⟩ + |010⟩ + |100⟩)/√3
- [ ] `test_ghz_classical_correlation`: GHZ shows all-or-nothing correlation

## 4.3 Entanglement Detection
- [ ] `test_concurrence_bell`: Bell states have concurrence = 1
- [ ] `test_concurrence_product`: Product states have concurrence = 0
- [ ] `test_entanglement_witness_bell`: Bell states violate appropriate witness
- [ ] `test_partial_transpose_criterion`: PPT criterion detects entanglement

## 4.4 Entanglement Manipulation
- [ ] `test_entanglement_preserved_under_local`: Local unitaries preserve entanglement
- [ ] `test_entanglement_cnot_creates`: CNOT on product state creates entanglement
- [ ] `test_entanglement_measuring_one_ qubit`: Measuring one qubit of Bell pair destroys entanglement

---

# 5. MEASUREMENT TESTS

## 5.1 Basic Measurement
- [ ] `test_measure_basis_state_zero`: Measuring |0⟩ gives 0 with probability 1
- [ ] `test_measure_basis_state_one`: Measuring |1⟩ gives 1 with probability 1
- [ ] `test_measure_superposition_random`: Superposition gives probabilistic results
- [ ] `test_measure_updates_state`: Measurement collapses state to outcome basis state
- [ ] `test_measure_mixed_basis`: Measuring in wrong basis gives incorrect results

## 5.2 Measurement Probabilities
- [ ] `test_measure_probability_computation`: Probability = |amplitude|²
- [ ] `test_measure_expectation_z`: ⟨Z⟩ = P(0) - P(1)
- [ ] `test_measure_statistics_converge`: Repeated measurements converge to theoretical distribution
- [ ] `test_measure_ Born rule`: Born rule holds for any state

## 5.3 Partial Measurement
- [ ] `test_partial_measure_first_qubit`: Measuring first qubit of pair collapses correctly
- [ ] `test_partial_measure_second_qubit`: Measuring second qubit works symmetrically
- [ ] `test_partial_trace_verification`: Reduced density matrix matches partial measurement

## 5.4 Measurement in Different Bases
- [ ] `test_measure_x_basis`: Measuring in X basis = Hadamard then Z basis
- [ ] `test_measure_y_basis`: Measuring in Y basis = S†H then Z basis
- [ ] `test_measure_random_basis`: Arbitrary basis measurement works
- [ ] `test_measure_bell_basis`: Bell basis measurement distinguishes entangled states

## 5.5 Repeated Measurement
- [ ] `test_repeated_measure_same`: Measuring same qubit twice gives same result
- [ ] `test_repeated_measure_doesnt_change`: Second measurement doesn't modify state

---

# 6. QUANTUM ALGORITHM TESTS

## 6.1 Deutsch-Jozsa Algorithm
- [ ] `test_deutsch_jozsa_constant_balanced`: Correctly identifies constant vs balanced
- [ ] `test_deutsch_jozsa_all_constant`: f(x) = 0 is constant
- [ ] `test_deutsch_jozsa_all_one`: f(x) = 1 is constant
- [ ] `test_deutsch_jozsa_parity`: f(x) = x₀ ⊕ x₁ is balanced
- [ ] `test_deutsch_jozsa_oracle_construction`: Oracle correctly implements function

## 6.2 Bernstein-Vazirani Algorithm
- [ ] `test_bernstein_vazirani_bit`: Finds hidden bit string
- [ ] `test_bernstein_vazirani_multi_bit`: Works for n-bit string
- [ ] `test_bernstein_vazirani_oracle`: Oracle correctly computes a·x

## 6.3 Simon's Algorithm
- [ ] `test_simon_algorithm_period`: Finds period of function
- [ ] `test_simon_algorithm_oracle`: Oracle satisfies f(x) = f(x ⊕ s)
- [ ] `test_simon_post_processing`: Classical post-processing finds s

## 6.4 Grover's Search Algorithm
- [ ] `test_grover_single_solution`: Finds single marked item with ~√N queries
- [ ] `test_grover_multiple_solutions`: Works with multiple marked items
- [ ] `test_grover_optimal_iterations`: Correct number of iterations for N items, M solutions
- [ ] `test_grover_amplitude_amplification`: Amplitude amplification works correctly
- [ ] `test_grover_oracle_construction`: Oracle correctly marks target state

## 6.5 Quantum Fourier Transform
- [ ] `test_qft_basis_states`: QFT|k⟩ produces correct output
- [ ] `test_qft_superposition`: QFT on superposition is correct
- [ ] `test_qft_inverse`: QFT⁻¹ = QFT†
- [ ] `test_qft_four_qubits`: QFT on 4 qubits produces correct phases
- [ ] `test_qft_phase_estimation`: QFT enables phase estimation

## 6.6 Quantum Phase Estimation
- [ ] `test_phase_estimation_eigenstate`: Estimates phase of eigenstate exactly
- [ ] `test_phase_estimation_superposition`: Handles superposition of eigenstates
- [ ] `test_phase_estimation_precision`: Higher precision with more bits
- [ ] `test_phase_estimation_oracle`: Oracle correctly applies controlled unitaries

## 6.7 Quantum Teleportation
- [ ] `test_teleportation_protocol`: Alice's qubits teleport to Bob
- [ ] `test_teleportation_entanglement`: Requires shared Bell pair
- [ ] `test_teleportation_classical_comm`: Needs 2 classical bits
- [ ] `test_teleportation_no_cloning`: Original state destroyed

## 6.8 Superdense Coding
- [ ] `test_superdense_two_bits`: Two classical bits from one qubit
- [ ] `test_superdense_encoding`: Correct encoding for each pair
- [ ] `test_superdense_bell_measurement`: Bell measurement at receiver

---

# 7. QUANTUM PROPERTIES THEOREMS

## 7.1 No-Cloning Theorem
- [ ] `test_no_cloning_impossible`: Cannot clone arbitrary quantum state
- [ ] `test_no_cloning_no_override`: Even with resources, cannot perfectly clone
- [ ] `test_no_teleportation_without_classical`: Teleportation needs classical channel

## 7.2 No-Signaling
- [ ] `test_no_signaling_faster_than_light`: Cannot use entanglement for FTL communication
- [ ] `test_measurement_outcome_independence`: Remote measurement doesn't affect local

## 7.3 Uncertainty Relations
- [ ] `test_heisenberg_uncertainty`: ΔX ΔP ≥ ℏ/2
- [ ] `test_information_disturbance`: Measuring disturbs conjugate observable
- [ ] `test_complementarity`: Cannot know both complementary quantities

## 7.4 Quantum Correlations
- [ ] `test_bell_inequality_violation`: CHSH > 2 for entangled state
- [ ] `test_classical_limit_chsh`: Classical correlations give CHSH ≤ 2

---

# 8. NOISE AND ERROR MODELS

## 8.1 Decoherence
- [ ] `test_depolarizing_channel`: Mixed state under depolarizing noise
- [ ] `test_amplitude_damping`: Energy loss in T1 relaxation
- [ ] `test_phase_damping`: Phase noise in T2 relaxation
- [ ] `test_thermal_state`: Thermal equilibrium state

## 8.2 Gate Errors
- [ ] `test_depolarizing_gate`: Gate with depolarizing error
- [ ] `test_amplitude_error`: Gate with rotation error
- [ ] `test_phase_error`: Gate with phase error
- [ ] `test_gate_fidelity`: Average gate fidelity measurement

## 8.3 Error Correction
- [ ] `test_bit_flip_code`: 3-qubit bit-flip code corrects single error
- [ ] `test_phase_flip_code`: 3-qubit phase-flip code
- [ ] `test_shor_code`: Shor code corrects both
- [ ] `test_syndrome_extraction`: Syndrome measurement identifies error

---

# 9. QUANTUM METRICS

## 9.1 Fidelity
- [ ] `test_fidelity_pure_same`: F(|ψ⟩, |ψ⟩) = 1
- [ ] `test_fidelity_pure_orthogonal`: F(|ψ⟩, |φ⟩) = 0 for orthogonal
- [ ] `test_fidelity_mixed_state`: Fidelity works with density matrices
- [ ] `test_fidelity_triangle_inequality`: F satisfies triangle inequality

## 9.2 Entropy
- [ ] `test_von_neumann_pure`: S(|ψ⟩⟨ψ|) = 0 for pure states
- [ ] `test_von_neumann_maximal`: S = log(d) for maximally mixed
- [ ] `test_entanglement_entropy_bipartite`: Reduced entropy measures entanglement

## 9.3 Distance Measures
- [ ] `test_trace_distance`: Trace distance is a metric
- [ ] `test_trace_distance_max`: Max distance between pure states = 1
- [ ] `test_fidelity_trace_relation`: F = 1 - D for pure states

---

# 10. CIRCUIT PROPERTIES

## 10.1 Circuit Structure
- [ ] `test_circuit_depth`: Depth calculation is correct
- [ ] `test_circuit_width`: Width = number of qubits
- [ ] `test_circuit_gate_count`: Gate count matches added gates
- [ ] `test_circuit_commutation`: Commuting gates can be reordered

## 10.2 Circuit Composition
- [ ] `test_circuit_compose`: Two circuits can be composed
- [ ] `test_circuit_inverse`: Circuit has correct inverse
- [ ] `test_circuit_identity`: Empty circuit = identity
- [ ] `test_circuit_parallel`: Parallel gates on different qubits

## 10.3 Circuit Transformations
- [ ] `test_circuit_draw`: Circuit renders to string/diagram
- [ ] `test_circuit_to_matrix`: Circuit compiles to unitary matrix
- [ ] `test_circuit_simplify`: Optimization pass simplifies circuit
- [ ] `test_circuit_decompose`: Decompose to basis gates

---

# 11. ADVANCED ALGORITHMS

## 11.1 Variational Algorithms
- [ ] `test_vqe_hamiltonian_expectation`: Calculates expectation correctly
- [ ] `test_vqe_parameter_optimization`: Gradient descent finds ground state
- [ ] `test_ansatz_variational`: Ansatz is variational (adjustable parameters)

## 11.2 Quantum Approximate Optimization
- [ ] `test_qaoa_cost_function`: QAOA evaluates cost function
- [ ] `test_qaoa_mixer_layer`: Mixer creates superposition
- [ ] `test_qaoa_parametrized_layers`: Layers are parameterized

## 11.3 Hamiltonian Simulation
- [ ] `test_trotter_decomposition`: First-order Trotter is accurate
- [ ] `test_trotter_error_bounds`: Error decreases with more steps
- [ ] `test_hamiltonian_simulation_exact`: Exact simulation for diagonal H

---

# 12. SPECIAL STATES AND TRANSFORMATIONS

## 12.1 Coherent States
- [ ] `test_coherent_state_gaussian`: Coherent states are Gaussian
- [ ] `test_coherent_state_displacement`: Displacement operator creates coherent state

## 12.2 Dicke States
- [ ] `test_dicke_state_k_excitations`: Dicke state |Dⁿₖ⟩ is symmetric
- [ ] `test_dicke_state_preparation`: Can prepare Dicke states

## 12.3 Cluster States
- [ ] `test_cluster_state_graph`: Cluster state from graph
- [ ] `test_cluster_measurements`: Measurements enable MBQC

---

# 13. NUMERICAL PRECISION

## 13.1 Complex Arithmetic
- [ ] `test_complex_normalization`: Complex numbers normalize correctly
- [ ] `test_complex_multiplication`: Complex multiplication is accurate
- [ ] `test_complex_conjugate_product`: |z|² = z*z

## 13.2 Floating Point
- [ ] `test_precision_small_amplitudes`: Handles very small amplitudes
- [ ] `test_precision_phase_preservation`: Preserves small phase differences
- [ ] `test_precision_numerical_stability`: No catastrophic cancellation

---

# 14. EDGE CASES

## 14.1 Boundary Conditions
- [ ] `test_zero_qubits`: Zero-qubit circuit (trivial)
- [ ] `test_single_qubit_circuit`: Single qubit works correctly
- [ ] `test_large_superposition`: Handles 2^30 dimensional state
- [ ] `test_identity_circuit`: All identity gates

## 14.2 Special Gates
- [ ] `test_global_phase`: Global phase unobservable
- [ ] `test_determinant_calculation`: Gate determinant always 1
- [ ] `test_adjoint_gate`: Adjoint computation is correct

## 14.3 Error Handling
- [ ] `test_invalid_qubit_index`: Rejects invalid qubit indices
- [ ] `test_negative_angles`: Handles negative rotation angles
- [ ] `test_empty_circuit`: Empty circuit is valid

---

# 15. INTEGRATION TESTS

## 15.1 Full Algorithm Verification
- [ ] `test_deutsch_jozsa_full`: Complete Deutsch-Jozsa from oracle to result
- [ ] `test_grover_full`: Complete Grover from oracle to solution
- [ ] `test_qft_full`: Complete QFT from input to output
- [ ] `test_teleportation_full`: Complete teleportation protocol

## 15.2 Comparison with Known Results
- [ ] `test_bell_inequality_value`: CHSH value for Bell state = 2√2
- [ ] `test_ghz_correlation_value`: GHZ correlation matches theory
- [ ] `test_shor_factorization`: Shor's factors small numbers

---

# 16. PERFORMANCE TARGETS

## 16.1 Simulation Limits
- [ ] `test_simulate_10_qubits`: Simulates 10 qubits (1024 states) efficiently
- [ ] `test_simulate_15_qubits`: Simulates 15 qubits (32768 states) reasonably
- [ ] `test_simulate_20_qubits`: 20 qubits achievable with optimizations
- [ ] `test_sparse_vs_dense`: Sparse representation faster for certain circuits

## 16.2 Gate Application Speed
- [ ] `test_single_qubit_gate_speed`: Single-qubit gates are fast
- [ ] `test_cnot_gate_speed`: CNOT is not much slower than single-qubit
- [ ] `test_parallel_gate_speed`: Multiple gates in sequence is efficient
