# Simulator and Testcase Guidelines

This file is the longer operational guide for running the simulator, generating testcases, and producing exact statevector ground truth for small instances.

## 1. Main Simulator

Build from `algorithm`:

```bash
cd algorithm
make
```

Run:

```bash
./code
```

Default run:

```text
Default setting? (y/n): y
test file name: QKE_n2
```

Custom run:

```text
Default setting? (y/n): n
Thread: 24
Fix n_shot ? (y/n): n
Delta: 0.2
Epsilon: 0.2
Use Empirical Variance? (y/n): y
Iteration: 1
test file name: QKE_n2
```

If `Fix n_shot ?` is `y`, the simulator uses the provided fixed shot count.

The simulator reads:

```text
../testcases/<name>/info.txt
../testcases/<name>/observable.txt
../testcases/<name>/circuits/*.txt
```

It writes:

```text
../data/details_<name>_<timestamp>.csv
```

## 2. Generate A Random Pauli Feature-Map QKE Testcase

`testcase_generator/RandomPauliFeatureMap.cpp` generates QKE testcases using:

$$
U(x) = \exp\left(i \sum_i x_i Z_i
           + i \sum_j (\pi - x_j)(\pi - x_{j+1}) Z_j Z_{j+1}\right)
$$

It samples:

$$
x_i \sim \text{Uniform}(0, 2\pi)
$$

and writes circuits for:

$$
V(x, x') = H^{\otimes n} U^\dagger(x) U(x') H^{\otimes n}
$$

The observable is:

```text
projector 0
```

which means $\ket{0^{\otimes n}}\bra{0^{\otimes n}}$.

Build:

```bash
cd testcase_generator
make
```

Usage:

```bash
./qke_gen <N_random_vectors> <n_qubits> [seed] [output_name] [pair_mode]
```

Example:

```bash
./qke_gen 10 64 5489 QKE_n64 full
```

Arguments:

```text
N_random_vectors   Number of random samples
n_qubits           Number of qubits
seed               Optional. Default: 5489
output_name        Optional. Default: QKE_RandomPauli_N<N>_n<n>.
pair_mode          Optional. full, ordered, or unordered. Default: full
```

Pair modes:

```text
full       Generates N^2 circuits: (i,j) for all i,j
ordered    Generates N(N-1) circuits: all i != j
unordered  Generates N(N-1)/2 circuits: only i < j
```

For `N=10` and `full`, pair order is row-major:

```text
(0,0), (0,1), ..., (0,9),
(1,0), (1,1), ..., (1,9),
...
(9,0), (9,1), ..., (9,9)
```

The generator writes:

```text
info.txt
observable.txt
metadata.txt
samples.csv
pairs.csv
circuits/*.txt
```

When run from any folder in repository, this writes `../testcases/QKE_n64/`. Run the generated testcase:

```bash
cd ../algorithm
printf "y\nQKE_n64\n" | ./code
```

Notes:

- The seed makes the random vectors reproducible.
- Change the seed for a new random dataset.
- `(i,j)` and `(j,i)` are different circuit files, but their exact kernel
  values should agree because the kernel is an overlap magnitude squared.
- `(i,i)` should evaluate to 1 in exact arithmetic.

## 3. Generate An E3LIN2 Testcase

`testcase_generator/E3LIN2.cpp` generates MaxE3LIN ansatz circuits for benchmarks.

Build from `testcase_generator`:

```bash
cd testcase_generator
make
```

Usage:

```bash
./qaoa_gen <N> <D> <Ngamma> [seed] [output_name]
```

Example:

```bash
./qaoa_gen 60 4 31 77 N60D4_custom
```

Arguments:

```text
N                  Number of qubits
D                  degree
Ngamma             Number of gammas
seed               Optional. Default: 19
output_name        Optional. Default:N<N>D<D>
```

When run from any folder in repository, a bare output name writes under `../testcases/`.
For example, the default command above creates:

```text
testcases/N60D4/
```

Then run:

```bash
cd ../algorithm
printf "y\nN60D4\n" | ./code
```

Notes:

- The seed makes the random vectors reproducible.
- `# of gammas` must be at least 2.
- The generator creates `N * D / 3` random `ZZZ` clauses.
- The observable is a Pauli sum with one `ZZZ` term per clause.

## 4. Statevector Ground Truth

`testcases/statevector_ground_truth.py` gives exact dense statevector results
for small testcases.

Default usage:

```bash
cd testcases
python3 statevector_ground_truth.py QKE_n5
```

Default output:

```text
testcases/QKE_n5/ground_truth_statevector.csv
```

Custom output:

```bash
python3 testcases/statevector_ground_truth.py testcases/QKE_n5 \
  --output /tmp/QKE_n5_ground_truth.csv
```

Print to terminal:

```bash
python3 testcases/statevector_ground_truth.py testcases/QKE_n5 --stdout
```

Only evaluate the first few circuits:

```bash
python3 testcases/statevector_ground_truth.py testcases/QKE_n5 --limit 5
```

The default limit is:

```text
--max-qubits 20
```

## 5. Writing Your Own Testcases

Each testcase must have:

```text
testcases/<name>/
  info.txt
  observable.txt
  circuits/
    0.txt
    1.txt
    ...
```

### info.txt

The first integer is the number of circuits:

```text
3
```

The simulator then reads:

```text
circuits/0.txt
circuits/1.txt
circuits/2.txt
```

### Circuit Files

Each circuit starts with:

```text
<n_qubits> <n_gates>
```

then one gate per line:

```text
2 4
h 0
u1 0 0.3
cx 0 1
h 1
```

Qubit indices are zero-based.

Support gate lines:

```text
x q
y q
z q
h q
s q
sdg q
sx q
sxdg q
t q
tdg q
u1 q theta
rz q theta
cx control target
cz control target
swap q0 q1
id q
delay q
u0 q
pauli STRING k q0 q1 ... q{k-1}
```

### Pauli-Sum Observable

The observable may optionally begin with one or more metadata lines:

```text
norm <value>
operator_norm <value>
```

If the norm metadata is omitted, the simulator computes the exact norm when the observable terms allow it and otherwise falls back to the Frobenius norm lower bound $\|\boldsymbol{c}\|_2$.

The observable body starts with the number of terms:

```text
2
1 pauli Z 1 0
-1 pauli ZZ 2 0 1
```

Each term has the form:

```text
<coefficient> pauli <pauli_string> <number_of_qubits> <qubit_0> ...
```

### Rank-1 Projector Observable

Projector onto $\ket{0^{\otimes n}}$:

```text
projector 0
```

Projector onto a Clifford-prepared stabilizer state:

```text
projector 2
h 0
cx 0 1
```

This means:

$$
\ket{\phi}\bra{\phi}, \text{ where } \ket{\phi} = CX_{0,1} H_0 \ket{0^{\otimes n}}
$$

The projector preparation circuit must be Clifford-only.