"""
Copyright 2026 Yunseo Hwang and Giwon Song.

SPDX-License-Identifier: Apache-2.0
"""

#!/usr/bin/env python3
"""Exact dense statevector evaluator for this repository's testcase format.

This is intended for small testcases only. It reads:
  <testcase>/info.txt
  <testcase>/observable.txt
  <testcase>/circuits/<id>.txt
and writes exact expectation values for each circuit.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
from typing import Iterable

import numpy as np


I2 = np.eye(2, dtype=np.complex128)
X = np.array([[0, 1], [1, 0]], dtype=np.complex128)
Y = np.array([[0, -1j], [1j, 0]], dtype=np.complex128)
Z = np.array([[1, 0], [0, -1]], dtype=np.complex128)
H = (1.0 / math.sqrt(2.0)) * np.array([[1, 1], [1, -1]], dtype=np.complex128)
S = np.array([[1, 0], [0, 1j]], dtype=np.complex128)
SDG = np.array([[1, 0], [0, -1j]], dtype=np.complex128)
T = np.array([[1, 0], [0, np.exp(1j * math.pi / 4.0)]], dtype=np.complex128)
TDG = np.array([[1, 0], [0, np.exp(-1j * math.pi / 4.0)]], dtype=np.complex128)
# Match the repository simulator convention: sx is applied as H S H.
SX = H @ S @ H
SXDG = H @ SDG @ H


def parse_complex_token(token: str) -> complex:
    token = token.strip()
    if token.startswith("(") and token.endswith(")"):
        body = token[1:-1]
        if "," in body:
            real, imag = body.split(",", 1)
            return complex(float(real), float(imag))
        return complex(body.replace("i", "j"))
    return complex(float(token), 0.0)


def apply_1q(state: np.ndarray, n: int, qubit: int, mat: np.ndarray) -> None:
    bit = 1 << qubit
    for i in range(1 << n):
        if i & bit:
            continue
        j = i | bit
        a = state[i]
        b = state[j]
        state[i] = mat[0, 0] * a + mat[0, 1] * b
        state[j] = mat[1, 0] * a + mat[1, 1] * b


def apply_cx(state: np.ndarray, n: int, control: int, target: int) -> None:
    cbit = 1 << control
    tbit = 1 << target
    for i in range(1 << n):
        if (i & cbit) and not (i & tbit):
            j = i | tbit
            state[i], state[j] = state[j], state[i]


def apply_cz(state: np.ndarray, n: int, control: int, target: int) -> None:
    mask = (1 << control) | (1 << target)
    for i in range(1 << n):
        if (i & mask) == mask:
            state[i] *= -1


def apply_swap(state: np.ndarray, n: int, q1: int, q2: int) -> None:
    if q1 == q2:
        return
    b1 = 1 << q1
    b2 = 1 << q2
    for i in range(1 << n):
        if bool(i & b1) == bool(i & b2):
            continue
        j = i ^ b1 ^ b2
        if i < j:
            state[i], state[j] = state[j], state[i]


def apply_ccx(state: np.ndarray, n: int, c1: int, c2: int, target: int) -> None:
    c1bit = 1 << c1
    c2bit = 1 << c2
    tbit = 1 << target
    for i in range(1 << n):
        if (i & c1bit) and (i & c2bit) and not (i & tbit):
            j = i | tbit
            state[i], state[j] = state[j], state[i]


def apply_ccz(state: np.ndarray, n: int, c1: int, c2: int, target: int) -> None:
    mask = (1 << c1) | (1 << c2) | (1 << target)
    for i in range(1 << n):
        if (i & mask) == mask:
            state[i] *= -1


def apply_pauli_string(state: np.ndarray, n: int, paulis: str, qubits: list[int]) -> np.ndarray:
    out = state.copy()
    for p, q in zip(paulis, qubits):
        if p == "I":
            continue
        if p == "X":
            apply_1q(out, n, q, X)
        elif p == "Y":
            apply_1q(out, n, q, Y)
        elif p == "Z":
            apply_1q(out, n, q, Z)
        elif p == "H":
            # Some observable tests use H as a Clifford observable atom.
            apply_1q(out, n, q, H)
        else:
            raise ValueError(f"unsupported Pauli/Clifford observable atom {p!r}")
    return out


def apply_gate(state: np.ndarray, n: int, parts: list[str]) -> None:
    gate = parts[0]
    if gate == "x":
        apply_1q(state, n, int(parts[1]), X)
    elif gate == "y":
        apply_1q(state, n, int(parts[1]), Y)
    elif gate == "z":
        apply_1q(state, n, int(parts[1]), Z)
    elif gate == "h":
        apply_1q(state, n, int(parts[1]), H)
    elif gate == "s":
        apply_1q(state, n, int(parts[1]), S)
    elif gate == "sdg":
        apply_1q(state, n, int(parts[1]), SDG)
    elif gate == "sx":
        apply_1q(state, n, int(parts[1]), SX)
    elif gate == "sxdg":
        apply_1q(state, n, int(parts[1]), SXDG)
    elif gate == "t":
        apply_1q(state, n, int(parts[1]), T)
    elif gate == "tdg":
        apply_1q(state, n, int(parts[1]), TDG)
    elif gate in {"u1", "p"}:
        theta = parse_complex_token(parts[2]).real
        apply_1q(state, n, int(parts[1]), np.array([[1, 0], [0, np.exp(1j * theta)]], dtype=np.complex128))
    elif gate == "rz":
        theta = parse_complex_token(parts[2]).real
        apply_1q(state, n, int(parts[1]), np.array([[np.exp(-0.5j * theta), 0], [0, np.exp(0.5j * theta)]], dtype=np.complex128))
    elif gate in {"cx", "CX"}:
        apply_cx(state, n, int(parts[1]), int(parts[2]))
    elif gate == "cz":
        apply_cz(state, n, int(parts[1]), int(parts[2]))
    elif gate == "swap":
        apply_swap(state, n, int(parts[1]), int(parts[2]))
    elif gate == "ccx":
        apply_ccx(state, n, int(parts[1]), int(parts[2]), int(parts[3]))
    elif gate == "ccz":
        apply_ccz(state, n, int(parts[1]), int(parts[2]), int(parts[3]))
    elif gate == "ecr":
        q0 = int(parts[1])
        q1 = int(parts[2])
        apply_1q(state, n, q0, S)
        apply_1q(state, n, q1, SDG)
        apply_1q(state, n, q1, H)
        apply_1q(state, n, q1, SDG)
        apply_cx(state, n, q0, q1)
        apply_1q(state, n, q0, X)
    elif gate == "pauli":
        paulis = parts[1]
        qsize = int(parts[2])
        qubits = [int(q) for q in parts[3:3 + qsize]]
        state[:] = apply_pauli_string(state, n, paulis, qubits)
    elif gate in {"id", "delay", "u0"}:
        return
    else:
        raise ValueError(f"unsupported gate {gate!r}")


def read_circuit(path: Path) -> tuple[int, list[list[str]]]:
    lines = [line.split() for line in path.read_text().splitlines() if line.strip()]
    n, gate_count = int(lines[0][0]), int(lines[0][1])
    ops = lines[1:]
    if len(ops) != gate_count:
        raise ValueError(f"{path}: header gate count {gate_count}, actual {len(ops)}")
    return n, ops


def run_circuit(path: Path, max_qubits: int) -> tuple[int, int, np.ndarray]:
    n, ops = read_circuit(path)
    if n > max_qubits:
        raise ValueError(f"{path}: n={n} exceeds --max-qubits={max_qubits}")
    state = np.zeros(1 << n, dtype=np.complex128)
    state[0] = 1.0
    for op in ops:
        apply_gate(state, n, op)
    return n, len(ops), state


def parse_observable(path: Path):
    tokens = [line.split() for line in path.read_text().splitlines() if line.strip()]
    first = tokens[0][0]
    if first in {"projector", "rank1_projector"}:
        op_count = int(tokens[0][1])
        return {"type": "projector", "coeff": 1.0, "ops": tokens[1:1 + op_count]}
    n_terms = int(first)
    terms = []
    cursor = 1
    for _ in range(n_terms):
        coeff = float(tokens[cursor][0])
        name = tokens[cursor][1]
        if name in {"projector", "rank1_projector"}:
            op_count = int(tokens[cursor][2])
            terms.append({"type": "projector", "coeff": coeff, "ops": tokens[cursor + 1:cursor + 1 + op_count]})
            cursor += 1 + op_count
        elif name == "pauli":
            paulis = tokens[cursor][2]
            qsize = int(tokens[cursor][3])
            qubits = [int(q) for q in tokens[cursor][4:4 + qsize]]
            terms.append({"type": "pauli", "coeff": coeff, "paulis": paulis, "qubits": qubits})
            cursor += 1
        else:
            raise ValueError(f"unsupported observable term {name!r}")
    return {"type": "sum", "terms": terms}


def projector_value(state: np.ndarray, n: int, coeff: float, prep_ops: Iterable[list[str]], max_qubits: int) -> complex:
    phi = np.zeros(1 << n, dtype=np.complex128)
    phi[0] = 1.0
    for op in prep_ops:
        apply_gate(phi, n, op)
    return coeff * abs(np.vdot(phi, state)) ** 2


def expectation(state: np.ndarray, n: int, observable, max_qubits: int) -> complex:
    if observable["type"] == "projector":
        return projector_value(state, n, observable["coeff"], observable["ops"], max_qubits)
    total = 0.0 + 0.0j
    for term in observable["terms"]:
        if term["type"] == "pauli":
            pstate = apply_pauli_string(state, n, term["paulis"], term["qubits"])
            total += term["coeff"] * np.vdot(state, pstate)
        elif term["type"] == "projector":
            total += projector_value(state, n, term["coeff"], term["ops"], max_qubits)
        else:
            raise ValueError(f"unsupported observable type {term['type']!r}")
    return total


def circuit_ids(testcase: Path) -> list[int]:
    info = testcase / "info.txt"
    first = info.read_text().split()[0]
    count = int(first)
    return list(range(count))


def main() -> None:
    parser = argparse.ArgumentParser(description="Exact dense statevector ground truth for small testcases")
    parser.add_argument("testcase", type=Path, help="testcase directory, e.g. testcases/N5D2")
    parser.add_argument(
        "--output",
        type=Path,
        help="CSV output path; defaults to <testcase>/ground_truth_statevector.csv",
    )
    parser.add_argument("--stdout", action="store_true", help="write CSV to stdout instead of a file")
    parser.add_argument("--max-qubits", type=int, default=20, help="refuse dense simulation above this size")
    parser.add_argument("--limit", type=int, help="only evaluate the first LIMIT circuits")
    args = parser.parse_args()

    observable = parse_observable(args.testcase / "observable.txt")
    ids = circuit_ids(args.testcase)
    if args.limit is not None:
        ids = ids[: args.limit]

    rows = []
    for cid in ids:
        circuit_path = args.testcase / "circuits" / f"{cid}.txt"
        n, gate_count, state = run_circuit(circuit_path, args.max_qubits)
        value = expectation(state, n, observable, args.max_qubits)
        rows.append({
            "circuit_id": cid,
            "n_qubits": n,
            "gate_count": gate_count,
            "expectation_real": f"{value.real:.17g}",
            "expectation_imag": f"{value.imag:.17g}",
        })

    fieldnames = ["circuit_id", "n_qubits", "gate_count", "expectation_real", "expectation_imag"]
    output_path = args.output or (args.testcase / "ground_truth_statevector.csv")
    if not args.stdout:
        with output_path.open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)
    else:
        import sys
        writer = csv.DictWriter(sys.stdout, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
