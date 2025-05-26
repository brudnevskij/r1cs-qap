# 💡 R1CS-to-QAP to Pinocchio Prover — A zk-SNARK Study Project

This project implements the core proving system pipeline used in zk-SNARKs, progressing from **Rank-1 Constraint Systems (R1CS)** to **Quadratic Arithmetic Programs (QAP)**, and ultimately toward the **Pinocchio** SNARK protocol.

---

## 🚀 Features

- ✅ R1CS representation with variables and constraints
- ✅ Conversion from R1CS → QAP using Lagrange interpolation
- ✅ Vanishing polynomial generation
- ✅ Witness satisfaction checking in QAP form
- ✅ Rust code with `arkworks`
- 🔜 Trusted setup (CRS) for Pinocchio
- 🔜 Polynomial commitment evaluation
- 🔜 Proof generation and verification (Pinocchio → Groth16)

---

## 🔬 Background

This project is based on the structure described in [Vitalik Buterin’s blog post on QAPs](https://medium.com/@VitalikButerin/quadratic-arithmetic-programs-from-zero-to-hero-f6d558cea649) and extends it with practical implementation goals:

> _"The idea behind zk-SNARKs is to transform computation into a form where proving and verifying correctness can be done with polynomials, pairings, and commitments."_

---

## 📐 Architecture Overview

```
  Circuit (x³ + x + 5 = 35)
        ↓
      R1CS
        ↓
      QAP (A(x), B(x), C(x), Z(x))
        ↓
  Witness Evaluation: A(s), B(s), C(s)
        ↓
     Compute H(x) = (A·B - C)/Z
        ↓
    Commitments & Pairings → [WIP]
```

---

## 🧱 Planned Extensions

- [ ] Pinocchio-style trusted setup
- [ ] Commitment to QAP evaluations at `s`
- [ ] Proof structure with G1/G2 elements
- [ ] Groth16 optimization (compressed proofs)
- [ ] CLI + serialization format (e.g., JSON circuit input)
- [ ] ZK option toggle (homomorphic blinding)
