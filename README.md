# 💡 R1CS-to-QAP to Pinocchio Prover — A zk-SNARK Study Project

This project implements the core proving pipeline used in zk-SNARKs, progressing from **Rank-1 Constraint Systems (R1CS)** to **Quadratic Arithmetic Programs (QAP)** to **Pinocchio**, and ultimately toward the **Groth16** SNARK protocol..

---

## 🚀 Features

- ✅ R1CS representation with variable and constraint definitions
- ✅ Conversion from R1CS → QAP using Lagrange interpolation
- ✅ Vanishing polynomial generation
- ✅ QAP witness satisfaction checking
- ✅ Pinocchio-style CRS generation (Evaluation & Verification Keys)
- ✅ zk-SNARK proof generation and verification (Pinocchio)
- ✅ Fully tested with toy circuits (e.g., cubic)

---

## 🔬 Background

This project is inspired by [Vitalik Buterin’s blog post on QAPs](https://medium.com/@VitalikButerin/quadratic-arithmetic-programs-from-zero-to-hero-f6d558cea649) and aims to provide a clean, modular Rust implementation using `arkworks` libraries.

> _"The idea behind zk-SNARKs is to transform computation into a form where proving and verifying correctness can be done with polynomials, pairings, and commitments."_

---

## 📐 Architecture Overview

```
  Circuit (e.g., x³ + x + 5 = 35)
        ↓
      R1CS
        ↓
      QAP (A(x), B(x), C(x), Z(x))
        ↓
  Witness Evaluation: A(s), B(s), C(s)
        ↓
     Compute H(x) = (A·B - C)/Z
        ↓
 Commitments & Pairing-based Verifier (Pinocchio)
```

---

## 🧱 Planned Extensions
- [ ] Transition to Groth16 (compressed CRS, single pairing check)

---

## 📦 Dependencies

- [`ark-ec`](https://docs.rs/ark-ec)
- [`ark-poly`](https://docs.rs/ark-poly)
- [`ark-ff`](https://docs.rs/ark-ff)
- [`ark-bls12-381`](https://docs.rs/ark-bls12-381)
- [`itertools`](https://docs.rs/itertools)

---

## 🧪 Running Tests

```bash
cargo test
```

Tests include:
- Basic cubic constraint system
- Proof soundness check (invalid proofs fail)
- CRS key generation checks

---

## 📚 Learning Resources

- [Vitalik Buterin – QAPs From Zero to Hero](https://medium.com/@VitalikButerin/quadratic-arithmetic-programs-from-zero-to-hero-f6d558cea649)
- [The MoonMath Manual](https://github.com/LeastAuthority/moonmath-manual)
- [zk-SNARKs in Rust by arkworks](https://arkworks.rs)

---

## 🛠️ Author

Built as a study project by a cryptography student, aiming to fully understand zkSNARKs from the ground up.

---