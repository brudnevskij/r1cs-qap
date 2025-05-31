# 💡 R1CS → QAP → Pinocchio → Groth16, A zk-SNARK Study Project

This project implements the core proving pipeline used in zk-SNARKs, progressing step-by-step from **Rank-1 Constraint Systems (R1CS)** to **Quadratic Arithmetic Programs (QAP)**, to **Pinocchio**, and finally toward the **Groth16** proving system.

---
## 🚀 Features

- ✅ R1CS representation with variable and constraint definitions
- ✅ Conversion from R1CS → QAP using Lagrange interpolation
- ✅ Pinocchio implementation
- ✅ Groth16 implementation
- ✅ Tested with toy circuits (e.g., cubic polynomial)

---

## 🔬 Background

This project was created as a deep learning exercise to understand and implement the Groth16 zk-SNARK proving system from the ground up. Starting with the formulation of arithmetic circuits and their transformation into Rank-1 Constraint Systems (R1CS), I then implemented the conversion to Quadratic Arithmetic Programs (QAP) using Lagrange interpolation.
From there, I built the Pinocchio protocol to understand how pairing-based zk-SNARKs work, and finally extended it to support the Groth16 protocol,
optimizing proof generation and verification with a compressed Common Reference String (CRS). The goal was to reconstruct the full proving pipeline almost from scratch using the Rust arkworks ecosystem,
with a focus on clarity, correctness, and hands-on understanding of each transformation layer.

> _"The idea behind zk-SNARKs is to transform computation into a form where proving and verifying correctness can be done with polynomials, pairings, and commitments."_

---

## 📐 Architecture Overview

```
  Arithmetic Circuit (e.g., x³ + x + 5 = 35)
        ↓
      R1CS
        ↓
      QAP (A(x), B(x), C(x), Z(x))
        ↓
  Witness Evaluation: A(s), B(s), C(s)
        ↓
     Compute H(x) = (A·B - C)/Z
        ↓
  → Pinocchio: Pairing-based proof with 3 pairings
  → Groth16: Optimized zk-SNARK with 1 pairing + compressed CRS
```

## 🧪 Running the Demo

This project includes a complete end-to-end example in `main.rs` for both the **Pinocchio** and **Groth16** zk-SNARK protocols. It builds a toy arithmetic circuit:

> x² - x + 132 = out

Then it:

1. Translates the circuit to R1CS
2. Converts R1CS → QAP
3. Generates a trusted setup (Pinocchio & Groth16)
4. Constructs a witness
5. Produces a proof
6. Verifies it

---

### 🛠 Running

```bash
cargo run
```

You should see output like:
```bash
✅ QAP is satisfied by witness
📦 Pinocchio setup...
🔏 Generating Pinocchio proof...
📄 Pinocchio proof: ACCEPTED ✅
📦 Groth16 setup...
🔏 Generating Groth16 proof...
📄 Groth16 proof: ACCEPTED ✅
```

---

## 📦 Dependencies

- [`ark-ec`](https://docs.rs/ark-ec)
- [`ark-poly`](https://docs.rs/ark-poly)
- [`ark-ff`](https://docs.rs/ark-ff)
- [`ark-bls12-381`](https://docs.rs/ark-bls12-381)
- [`itertools`](https://docs.rs/itertools)

---


---

## 📚 Learning Resources

- [Vitalik Buterin – QAPs From Zero to Hero](https://medium.com/@VitalikButerin/quadratic-arithmetic-programs-from-zero-to-hero-f6d558cea649)
- [The MoonMath Manual](https://github.com/LeastAuthority/moonmath-manual)
- [arkworks zk-SNARKs in Rust](https://arkworks.rs)

---