use crate::qap::QAP;
use crate::r1cs::R1CS;
use crate::r1cs::VariableType::{Intermediate, Private, Public};
use ark_bls12_381::{Bls12_381, Fr};
use ark_ec::PrimeGroup;
use ark_ff::{One, UniformRand, Zero};

use crate::groth16::{
    CRS as Groth16CRS, Proof as Groth16Proof, TrapDoor as GrothTrapdoor, TrapDoor,
    verify as verify_groth16,
};
use crate::pinocchio::{
    CRS as PinocchioCRS, Proof as PinocchioProof, Trapdoor as PinocchioTrapdoor,
    verify as verify_pinocchio,
};

mod groth16;
mod pinocchio;
mod qap;
mod r1cs;

/// Returns circuit for x² - x + 132 = out
fn get_square_circuit() -> R1CS<Fr> {
    let mut r1cs = R1CS::new();
    r1cs.add_variable("x".to_string(), Private);
    r1cs.add_variable("x_sq".to_string(), Public);
    r1cs.add_variable("sym_1".to_string(), Intermediate);
    r1cs.add_variable("~out".to_string(), Public);

    // Constraint 1: x * x = x_sq
    r1cs.add_constraint(
        vec![Fr::zero(), Fr::zero(), Fr::zero(), Fr::one(), Fr::zero()],
        vec![Fr::zero(), Fr::zero(), Fr::zero(), Fr::one(), Fr::zero()],
        vec![Fr::zero(), Fr::one(), Fr::zero(), Fr::zero(), Fr::zero()],
    );

    // Constraint 2: x_sq - x = sym_1
    r1cs.add_constraint(
        vec![Fr::zero(), Fr::one(), Fr::zero(), Fr::from(-1), Fr::zero()],
        vec![Fr::one(), Fr::zero(), Fr::zero(), Fr::zero(), Fr::zero()],
        vec![Fr::zero(), Fr::zero(), Fr::zero(), Fr::zero(), Fr::one()],
    );

    // Constraint 3: sym_1 + 132 = ~out
    r1cs.add_constraint(
        vec![Fr::from(132), Fr::zero(), Fr::zero(), Fr::zero(), Fr::one()],
        vec![Fr::one(), Fr::zero(), Fr::zero(), Fr::zero(), Fr::zero()],
        vec![Fr::zero(), Fr::zero(), Fr::one(), Fr::zero(), Fr::zero()],
    );

    r1cs
}

fn main() {
    println!("📐 Starting zk-SNARK pipeline demo...");

    // 1. Build R1CS
    let r1cs = get_square_circuit();
    println!("Proving:\n{}", r1cs);

    // 2. Convert to QAP
    let qap = QAP::from(&r1cs);

    // 3. Build witness
    let x = Fr::from(147);
    let x_sq = x * x;
    let sym_1 = x_sq - x;
    let out = sym_1 + Fr::from(132);
    let witness = vec![Fr::one(), x_sq, out, x, sym_1];
    let public_input = vec![Fr::one(), x_sq, out];

    assert!(qap.is_satisfied(&witness), "❌ QAP not satisfied!");
    println!("✅ QAP is satisfied by witness");

    let rng = &mut ark_std::test_rng();

    // ==========================
    // 4. Pinocchio proving
    // ==========================
    println!("📦 Pinocchio setup...");
    let toxic_pino = PinocchioTrapdoor {
        rv: Fr::rand(rng),
        rw: Fr::rand(rng),
        s: Fr::rand(rng),
        alpha_v: Fr::rand(rng),
        alpha_w: Fr::rand(rng),
        alpha_y: Fr::rand(rng),
        beta: Fr::rand(rng),
        gamma: Fr::rand(rng),
    };
    let crs_pino = PinocchioCRS::<Bls12_381>::new(
        &qap,
        &toxic_pino,
        <Bls12_381 as ark_ec::pairing::Pairing>::G1::generator(),
        <Bls12_381 as ark_ec::pairing::Pairing>::G2::generator(),
    );

    println!("🔏 Generating Pinocchio proof...");
    let proof_pino = PinocchioProof::<Bls12_381>::new(&crs_pino.evaluation_key, &witness, &qap);

    println!("🔍 Verifying Pinocchio proof...");
    let ok_pino = verify_pinocchio(proof_pino, crs_pino.verification_key, &public_input);
    println!(
        "📄 Pinocchio proof: {}",
        if ok_pino {
            "ACCEPTED ✅"
        } else {
            "REJECTED ❌"
        }
    );

    // ==========================
    // 5. Groth16 proving
    // ==========================
    println!("\n📦 Groth16 setup...");
    let toxic_groth = GrothTrapdoor {
        alpha: Fr::rand(rng),
        beta: Fr::rand(rng),
        gamma: Fr::rand(rng),
        delta: Fr::rand(rng),
        tau: Fr::rand(rng),
    };
    let g1 = <Bls12_381 as ark_ec::pairing::Pairing>::G1::generator();
    let g2 = <Bls12_381 as ark_ec::pairing::Pairing>::G2::generator();

    let crs_groth = Groth16CRS::<Bls12_381>::new(g1, g2, &toxic_groth, &qap);

    println!("🔏 Generating Groth16 proof...");
    let proof_groth = Groth16Proof::new(Fr::rand(rng), Fr::rand(rng), &crs_groth, &qap, &witness);

    println!("🔍 Verifying Groth16 proof...");
    let ok_groth = verify_groth16(&crs_groth, &proof_groth, &public_input);
    println!(
        "📄 Groth16 proof: {}",
        if ok_groth {
            "ACCEPTED ✅"
        } else {
            "REJECTED ❌"
        }
    );
}
