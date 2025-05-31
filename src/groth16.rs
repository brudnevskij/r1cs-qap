use crate::qap::QAP;
use ark_ec::pairing::Pairing;
use ark_ff::Field;
use ark_poly::Polynomial;
use itertools::izip;

struct TrapDoor<F> {
    alpha: F,
    beta: F,
    gamma: F,
    delta: F,
    tau: F,
}

struct CRS<E: Pairing> {
    proving_key: ProvingKey<E>,
    verification_key: VerificationKey<E>,
}

/// Proving key for the Groth16 zk-SNARK protocol.
struct ProvingKey<E: Pairing> {
    /// Commitment to the trapdoor element α · G1
    alpha: E::G1Affine,

    /// Commitment to the trapdoor element β · G1
    beta: E::G1Affine,

    /// Commitment to the trapdoor element δ · G1
    delta: E::G1Affine,

    /// Powers of τ · G1, i.e., [τ⁰·G1, τ¹·G1, ..., τⁿ·G1]
    /// Used to evaluate the polynomials A(x), B(x), and C(x) in the QAP
    tau_powers: Vec<E::G1Affine>,

    /// Commitments to witness polynomials evaluated at τ, divided by δ
    /// (β·A_i(τ) + α·B_i(τ) + C_i(τ)) / δ · G1
    committed_witnesses: Vec<E::G1Affine>,

    /// Commitments to the quotient polynomial powers: τ^i · t(τ) / δ · G1
    h_query: Vec<E::G1Affine>,
}

impl<E: Pairing> ProvingKey<E> {
    pub fn new(
        g: E::G1,
        trap_door: TrapDoor<E::ScalarField>,
        qap: QAP<E::ScalarField>,
    ) -> ProvingKey<E> {
        let TrapDoor {
            alpha,
            beta,
            delta,
            tau,
            ..
        } = trap_door;

        let t_at_tau = qap.target_poly.evaluate(&tau);

        let tau_powers: Vec<_> = (0..=qap.target_poly.degree() - 1)
            .map(|i| {
                let tau_i = tau.pow([i as u64]);
                (g * tau_i).into()
            })
            .collect();

        let committed_witnesses: Vec<_> = izip!(qap.a.iter(), qap.b.iter(), qap.c.iter())
            .skip(qap.public_variables_count + 1)
            .map(|(a, b, c)| {
                let a_tau = beta * a.evaluate(&tau);
                let b_tau = alpha * b.evaluate(&tau);
                let c_tau = c.evaluate(&tau);
                let coeff = (a_tau + b_tau + c_tau) / delta;
                (g * coeff).into()
            })
            .collect();

        let h_query: Vec<_> = (0..=qap.target_poly.degree() - 2)
            .map(|i| {
                let scaled = tau.pow([i as u64]) * t_at_tau / delta;
                (g * scaled).into()
            })
            .collect();

        Self {
            alpha: (g * trap_door.alpha).into(),
            beta: (g * trap_door.beta).into(),
            delta: (g * trap_door.delta).into(),
            tau_powers,
            committed_witnesses,
            h_query,
        }
    }
}

/// Verification key for the Groth16 zk-SNARK protocol.
struct VerificationKey<E: Pairing> {
    /// Commitment to the trapdoor element α · G1
    alpha_g1: E::G1Affine,

    /// Commitment to the trapdoor element β · G1
    beta_g2: E::G2Affine,

    /// Commitment to the trapdoor element δ · G1
    delta_g2: E::G2Affine,

    /// Commitment to the trapdoor element γ · G2
    gamma_g2: E::G2Affine,

    /// Encoded public input polynomials:
    /// (β·A_i(τ) + α·B_i(τ) + C_i(τ)) / γ · G1
    gamma_abc_g1: Vec<E::G1Affine>,
}

impl<E: Pairing> VerificationKey<E> {
    pub fn new(
        g1: E::G1,
        g2: E::G2,
        trap_door: &TrapDoor<E::ScalarField>,
        qap: &QAP<E::ScalarField>,
    ) -> VerificationKey<E> {
        let TrapDoor {
            alpha,
            beta,
            gamma,
            delta,
            tau,
        } = trap_door;

        let gamma_abc_g1 = izip!(qap.a.iter(), qap.b.iter(), qap.c.iter())
            .take(qap.public_variables_count + 1)
            .map(|(a, b, c)| {
                let a_tau = *beta * a.evaluate(tau);
                let b_tau = *alpha * b.evaluate(tau);
                let c_tau = c.evaluate(tau);
                let coeff = (a_tau + b_tau + c_tau) / *gamma;
                (g1 * coeff).into()
            })
            .collect();

        Self {
            alpha_g1: (g1 * alpha).into(),
            beta_g2: (g2 * beta).into(),
            delta_g2: (g2 * delta).into(),
            gamma_g2: (g2 * gamma).into(),
            gamma_abc_g1,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::r1cs::R1CS;
    use crate::r1cs::VariableType::{Intermediate, Private, Public};
    use ark_bls12_381::{Bls12_381, Fr as F, Fr};
    use ark_ec::PrimeGroup;
    use ark_ff::{One, Zero};
    use ark_std::UniformRand;
    use rand::thread_rng;

    fn cubic_constraint_system() -> QAP<Fr> {
        // Create constraints for x**3 + x + 5 = 35
        let mut r1cs = R1CS::<Fr>::new();
        r1cs.add_variable("x".to_string(), Private);
        r1cs.add_variable("x_sq".to_string(), Intermediate);
        r1cs.add_variable("x_cb".to_string(), Intermediate);
        r1cs.add_variable("sym_1".to_string(), Intermediate);
        r1cs.add_variable("~out".to_string(), Public);

        // vars = [~one, ~out, x, x_sq, x_cb, sym_1]
        // x*x = x_sq
        let a = vec![
            Fr::zero(),
            Fr::zero(),
            Fr::one(),
            Fr::zero(),
            Fr::zero(),
            Fr::zero(),
        ];
        let b = vec![
            Fr::zero(),
            Fr::zero(),
            Fr::one(),
            Fr::zero(),
            Fr::zero(),
            Fr::zero(),
        ];
        let c = vec![
            Fr::zero(),
            Fr::zero(),
            Fr::zero(),
            Fr::one(),
            Fr::zero(),
            Fr::zero(),
        ];
        r1cs.add_constraint(a, b, c);

        // x_sq * x = x_cb
        let a = vec![
            Fr::zero(),
            Fr::zero(),
            Fr::zero(),
            Fr::one(),
            Fr::zero(),
            Fr::zero(),
        ];
        let b = vec![
            Fr::zero(),
            Fr::zero(),
            Fr::one(),
            Fr::zero(),
            Fr::zero(),
            Fr::zero(),
        ];
        let c = vec![
            Fr::zero(),
            Fr::zero(),
            Fr::zero(),
            Fr::zero(),
            Fr::one(),
            Fr::zero(),
        ];
        r1cs.add_constraint(a, b, c);

        // x_cb + x = sym_1
        let a = vec![
            Fr::zero(),
            Fr::zero(),
            Fr::one(),
            Fr::zero(),
            Fr::one(),
            Fr::zero(),
        ];
        let b = vec![
            Fr::one(),
            Fr::zero(),
            Fr::zero(),
            Fr::zero(),
            Fr::zero(),
            Fr::zero(),
        ];
        let c = vec![
            Fr::zero(),
            Fr::zero(),
            Fr::zero(),
            Fr::zero(),
            Fr::zero(),
            Fr::one(),
        ];
        r1cs.add_constraint(a, b, c);

        // sym_1 + 5 = out
        let a = vec![
            Fr::from(5),
            Fr::zero(),
            Fr::zero(),
            Fr::zero(),
            Fr::zero(),
            Fr::one(),
        ];
        let b = vec![
            Fr::one(),
            Fr::zero(),
            Fr::zero(),
            Fr::zero(),
            Fr::zero(),
            Fr::zero(),
        ];
        let c = vec![
            Fr::zero(),
            Fr::one(),
            Fr::zero(),
            Fr::zero(),
            Fr::zero(),
            Fr::zero(),
        ];
        r1cs.add_constraint(a, b, c);
        QAP::from(&r1cs)
    }

    #[test]
    fn test_proving_key_generation() {
        let mut rng = thread_rng();

        let g = <Bls12_381 as Pairing>::G1::generator();

        // Random trapdoor values
        let trapdoor = TrapDoor {
            alpha: F::rand(&mut rng),
            beta: F::rand(&mut rng),
            gamma: F::rand(&mut rng),
            delta: F::rand(&mut rng),
            tau: F::rand(&mut rng),
        };
        let qap = cubic_constraint_system();
        let pk = ProvingKey::<Bls12_381>::new(g, trapdoor, qap);

        // Expect tau_powers = degree(target_poly) = 4
        assert_eq!(pk.tau_powers.len(), 4, "Incorrect number of tau powers");

        // Expect h_query = degree(target_poly) - 1 = 3
        assert_eq!(pk.h_query.len(), 3, "Incorrect number of h_query elements");

        // Expect committed_witnesses = #vars - public - 1 (for 1-based indexing)
        assert_eq!(
            pk.committed_witnesses.len(),
            4,
            "Incorrect witness commitments"
        );
    }

    #[test]
    fn test_verification_key_generation() {
        let mut rng = thread_rng();

        // Generator elements
        let g1 = <Bls12_381 as Pairing>::G1::generator();
        let g2 = <Bls12_381 as Pairing>::G2::generator();

        // Random trapdoor for trusted setup
        let trapdoor = TrapDoor {
            alpha: F::rand(&mut rng),
            beta: F::rand(&mut rng),
            gamma: F::rand(&mut rng),
            delta: F::rand(&mut rng),
            tau: F::rand(&mut rng),
        };

        // Build cubic constraint system → QAP
        let qap = cubic_constraint_system();

        // Construct verification key
        let vk = VerificationKey::<Bls12_381>::new(g1, g2, &trapdoor, &qap);

        // Check that the number of gamma_abc elements matches number of public vars + 1 (for ~1)
        let expected_gamma_abc_len = qap.public_variables_count + 1;
        assert_eq!(
            vk.gamma_abc_g1.len(),
            expected_gamma_abc_len,
            "gamma_abc_g1 should have one entry per public variable + 1"
        );
    }

}
