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
    alpha: E::G2Affine,

    /// Commitment to the trapdoor element β · G1
    beta: E::G2Affine,

    /// Commitment to the trapdoor element δ · G1
    delta: E::G2Affine,

    /// Encoded public input polynomials:
    /// (β·A_i(τ) + α·B_i(τ) + C_i(τ)) / γ · G1
    committed_statements: Vec<E::G1Affine>,
}
