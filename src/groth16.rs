use crate::qap::QAP;
use ark_ec::CurveGroup;
use ark_ec::pairing::Pairing;
use ark_ff::{Field, Zero};
use ark_poly::Polynomial;
use ark_poly::univariate::DensePolynomial;
use itertools::izip;

struct TrapDoor<F> {
    alpha: F,
    beta: F,
    gamma: F,
    delta: F,
    tau: F,
}

struct CRS<E: Pairing> {
    sigma1: Sigma1<E>,
    sigma2: Sigma2<E>,
}

impl<E: Pairing> CRS<E> {
    pub fn new(
        g1: E::G1,
        g2: E::G2,
        trap_door: &TrapDoor<E::ScalarField>,
        qap: &QAP<E::ScalarField>,
    ) -> Self {
        Self {
            sigma1: Sigma1::new(g1, trap_door, qap),
            sigma2: Sigma2::new(g2, trap_door, qap),
        }
    }
}

/// Proving key for the Groth16 zk-SNARK protocol.
struct Sigma1<E: Pairing> {
    /// Commitment to the trapdoor element α · G1
    alpha: E::G1Affine,

    /// Commitment to the trapdoor element β · G1
    beta: E::G1Affine,

    /// Commitment to the trapdoor element δ · G1
    delta: E::G1Affine,

    /// Powers of τ · G1, i.e., [τ⁰·G1, τ¹·G1, ..., τⁿ·G1]
    /// Used to evaluate the polynomials A(x), B(x), and C(x) in the QAP
    tau_powers: Vec<E::G1Affine>,

    /// Commitments to statements polynomials evaluated at τ, divided by γ
    committed_statements: Vec<E::G1Affine>,

    /// Commitments to witness polynomials evaluated at τ, divided by δ
    /// (β·A_i(τ) + α·B_i(τ) + C_i(τ)) / δ · G1
    committed_witnesses: Vec<E::G1Affine>,

    /// Commitments to the quotient polynomial powers: τ^i · t(τ) / δ · G1
    h_query: Vec<E::G1Affine>,
}

impl<E: Pairing> Sigma1<E> {
    pub fn new(
        g: E::G1,
        trap_door: &TrapDoor<E::ScalarField>,
        qap: &QAP<E::ScalarField>,
    ) -> Sigma1<E> {
        let TrapDoor {
            alpha,
            beta,
            gamma,
            delta,
            tau,
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
                let a_tau = *beta * a.evaluate(&tau);
                let b_tau = *alpha * b.evaluate(&tau);
                let c_tau = c.evaluate(&tau);
                let coeff = (a_tau + b_tau + c_tau) / *delta;
                (g * coeff).into()
            })
            .collect();

        let committed_statements = izip!(qap.a.iter(), qap.b.iter(), qap.c.iter())
            .take(qap.public_variables_count + 1)
            .map(|(a, b, c)| {
                let a_tau = *beta * a.evaluate(tau);
                let b_tau = *alpha * b.evaluate(tau);
                let c_tau = c.evaluate(tau);
                let coeff = (a_tau + b_tau + c_tau) / *gamma;
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
            committed_statements,
            committed_witnesses,
            h_query,
        }
    }
}

/// Verification key for the Groth16 zk-SNARK protocol.
struct Sigma2<E: Pairing> {
    /// Commitment to the trapdoor element α · G1
    alpha: E::G2Affine,

    /// Commitment to the trapdoor element β · G1
    beta: E::G2Affine,

    /// Commitment to the trapdoor element δ · G1
    delta: E::G2Affine,

    /// Commitment to the trapdoor element γ · G2
    gamma: E::G2Affine,

    tau_powers: Vec<E::G2Affine>,
}

impl<E: Pairing> Sigma2<E> {
    pub fn new(
        g: E::G2,
        trap_door: &TrapDoor<E::ScalarField>,
        qap: &QAP<E::ScalarField>,
    ) -> Sigma2<E> {
        let TrapDoor {
            alpha,
            beta,
            gamma,
            delta,
            tau,
        } = trap_door;

        let tau_powers = (0..=qap.target_poly.degree() - 1)
            .map(|i| (g * tau.pow([i as u64])).into())
            .collect();

        Self {
            alpha: (g * alpha).into(),
            beta: (g * beta).into(),
            delta: (g * delta).into(),
            gamma: (g * gamma).into(),
            tau_powers,
        }
    }
}

pub fn naive_msm<G: CurveGroup>(bases: &[G::Affine], scalars: &[G::ScalarField]) -> G {
    izip!(scalars, bases).map(|(s, b)| *b * s).sum()
}

struct Proof<E: Pairing> {
    a: E::G1Affine,
    b: E::G2Affine,
    c: E::G1Affine,
}

impl<E: Pairing> Proof<E> {
    pub fn new(
        r: E::ScalarField,
        s: E::ScalarField,
        crs: &CRS<E>,
        qap: &QAP<E::ScalarField>,
        assignment: &[E::ScalarField],
    ) -> Proof<E> {
        let CRS { sigma1, sigma2 } = crs;

        // 1. Compute A = α + Σ a_i * u_i(τ) + rδ
        let a_query = qap
            .a
            .iter()
            .zip(assignment)
            .map(|(poly, coeff)| {
                let eval = Self::eval_poly_at_tau::<E::G1>(poly.clone(), &sigma1.tau_powers);
                eval * *coeff
            })
            .sum::<E::G1>();

        let a = sigma1.alpha + a_query + sigma1.delta * r;

        // 2. Compute B in G2: β + Σ b_i * v_i(τ) + sδ
        let b_query_g2 = qap
            .b
            .iter()
            .zip(assignment)
            .map(|(poly, coeff)| {
                let eval = Self::eval_poly_at_tau::<E::G2>(poly.clone(), &sigma2.tau_powers);
                eval * *coeff
            })
            .sum::<E::G2>();

        let b = sigma2.beta + b_query_g2 + sigma2.delta * s;

        // 3. Compute B in G1 (needed for C): β + Σ b_i * v_i(τ) + sδ
        let b_query_g1 = qap
            .b
            .iter()
            .zip(assignment)
            .map(|(poly, coeff)| {
                let eval = Self::eval_poly_at_tau::<E::G1>(poly.clone(), &sigma1.tau_powers);
                eval * *coeff
            })
            .sum::<E::G1>();

        let b_c = sigma1.beta + b_query_g1 + sigma1.delta * s;

        // 4. Compute H(τ)
        let h_poly = qap.compute_h(assignment);
        let h_term = Self::eval_poly_at_tau::<E::G1>(h_poly, &sigma1.h_query);

        // 5. Compute witness commitments
        let witness_scalars: Vec<_> = assignment
            .iter()
            .skip(qap.public_variables_count + 1)
            .cloned()
            .collect();

        let witness_term: E::G1 = naive_msm(&sigma1.committed_witnesses, &witness_scalars);

        // 6. Compute C = witness_term + h_term + A·s + B_c·s - (r·s·δ)
        let c: E::G1 = witness_term + h_term + (a * s) + (b_c * s) - (sigma1.delta * r * s);

        Proof {
            a: a.into(),
            b: b.into(),
            c: c.into(),
        }
    }

    fn eval_poly_at_tau<G: CurveGroup>(
        poly: DensePolynomial<G::ScalarField>,
        tau_powers: &[G::Affine],
    ) -> G {
        izip!(poly.coeffs.iter(), tau_powers)
            .map(|(&x, &tau)| tau * x)
            .sum()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::r1cs::R1CS;
    use crate::r1cs::VariableType::{Intermediate, Private, Public};
    use ark_bls12_381::{Bls12_381, Fr as F, Fr};
    use ark_ec::{AffineRepr, PrimeGroup};
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
    fn test_sigma1_generation() {
        let mut rng = thread_rng();
        let g1 = <Bls12_381 as Pairing>::G1::generator();

        let trapdoor = TrapDoor {
            alpha: Fr::rand(&mut rng),
            beta: Fr::rand(&mut rng),
            gamma: Fr::rand(&mut rng),
            delta: Fr::rand(&mut rng),
            tau: Fr::rand(&mut rng),
        };

        let qap = cubic_constraint_system();
        let sigma1 = Sigma1::<Bls12_381>::new(g1, &trapdoor, &qap);

        assert_eq!(
            sigma1.tau_powers.len(),
            qap.target_poly.degree(),
            "tau powers length should equal target poly degree"
        );

        assert_eq!(
            sigma1.h_query.len(),
            qap.target_poly.degree() - 1,
            "h_query length should equal target poly degree - 1"
        );

        assert_eq!(
            sigma1.committed_statements.len(),
            qap.public_variables_count + 1,
            "committed_statements should include ~1 and all public vars"
        );

        let expected_witnesses = qap.a.len() - (qap.public_variables_count + 1);
        assert_eq!(
            sigma1.committed_witnesses.len(),
            expected_witnesses,
            "committed_witnesses length incorrect"
        );
    }

    #[test]
    fn test_sigma2_generation() {
        let mut rng = thread_rng();
        let g2 = <Bls12_381 as Pairing>::G2::generator();

        let trapdoor = TrapDoor {
            alpha: Fr::rand(&mut rng),
            beta: Fr::rand(&mut rng),
            gamma: Fr::rand(&mut rng),
            delta: Fr::rand(&mut rng),
            tau: Fr::rand(&mut rng),
        };

        let qap = cubic_constraint_system();
        let sigma2 = Sigma2::<Bls12_381>::new(g2, &trapdoor, &qap);

        assert_eq!(
            sigma2.tau_powers.len(),
            qap.target_poly.degree(),
            "sigma2 tau powers should match target poly degree"
        );
    }

    #[test]
    fn test_groth16_proof_generation() {
        use ark_bls12_381::{Bls12_381, Fr as F};
        use ark_ec::{pairing::Pairing, CurveGroup};
        use ark_std::UniformRand;
        use rand::thread_rng;

        // 1. Generate a QAP from your example constraint system
        let qap = cubic_constraint_system(); // x^3 + x + 5 = 35

        // 2. Build the assignment vector
        // assignment = [1, 35, x, x^2, x^3, x^3 + x]
        let x = F::from(3u32);
        let x_sq = x * x;
        let x_cb = x_sq * x;
        let sym_1 = x_cb + x;
        let out = sym_1 + F::from(5);
        let assignment = vec![F::one(), out, x, x_sq, x_cb, sym_1];

        // 3. Create trapdoor (toxic waste)
        let mut rng = thread_rng();
        let trapdoor = TrapDoor {
            alpha: F::rand(&mut rng),
            beta: F::rand(&mut rng),
            gamma: F::rand(&mut rng),
            delta: F::rand(&mut rng),
            tau: F::rand(&mut rng),
        };

        // 4. Generate CRS
        let g1 = <Bls12_381 as Pairing>::G1::generator();
        let g2 = <Bls12_381 as Pairing>::G2::generator();
        let crs = CRS::<Bls12_381>::new(g1, g2, &trapdoor, &qap);

        // 5. Generate random r, s for proof
        let r = F::rand(&mut rng);
        let s = F::rand(&mut rng);

        // 6. Generate the proof
        let proof = Proof::<Bls12_381>::new(r, s, &crs, &qap, &assignment);

        // 7. Assertions
        assert!(!proof.a.is_zero(), "A is zero");
        assert!(!proof.b.is_zero(), "B is zero");
        assert!(!proof.c.is_zero(), "C is zero");

        // NOTE: You could add a pairing check here later once a verifier is implemented.
        println!("Proof successfully generated.");
    }
}
