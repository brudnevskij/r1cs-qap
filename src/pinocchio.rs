use crate::qap::QAP;
use ark_ec::{CurveGroup, pairing::Pairing};
use ark_ff::Field;
use ark_poly::Polynomial;
use itertools::{Itertools, izip};

/// EvaluationKey contains all prover-side CRS elements.
/// These are used to construct a proof for a given QAP instance.
pub struct EvaluationKey<G: CurveGroup> {
    /// Commitments to v_k(s) for intermediate (non-public) wires, randomized with r_v:
    /// g^{r_v · v_k(s)}
    pub v_s: Vec<G>,

    /// Commitments to w_k(s): g^{r_w · w_k(s)}
    pub w_s: Vec<G>,

    /// Commitments to y_k(s): g^{r_y · y_k(s)} where r_y = r_v · r_w
    pub y_s: Vec<G>,

    /// Commitments to α_v · v_k(s): g^{r_v · α_v · v_k(s)}
    pub alpha_v_s: Vec<G>,

    /// Commitments to α_w · w_k(s): g^{r_w · α_w · w_k(s)}
    pub alpha_w_s: Vec<G>,

    /// Commitments to α_y · y_k(s): g^{r_y · α_y · y_k(s)} where r_y = r_v · r_w
    pub alpha_y_s: Vec<G>,

    /// Commitments to s^i: g^{s^i} for i in 0..degree(t(x))
    pub powers_of_s: Vec<G>,

    /// Combined QAP evaluation for β randomizer:
    /// g^{β(r_v · v_k(s) + r_w · w_k(s) + r_y · y_k(s))} for k ∈ intermediate wires
    pub beta_sum_g: Vec<G>,
}

impl<G: CurveGroup> EvaluationKey<G> {
    // TODO: consider using just Polynomial for QAP_polys
    fn new(g: G, toxic_waste: &Trapdoor<G::ScalarField>, qap: &QAP<G::ScalarField>) -> Self {
        let mut v_s = Vec::<G>::new();
        let mut w_s = Vec::<G>::new();
        let mut y_s = Vec::<G>::new();
        let mut alpha_v_s = Vec::<G>::new();
        let mut alpha_w_s = Vec::<G>::new();
        let mut alpha_y_s = Vec::<G>::new();
        let mut beta_sum_g = Vec::<G>::new();
        let mut powers_of_s = Vec::<G>::new();

        let ry = toxic_waste.rv * toxic_waste.rw;
        let g_v = g * toxic_waste.rv;
        let g_w = g * toxic_waste.rw;
        let g_y = g * ry;

        for (i, (v_poly, w_poly, y_poly)) in
            izip!(qap.a.iter(), qap.b.iter(), qap.c.iter()).enumerate()
        {
            // skipping public input variables
            if i <= qap.public_variables_count {
                continue;
            }

            let v_at_s_eval = v_poly.evaluate(&toxic_waste.s);
            v_s.push(g_v * v_at_s_eval);
            alpha_v_s.push(g_v * v_at_s_eval * toxic_waste.alpha_v);

            let w_at_s_eval = w_poly.evaluate(&toxic_waste.s);
            w_s.push(g_w * w_at_s_eval);
            alpha_w_s.push(g_w * w_at_s_eval * toxic_waste.alpha_w);

            let y_at_s_eval = y_poly.evaluate(&toxic_waste.s);
            y_s.push(g_y * y_at_s_eval);
            alpha_y_s.push(g_y * y_at_s_eval * toxic_waste.alpha_y);

            let beta = (g_v * &toxic_waste.beta * v_at_s_eval)
                + (g_w * &toxic_waste.beta * w_at_s_eval)
                + (g_y * &toxic_waste.beta * y_at_s_eval);
            beta_sum_g.push(beta);
        }

        for i in 0..=qap.target_poly.degree() {
            powers_of_s.push(g * &toxic_waste.s.pow([i as u64]));
        }

        Self {
            v_s,
            w_s,
            y_s,
            alpha_v_s,
            alpha_w_s,
            alpha_y_s,
            powers_of_s,
            beta_sum_g,
        }
    }
}

// TODO: consider scalar field or smth
pub struct Trapdoor<F> {
    pub rv: F,
    pub rw: F,
    pub s: F,
    pub alpha_v: F,
    pub alpha_w: F,
    pub alpha_y: F,
    pub beta: F,
    pub gamma: F,
}

/// VerificationKey contains all verifier-side CRS elements.
/// These are used to verify a proof generated from a QAP instance.
pub struct VerificationKey<G: CurveGroup> {
    /// Generator of the group: g ∈ G1
    pub g: G,

    /// Commitment to α_v: g^{α_v}
    pub g_alpha_v: G,

    /// Commitment to α_w: g^{α_w}
    pub g_alpha_w: G,

    /// Commitment to α_y: g^{α_y}
    pub g_alpha_y: G,

    /// Commitment to γ: g^{γ}
    pub g_gamma: G,

    /// Commitment to β * γ: g^{β * γ}
    pub g_beta_gamma: G,

    /// Commitment to the target polynomial t(s): g^{y(s)} for public inputs
    pub g_y_target_s: G,

    /// Commitments to the evaluated polynomials for each public input variable:
    /// Each tuple is:
    /// (
    ///   g^{r_v · v_i(s)},
    ///   g^{r_w · w_i(s)},
    ///   g^{r_y · y_i(s)}, where r_y = r_v · r_w
    /// )
    pub committed_input_polynomials: Vec<(G, G, G)>,
}

impl<G: CurveGroup> VerificationKey<G> {
    pub fn new(
        g: G,
        toxic_waste: &Trapdoor<G::ScalarField>,
        qap: &QAP<G::ScalarField>,
    ) -> VerificationKey<G> {
        let g_v = g * toxic_waste.rv;
        let g_w = g * toxic_waste.rw;
        let r_y = toxic_waste.rv * toxic_waste.rw;
        let g_y = g * r_y;

        let g_alpha_v = g * toxic_waste.alpha_v;
        let g_alpha_w = g * toxic_waste.alpha_w;
        let g_alpha_y = g * toxic_waste.alpha_y;
        let g_gamma = g * toxic_waste.gamma;
        let g_beta_gamma = g * toxic_waste.beta * toxic_waste.gamma;
        let g_y_target_s = g_y * qap.target_poly.evaluate(&toxic_waste.s);
        let committed_input_polynomials = izip!(
            qap.a.iter().take(qap.public_variables_count + 1),
            qap.b.iter().take(qap.public_variables_count + 1),
            qap.c.iter().take(qap.public_variables_count + 1)
        )
        .map(|(v_poly, w_poly, y_poly)| {
            let v = g_v * v_poly.evaluate(&toxic_waste.s);
            let w = g_w * w_poly.evaluate(&toxic_waste.s);
            let y = g_y * y_poly.evaluate(&toxic_waste.s);
            (v, w, y)
        })
        .collect();

        Self {
            g,
            g_alpha_v,
            g_alpha_w,
            g_alpha_y,
            g_gamma,
            g_beta_gamma,
            g_y_target_s,
            committed_input_polynomials,
        }
    }
}

struct CRS<G: CurveGroup> {
    evaluation_key: EvaluationKey<G>,
    verification_key: VerificationKey<G>,
}
impl<G: CurveGroup> CRS<G> {
    fn new(qap: &QAP<G::ScalarField>, toxic_waste: &Trapdoor<G::ScalarField>, g: G) -> Self {
        let evaluation_key = EvaluationKey::<G>::new(g, &toxic_waste, &qap);
        let verification_key = VerificationKey::<G>::new(g, &toxic_waste, &qap);
        Self {
            evaluation_key,
            verification_key,
        }
    }
}

/// `Proof` contains all prover-generated values to be sent to the verifier
/// in a Pinocchio-style zk-SNARK protocol.
pub struct Proof<G: CurveGroup> {
    /// Commitment to v(s): g^{v(s)}
    pub g_v_s: G,

    /// Commitment to w(s): g^{w(s)}
    pub g_w_s: G,

    /// Commitment to y(s): g^{y(s)} = v(s) * w(s)
    pub g_y_s: G,

    /// Commitment to the quotient polynomial h(s): g^{h(s)}
    /// where h(x) = (v(x)w(x) - y(x)) / t(x)
    pub g_h_s: G,

    /// Commitment to α_v · v(s): g^{α_v · v(s)}
    pub g_v_alpha_v_s: G,

    /// Commitment to α_w · w(s): g^{α_w · w(s)}
    pub g_w_alpha_w_s: G,

    /// Commitment to α_y · y(s): g^{α_y · y(s)}
    pub g_y_alpha_y_s: G,

    /// Commitment to β(v(s) + w(s) + y(s)) — combined using CRS values
    pub g_combined: G,
}

impl<G: CurveGroup> Proof<G> {
    pub fn new(
        evaluation_key: &EvaluationKey<G>,
        witness: &[G::ScalarField],
        qap: &QAP<G::ScalarField>,
    ) -> Self {
        let mut g_v_s = G::zero();
        let mut g_w_s = G::zero();
        let mut g_y_s = G::zero();
        let mut g_v_alpha_v_s = G::zero();
        let mut g_w_alpha_w_s = G::zero();
        let mut g_y_alpha_y_s = G::zero();
        let mut g_combined = G::zero();

        for (witness_i, &v, &w, &y, &alpha_v, &alpha_w, &alpha_y, &beta_sum) in izip!(
            witness.iter().skip(qap.public_variables_count + 1),
            evaluation_key.v_s.iter(),
            evaluation_key.w_s.iter(),
            evaluation_key.y_s.iter(),
            evaluation_key.alpha_v_s.iter(),
            evaluation_key.alpha_w_s.iter(),
            evaluation_key.alpha_y_s.iter(),
            evaluation_key.beta_sum_g.iter(),
        ) {
            g_v_s += v * witness_i;
            g_w_s += w * witness_i;
            g_y_s += y * witness_i;
            g_v_alpha_v_s += alpha_v * witness_i;
            g_w_alpha_w_s += alpha_w * witness_i;
            g_y_alpha_y_s += alpha_y * witness_i;
            g_combined += beta_sum * witness_i;
        }

        let h = qap.compute_h(witness);
        let mut g_h_s = G::zero();
        for (c, &g_s) in izip!(h.coeffs.iter(), evaluation_key.powers_of_s.iter()) {
            g_h_s += g_s * c;
        }

        Self {
            g_v_s,
            g_w_s,
            g_y_s,
            g_h_s,
            g_v_alpha_v_s,
            g_w_alpha_w_s,
            g_y_alpha_y_s,
            g_combined,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::r1cs::{
        R1CS,
        VariableType::{Intermediate, Private, Public},
    };
    use ark_bls12_381::{Fr, G1Projective};
    use ark_ec::PrimeGroup;
    use ark_ff::{One, UniformRand, Zero};

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

    fn dummy_trapdoor() -> Trapdoor<Fr> {
        Trapdoor {
            rv: Fr::from(2u64),
            rw: Fr::from(3u64),
            s: Fr::from(5u64),
            alpha_v: Fr::from(7u64),
            alpha_w: Fr::from(11u64),
            alpha_y: Fr::from(13u64),
            beta: Fr::from(17u64),
            gamma: Fr::from(19u64),
        }
    }

    #[test]
    fn test_eval_key_generation() {
        let qap = cubic_constraint_system();
        let trapdoor = dummy_trapdoor();
        let g = G1Projective::generator();
        let eval_key = EvaluationKey::<G1Projective>::new(g, &trapdoor, &qap);

        // Very basic sanity check — lengths must match
        assert!(!eval_key.v_s.is_empty());
        assert_eq!(eval_key.v_s.len(), eval_key.alpha_v_s.len());
        assert!(!eval_key.w_s.is_empty());
        assert_eq!(eval_key.w_s.len(), eval_key.alpha_w_s.len());
        assert!(!eval_key.y_s.is_empty());
        assert_eq!(eval_key.y_s.len(), eval_key.alpha_y_s.len());
        assert_eq!(eval_key.powers_of_s.len(), 5); // target poly is degree 4
    }

    #[test]
    fn test_verification_key_generation() {
        let qap = cubic_constraint_system();
        let trapdoor = dummy_trapdoor();
        let g = G1Projective::generator();

        let vk = VerificationKey::new(g, &trapdoor, &qap);

        assert_eq!(vk.committed_input_polynomials.len(), 2); // [~one, ~out]
        assert!(!vk.g.is_zero());
        assert!(!vk.g_y_target_s.is_zero());
    }

    #[test]
    fn test_proof_generation_outputs_nonzero_commitments() {
        let qap = cubic_constraint_system();
        let trapdoor = dummy_trapdoor();
        let g = G1Projective::generator();
        let crs = CRS::<G1Projective>::new(&qap, &trapdoor, g);

        // x = 2 → x^3 + x + 5 = 8 + 2 + 5 = 15
        // So ~out = 15
        let x = Fr::from(2u64);
        let x_sq = x * x;
        let x_cb = x_sq * x;
        let sym_1 = x_cb + x;
        let out = sym_1 + Fr::from(5u64);
        let witness = vec![Fr::one(), out, x, x_sq, x_cb, sym_1];
        assert!(qap.is_satisfied(&witness));

        let proof = Proof::<G1Projective>::new(&crs.evaluation_key, &witness, &qap);

        assert!(!proof.g_v_s.is_zero());
        assert!(!proof.g_w_s.is_zero());
        assert!(!proof.g_y_s.is_zero());
        assert!(!proof.g_h_s.is_zero());
        assert!(!proof.g_v_alpha_v_s.is_zero());
        assert!(!proof.g_w_alpha_w_s.is_zero());
        assert!(!proof.g_y_alpha_y_s.is_zero());
        assert!(!proof.g_combined.is_zero());
    }

    #[test]
    fn test_proof_generation_consistency_of_sizes() {
        let qap = cubic_constraint_system();
        let trapdoor = dummy_trapdoor();
        let g = G1Projective::generator();
        let ek = EvaluationKey::<G1Projective>::new(g, &trapdoor, &qap);

        // Valid witness for the R1CS: x = 2
        let x = Fr::from(2u64);
        let x_sq = x * x;
        let x_cb = x_sq * x;
        let sym_1 = x_cb + x;
        let out = sym_1 + Fr::from(5u64);
        let witness = vec![Fr::one(), out, x, x_sq, x_cb, sym_1, out];

        let proof = Proof::<G1Projective>::new(&ek, &witness, &qap);

        // This test isn't about values, but rather structural integrity
        let g1_identity = G1Projective::zero();
        assert_ne!(proof.g_v_s, g1_identity);
        assert_ne!(proof.g_combined, g1_identity);
    }
}
