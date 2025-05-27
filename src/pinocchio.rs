use crate::qap::QAP;
use ark_ec::{CurveGroup, pairing::Pairing};
use ark_ff::Field;
use ark_poly::Polynomial;
use ark_poly::univariate::DensePolynomial;
use itertools::izip;

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
    fn new<F: Field>(
        g: G,
        toxic_waste: Trapdoor<G::ScalarField>,
        qap: QAP<G::ScalarField>,
    ) -> Self {
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

struct CRS<G: CurveGroup> {
    evaluation_key: EvaluationKey<G>,
}
impl<G: CurveGroup> CRS<G> {
    fn new(qap: QAP<G::ScalarField>, toxic_waste: Trapdoor<G::ScalarField>, g: G) -> Self {
        let evaluation_key = EvaluationKey::<G>::new::<G::ScalarField>(g, toxic_waste, qap);
        Self{
            evaluation_key ,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_bls12_381::{Bls12_381, G1Projective, Fr};
    use ark_ec::PrimeGroup;
    use ark_ff::UniformRand;
    use ark_poly::DenseUVPolynomial;
    use ark_std::test_rng;

    #[test]
    fn test_eval_key_generation() {
        let mut rng = test_rng();

        // Fake 3-variable QAP for testing
        let a = vec![
            DensePolynomial::from_coefficients_vec(vec![Fr::from(1)]), // 1
            DensePolynomial::from_coefficients_vec(vec![Fr::from(2)]), // 2
            DensePolynomial::from_coefficients_vec(vec![Fr::from(3)]), // 3
        ];
        let b = a.clone();
        let c = a.clone();

        let target_poly = DensePolynomial::from_coefficients_vec(vec![Fr::from(1)]);

        let qap = QAP {
            a,
            b,
            c,
            target_poly,
            public_variables_count: 1, // skip 0 and 1
        };

        let toxic = Trapdoor {
            rv: Fr::rand(&mut rng),
            rw: Fr::rand(&mut rng),
            s: Fr::rand(&mut rng),
            alpha_v: Fr::rand(&mut rng),
            alpha_w: Fr::rand(&mut rng),
            alpha_y: Fr::rand(&mut rng),
            beta: Fr::rand(&mut rng),
            gamma: Fr::rand(&mut rng),
        };

        let g = G1Projective::generator();

        let eval_key = EvaluationKey::<G1Projective>::new::<Fr>(g, toxic, qap);

        // Very basic sanity check — lengths must match
        assert!(!eval_key.v_s.is_empty());
        assert_eq!(eval_key.v_s.len(), eval_key.alpha_v_s.len());
        assert_eq!(eval_key.powers_of_s.len(), 1); // target poly is degree 1
    }
}
