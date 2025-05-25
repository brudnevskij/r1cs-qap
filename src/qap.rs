use crate::r1cs::R1CS;
use ark_ff::Field;
use ark_poly::DenseUVPolynomial;
use ark_poly::Polynomial;
use ark_poly::univariate::DensePolynomial;
use std::ops::Add;

struct QAP<F: Field> {
    target_poly: DensePolynomial<F>,
    a: Vec<DensePolynomial<F>>,
    b: Vec<DensePolynomial<F>>,
    c: Vec<DensePolynomial<F>>,
}
impl<F: Field> QAP<F> {
    pub fn new(target_poly: DensePolynomial<F>) -> QAP<F> {
        Self {
            target_poly,
            a: vec![],
            b: vec![],
            c: vec![],
        }
    }
}

pub fn get_vanishing_polynomial<F: Field>(roots: Vec<F>) -> DensePolynomial<F> {
    let mut polynomial = DensePolynomial::from_coefficients_vec(vec![F::one()]);

    for root in roots {
        let factor_polynomial = DensePolynomial::from_coefficients_vec(vec![root.neg(), F::one()]);
        polynomial = polynomial.naive_mul(&factor_polynomial);
    }
    polynomial
}

impl<F: Field> From<&R1CS<F>> for QAP<F> {
    fn from(value: &R1CS<F>) -> Self {
        let mut a = vec![];
        let mut b = vec![];
        let mut c = vec![];

        for var_index in 0..value.variables.len() {
            let mut a_domain_points = vec![];
            let mut a_values = vec![];

            let mut b_domain_points = vec![];
            let mut b_values = vec![];

            let mut c_domain_points = vec![];
            let mut c_values = vec![];

            for (i, constraint) in value.constraints.iter().enumerate() {
                a_domain_points.push(F::from(i as i32));
                b_domain_points.push(F::from(i as i32));
                c_domain_points.push(F::from(i as i32));

                a_values.push(constraint.a[var_index].clone());
                b_values.push(constraint.b[var_index].clone());
                c_values.push(constraint.c[var_index].clone());
            }
            a.push(interpolate_lagrange(&a_domain_points, &a_values));
            b.push(interpolate_lagrange(&b_domain_points, &b_values));
            c.push(interpolate_lagrange(&c_domain_points, &c_values));
        }

        let target_poly = get_vanishing_polynomial(
            (1..=value.constraints.len())
                .map(|n| F::from(n as i32))
                .collect(),
        );
        Self {
            target_poly,
            a,
            b,
            c,
        }
    }
}

fn interpolate_lagrange<F: Field>(xs: &[F], ys: &[F]) -> DensePolynomial<F> {
    let mut polynomial = DensePolynomial::from_coefficients_vec(vec![F::zero()]);

    for (i, &xi) in xs.iter().enumerate() {
        let mut buffer_poly = DensePolynomial::from_coefficients_vec(vec![F::one()]);

        for (j, &xj) in xs.iter().enumerate() {
            if i == j {
                continue;
            }

            let denominator = xi - xj;
            let denominator_inv = denominator.inverse().unwrap();
            let numerator_poly = DensePolynomial::from_coefficients_vec(vec![xj.neg(), F::one()]);
            buffer_poly = buffer_poly.naive_mul(&numerator_poly);
            buffer_poly = buffer_poly.naive_mul(&DensePolynomial::from_coefficients_vec(vec![
                denominator_inv,
            ]));
        }

        buffer_poly = buffer_poly.naive_mul(&DensePolynomial::from_coefficients_vec(vec![ys[i]]));
        polynomial = polynomial + buffer_poly;
    }

    polynomial
}

#[cfg(test)]
mod test {
    use ark_ff::Field;
    use ark_poly::Polynomial;
    use super::get_vanishing_polynomial;
    #[test]
    fn test_get_vanishing_polynomial() {
        use ark_bls12_381::Fr;

        let roots = vec![Fr::from(1), Fr::from(2), Fr::from(3)];
        let poly = get_vanishing_polynomial(roots);
        let expected = vec![
            Fr::from(-6i32), // constant term
            Fr::from(11u32), // x
            Fr::from(-6i32), // x^2
            Fr::from(1u32),  // x^3
        ];

        assert_eq!(poly.coeffs, expected);
    }

    #[test]
    fn test_lagrange_interpolation() {
        use super::interpolate_lagrange;
        use ark_bls12_381::Fr;

        let xs = vec![Fr::from(1), Fr::from(2), Fr::from(3)];
        let ys = vec![Fr::from(2), Fr::from(5), Fr::from(10)]; // Correct: x^2 + 1

        let poly = interpolate_lagrange(&xs, &ys);
        let expected = vec![Fr::from(1u32), Fr::from(0u32), Fr::from(1u32)];

        assert_eq!(poly.coeffs, expected);
    }
    fn unit<F: Field>(index: usize, len: usize) -> Vec<F> {
        let mut vec = vec![F::zero(); len];
        vec[index] = F::one();
        vec
    }

    fn add<F: Field>(a: Vec<F>, b: Vec<F>) -> Vec<F> {
        a.into_iter().zip(b).map(|(x, y)| x + y).collect()
    }

    #[test]
    fn test_r1cs_to_qap_structure() {
        use ark_bls12_381::Fr;
        use crate::r1cs::R1CS;
        use super::QAP;

        let mut r1cs = R1CS::<Fr>::new();

        let x = r1cs.add_variable("x".to_string());
        let x_sq = r1cs.add_variable("x_sq".to_string());
        let x_cb = r1cs.add_variable("x_cb".to_string());
        let sym_1 = r1cs.add_variable("sym_1".to_string());
        let out = r1cs.add_variable("~out".to_string());

        // x * x = x_sq
        r1cs.add_constraint(unit(x, r1cs.variables.len()), unit(x, r1cs.variables.len()), unit(x_sq, r1cs.variables.len()));

        // x_sq * x = x_cb
        r1cs.add_constraint(unit(x_sq, r1cs.variables.len()), unit(x, r1cs.variables.len()), unit(x_cb, r1cs.variables.len()));

        // x_cb + x = sym_1
        let a3 = add(unit(x_cb, r1cs.variables.len()), unit(x, r1cs.variables.len()));
        r1cs.add_constraint(a3, unit(0, r1cs.variables.len()), unit(sym_1, r1cs.variables.len()));

        // sym_1 + 5 = out
        let mut a4 = unit(sym_1, r1cs.variables.len());
        a4[0] = Fr::from(5); // constant term
        r1cs.add_constraint(a4, unit(0, r1cs.variables.len()), unit(out, r1cs.variables.len()));

        // Now convert to QAP
        let qap = QAP::from(&r1cs);

        // Check lengths
        assert_eq!(qap.a.len(), r1cs.variables.len());
        assert_eq!(qap.b.len(), r1cs.variables.len());
        assert_eq!(qap.c.len(), r1cs.variables.len());

        // Check target_poly degree (should be 4 for 4 constraints)
        assert_eq!(qap.target_poly.degree(), 4);
    }

}
