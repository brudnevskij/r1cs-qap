use crate::r1cs::R1CS;
use ark_ff::Field;
use ark_poly::DenseUVPolynomial;
use ark_poly::Polynomial;
use ark_poly::univariate::DensePolynomial;


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
        let nr = root.neg();
        let factor_polynomial = DensePolynomial::from_coefficients_vec(vec![nr, F::one()]);
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

        let target_poly =
            get_vanishing_polynomial((1..value.constraints.len()).map(|n| F::from(n as i32)).collect());
        Self {
            target_poly,
            a,
            b,
            c,
        }
    }
}

fn interpolate_lagrange<F: Field>(xs: &[F], ys: &[F]) -> DensePolynomial<F> {
    DensePolynomial::from_coefficients_vec(vec![F::one()])
}
#[cfg(test)]
mod test {
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
}
