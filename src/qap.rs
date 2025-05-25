use ark_ff::Field;
use ark_poly::DenseUVPolynomial;
use ark_poly::univariate::DensePolynomial;
use ark_poly::polynomial::Polynomial;
use ark_std::ops::Neg;

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


#[cfg(test)]
mod test {
    use super::get_vanishing_polynomial;
    #[test]
    fn test_get_vanishing_polynomial() {
        use ark_bls12_381::Fr;

        let roots = vec![Fr::from(1), Fr::from(2), Fr::from(3)];
        let poly = get_vanishing_polynomial(roots);
        let expected = vec![
            Fr::from(-6i32),  // constant term
            Fr::from(11u32),  // x
            Fr::from(-6i32),  // x^2
            Fr::from(1u32),   // x^3
        ];

        assert_eq!(poly.coeffs, expected);
    }
}