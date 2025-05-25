use ark_ff::Field;
use std::fmt::{Display, Formatter};

pub struct Constraint<F: Field> {
    pub a: Vec<F>,
    pub b: Vec<F>,
    pub c: Vec<F>,
}

struct R1CS<F: Field> {
    pub constraints: Vec<Constraint<F>>,
    pub variables: Vec<String>,
}

impl<F: Field> R1CS<F> {
    pub fn new() -> Self {
        Self {
            constraints: vec![],
            variables: vec!["~one".to_string()],
        }
    }

    pub fn add_variable(&mut self, name: String) -> usize {
        self.variables.push(name);
        self.variables.len() - 1
    }

    pub fn add_constraint(&mut self, a: Vec<F>, b: Vec<F>, c: Vec<F>) {
        self.constraints.push(Constraint { a, b, c });
    }
}

fn format_side<F:Field>(coefficients: &[F], variables: &Vec<String>) -> String {
    let mut side = vec![];

    if !coefficients[0].is_zero() {
        side.push(coefficients[0].to_string());
    }

    for (i, var) in variables.iter().enumerate().skip(1) {
        let coeff = &coefficients[i];
        if !coeff.is_zero() {
            if coeff == &F::ONE {
                side.push(var.clone());
            } else {
                side.push(format!("{}*{}", coeff, var));
            }
        }
    }
    format!("({})", side.join(" + "))
}

impl<F: Field> Display for R1CS<F> {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        let mut r1cs_equations: Vec<String> = vec![];
        for constraint in &self.constraints {
            let a = format_side(&constraint.a, &self.variables);
            let b = format_side(&constraint.b, &self.variables);
            let c = format_side(&constraint.c, &self.variables);
            r1cs_equations.push(format!("{} * {} = {}", a, b, c));
        }
        write!(f,"{}", r1cs_equations.join("\n"))
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use ark_ff::{Zero, One};
    use ark_bls12_381::Fr;

    #[test]
    fn test_r1cs_display() {
        // Create constraints for x**3 + x + 5 = 35
        let mut r1cs = R1CS::<Fr>::new();
        r1cs.add_variable("x".to_string());
        r1cs.add_variable("x_sq".to_string());
        r1cs.add_variable("x_cb".to_string());
        r1cs.add_variable("sym_1".to_string());
        r1cs.add_variable("~out".to_string());


        // vars = [~one, x, x_sq, x_cb, sym_1, ~out]
        // x*x = x_sq
        let  a = vec![Fr::zero(), Fr::one(), Fr::zero(), Fr::zero(), Fr::zero(), Fr::zero()];
        let  b = vec![Fr::zero(), Fr::one(), Fr::zero(), Fr::zero(), Fr::zero(), Fr::zero()];
        let  c = vec![Fr::zero(), Fr::zero(), Fr::one(), Fr::zero(), Fr::zero(), Fr::zero()];
        r1cs.add_constraint(a, b, c);

        // x_sq * x = x_cb
        let  a = vec![Fr::zero(), Fr::zero(), Fr::one(), Fr::zero(), Fr::zero(), Fr::zero()];
        let  b = vec![Fr::zero(), Fr::one(), Fr::zero(), Fr::zero(), Fr::zero(), Fr::zero()];
        let  c = vec![Fr::zero(), Fr::zero(), Fr::zero(), Fr::one(), Fr::zero(), Fr::zero()];
        r1cs.add_constraint(a, b, c);

        // x_cb + x = sym_1
        let  a = vec![Fr::zero(), Fr::one(), Fr::zero(), Fr::one(), Fr::zero(), Fr::zero()];
        let  b = vec![Fr::one(), Fr::zero(), Fr::zero(), Fr::zero(), Fr::zero(), Fr::zero()];
        let  c = vec![Fr::zero(), Fr::zero(), Fr::zero(), Fr::zero(), Fr::one(), Fr::zero()];
        r1cs.add_constraint(a, b, c);

        // sym_1 + 5 = out
        let  a = vec![Fr::from(5), Fr::zero(), Fr::zero(), Fr::zero(), Fr::one(), Fr::zero()];
        let  b = vec![Fr::one(), Fr::zero(), Fr::zero(), Fr::zero(), Fr::zero(), Fr::zero()];
        let  c = vec![Fr::zero(), Fr::zero(), Fr::zero(), Fr::zero(), Fr::zero(), Fr::one()];
        r1cs.add_constraint(a, b, c);

        let display_output = format!("{}", r1cs);
        println!("{}", display_output);

        assert!(display_output.contains("(x) * (x) = (x_sq)"));
        assert!(display_output.contains("(x_sq) * (x) = (x_cb)"));
        assert!(display_output.contains("(x + x_cb) * (1) = (sym_1)"));
        assert!(display_output.contains("(5 + sym_1) * (1) = (~out)"));
    }
}
