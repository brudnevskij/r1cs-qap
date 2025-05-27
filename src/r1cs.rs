use ark_ff::Field;
use std::fmt::{Display, Formatter};

pub struct Constraint<F: Field> {
    pub a: Vec<F>,
    pub b: Vec<F>,
    pub c: Vec<F>,
}

pub enum VariableType {
    Public,
    Intermediate,
    Private,
}

pub struct R1CS<F: Field> {
    pub constraints: Vec<Constraint<F>>,
    pub variables: Vec<String>,
    pub public_variables_count: usize,
}

impl<F: Field> R1CS<F> {
    pub fn new() -> Self {
        Self {
            constraints: vec![],
            variables: vec!["~one".to_string()],
            public_variables_count: 0,
        }
    }

    pub fn add_variable(&mut self, name: String, variable_type: VariableType) -> usize {
        match variable_type {
            VariableType::Public => {
                self.public_variables_count += 1;
                self.variables.insert(self.public_variables_count, name);
                self.variables.len() - 1
            }
            _ => {
                self.variables.push(name);
                self.variables.len() - 1
            }
        }
    }

    pub fn add_constraint(&mut self, a: Vec<F>, b: Vec<F>, c: Vec<F>) {
        self.constraints.push(Constraint { a, b, c });
    }

    pub fn is_satisfied(&self, witness: &[F]) -> bool {
        for constraint in &self.constraints {
            let a = inner_product(&constraint.a, witness);
            let b = inner_product(&constraint.b, witness);
            let c = inner_product(&constraint.c, witness);

            if a * b != c {
                return false;
            }
        }
        true
    }
}

fn inner_product<F: Field>(a: &[F], b: &[F]) -> F {
    a.iter().zip(b).map(|(x, y)| (*x) * (*y)).sum()
}

fn format_side<F: Field>(coefficients: &[F], variables: &Vec<String>) -> String {
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
        write!(f, "{}", r1cs_equations.join("\n"))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::r1cs::VariableType::{Intermediate, Private, Public};
    use ark_bls12_381::Fr;
    use ark_ff::{One, Zero};

    fn cubic_constraint_system() -> R1CS<Fr> {
        // Create constraints for x**3 + x + 5 = 35
        let mut r1cs = R1CS::<Fr>::new();
        r1cs.add_variable("x".to_string(), Public);
        r1cs.add_variable("x_sq".to_string(), Public);
        r1cs.add_variable("x_cb".to_string(), Public);
        r1cs.add_variable("sym_1".to_string(), Public);
        r1cs.add_variable("~out".to_string(), Public);

        // vars = [~one, x, x_sq, x_cb, sym_1, ~out]
        // x*x = x_sq
        let a = vec![
            Fr::zero(),
            Fr::one(),
            Fr::zero(),
            Fr::zero(),
            Fr::zero(),
            Fr::zero(),
        ];
        let b = vec![
            Fr::zero(),
            Fr::one(),
            Fr::zero(),
            Fr::zero(),
            Fr::zero(),
            Fr::zero(),
        ];
        let c = vec![
            Fr::zero(),
            Fr::zero(),
            Fr::one(),
            Fr::zero(),
            Fr::zero(),
            Fr::zero(),
        ];
        r1cs.add_constraint(a, b, c);

        // x_sq * x = x_cb
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
            Fr::one(),
            Fr::zero(),
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

        // x_cb + x = sym_1
        let a = vec![
            Fr::zero(),
            Fr::one(),
            Fr::zero(),
            Fr::one(),
            Fr::zero(),
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
            Fr::one(),
            Fr::zero(),
        ];
        r1cs.add_constraint(a, b, c);

        // sym_1 + 5 = out
        let a = vec![
            Fr::from(5),
            Fr::zero(),
            Fr::zero(),
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
        r1cs
    }

    #[test]
    fn test_r1cs_display() {
        let r1cs = cubic_constraint_system();
        let display_output = format!("{}", r1cs);
        println!("{}", display_output);

        assert!(display_output.contains("(x) * (x) = (x_sq)"));
        assert!(display_output.contains("(x_sq) * (x) = (x_cb)"));
        assert!(display_output.contains("(x + x_cb) * (1) = (sym_1)"));
        assert!(display_output.contains("(5 + sym_1) * (1) = (~out)"));
    }

    #[test]
    fn test_cubic_constraint_system_satisfied() {
        use ark_bls12_381::Fr;
        use ark_ff::Field;

        let r1cs = cubic_constraint_system();

        // witness = [1, 3, 9, 27, 30, 35]
        let witness = vec![
            Fr::from(1u32),
            Fr::from(3u32),
            Fr::from(9u32),
            Fr::from(27u32),
            Fr::from(30u32),
            Fr::from(35u32),
        ];

        assert!(
            r1cs.is_satisfied(&witness),
            "R1CS should be satisfied by this witness"
        );
    }

    #[test]
    fn test_add_variable_public_insert_position() {
        let mut r1cs = R1CS::<Fr>::new();

        let idx_x = r1cs.add_variable("x".to_string(), Public);
        let idx_y = r1cs.add_variable("y".to_string(), Public);

        // Expected: ["~one", "x", "y"]
        assert_eq!(r1cs.variables[idx_x], "x");
        assert_eq!(r1cs.variables[idx_y], "y");
        assert_eq!(idx_x, 1);
        assert_eq!(idx_y, 2);
        assert_eq!(r1cs.public_variables_count, 2);
    }

    #[test]
    fn test_add_variable_private_append_position() {
        let mut r1cs = R1CS::<Fr>::new();

        r1cs.add_variable("x".to_string(), Public); // idx 1
        r1cs.add_variable("y".to_string(), Public); // idx 2
        let idx_z = r1cs.add_variable("z".to_string(), Private); // should be 3

        assert_eq!(r1cs.variables[idx_z], "z");
        assert_eq!(idx_z, 3);
        assert_eq!(r1cs.variables.len(), 4);
    }

    #[test]
    fn test_add_variable_intermediate_append_position() {
        let mut r1cs = R1CS::<Fr>::new();

        r1cs.add_variable("a".to_string(), Public); // idx 1
        let idx_tmp = r1cs.add_variable("tmp".to_string(), Intermediate); // should be 2

        assert_eq!(r1cs.variables[idx_tmp], "tmp");
        assert_eq!(idx_tmp, 2);
    }

    #[test]
    fn test_public_variables_count_tracking() {
        let mut r1cs = R1CS::<Fr>::new();
        assert_eq!(r1cs.public_variables_count, 0);

        r1cs.add_variable("x".to_string(), Public);
        assert_eq!(r1cs.public_variables_count, 1);

        r1cs.add_variable("y".to_string(), Private);
        assert_eq!(r1cs.public_variables_count, 1);

        r1cs.add_variable("z".to_string(), Public);
        assert_eq!(r1cs.public_variables_count, 2);
    }
}
