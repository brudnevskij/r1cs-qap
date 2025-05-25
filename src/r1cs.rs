use ark_ff::Field;
use std::fmt::{Display, Formatter};

pub struct Constraint<F: Field> {
    pub a: Vec<F>,
    pub b: Vec<F>,
    pub c: Vec<F>,
}

struct R1CS<F: Field> {
    pub constraints: Vec<Constraint<F>>,
    pub num_variables: usize,
    pub variables: Vec<String>,
}

impl<F: Field> R1CS<F> {
    pub fn new() -> Self {
        Self {
            constraints: vec![],
            num_variables: 0,
            variables: vec!["~one".to_string()],
        }
    }

    pub fn add_variable(&mut self, name: String) -> usize {
        self.variables.push(name);
        self.num_variables += 1;
        self.num_variables
    }

    pub fn add_constraint(&mut self, a: Vec<F>, b: Vec<F>, c: Vec<F>) {
        self.constraints.push(Constraint { a, b, c });
    }
}

impl<F: Field> Display for R1CS<F> {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        let mut r1cs_equations: Vec<String> = vec![];
        for constraint in &self.constraints {
            let mut a_buffer = String::from("(");
            let mut b_buffer = String::from("(");
            let mut c_buffer = String::from("(");

            // adding constants if present
            if constraint.a[0] != F::ZERO {
                a_buffer.push_str(&&constraint.a[0].to_string());
            }
            if constraint.b[0] != F::ZERO {
                b_buffer.push_str(&constraint.b[0].to_string());
            }
            if constraint.c[0] != F::ZERO {
                c_buffer.push_str(&constraint.c[0].to_string());
            }
            // handling rest of the variables
            for (i, variable) in self.variables.iter().enumerate().skip(1){
                // A
                let mut a_side = String::new();
                let coeff:F = constraint.a[i];

                if !coeff.is_zero() {
                    a_side = if coeff == F::ONE {
                        variable.clone()
                    } else {
                        format!("{}*{}", coeff, variable)
                    };

                    if a_buffer == "(" {
                        a_buffer.push_str(&a_side);
                    } else {
                        a_buffer.push_str(&format!(" + {}", a_side));
                    }
                }

                // B
                let mut b_side = String::new();
                let coeff: F = constraint.b[i];

                if !coeff.is_zero() {
                    b_side = if coeff == F::ONE {
                        variable.clone()
                    } else {
                        format!("{}*{}", coeff, variable)
                    };

                    if b_buffer == "(" {
                        b_buffer.push_str(&b_side);
                    } else {
                        b_buffer.push_str(&format!(" + {}", b_side));
                    }
                }

                // C
                let mut c_side = String::new();
                let coeff: F = constraint.c[i];

                if !coeff.is_zero() {
                    c_side = if coeff == F::ONE {
                        variable.clone()
                    } else {
                        format!("{}*{}", coeff, variable)
                    };

                    if c_buffer == "(" {
                        c_buffer.push_str(&c_side);
                    } else {
                        c_buffer.push_str(&format!(" + {}", c_side));
                    }
                }
            }

            r1cs_equations.push(format!("{}) * {}) = {})", a_buffer, b_buffer, c_buffer));
        }
        write!(f,"{}", r1cs_equations.join("\n"))
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use ark_ff::{Field, Zero, One};
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
