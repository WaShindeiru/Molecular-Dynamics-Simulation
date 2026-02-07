use nalgebra::Vector3;

pub fn cos_from_vec(r_ij: &Vector3<f64>, r_ik: &Vector3<f64>) -> f64 {
    let rij_length = r_ij.magnitude();
    let rik_length = r_ik.magnitude();

    if rij_length == 0. || rik_length == 0. {
        return 0.
    }

    let result = r_ij.dot(r_ik) / (rij_length * rik_length);

    assert!(result >= -1. - 1e-12 && result <= 1. + 1e-12);

    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_perpendicular() {
        let a = Vector3::new(0., 4., 0.);
        let b = Vector3::new(0., 0., 10.);

        assert!(cos_from_vec(&a, &b).abs() == 0.);
    }

    #[test]
    fn test_45() {
        let a = Vector3::new(4., 0., 0.);
        let b = Vector3::new(4., 4., 0.);

        assert!((cos_from_vec(&a, &b) - (2.0_f64.sqrt() / 2.)).abs() < 1e-14);
    }
}