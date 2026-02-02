fn main() {
    // Test is_identity_metric logic
    let metric_orig: [[f64; 3]; 3] = [[2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 2.0]]; // For lattice vectors of length sqrt(2)
    let metric_rot: [[f64; 3]; 3] = [[2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 2.0]];

    let symprec = 1e-5;
    let angle_symprec = -1.0;

    // Manual check
    let mut length_orig = [0.0; 3];
    let mut length_rot = [0.0; 3];
    for i in 0..3 {
        length_orig[i] = metric_orig[i][i].sqrt();
        length_rot[i] = metric_rot[i][i].sqrt();
        let diff = (length_orig[i] - length_rot[i]).abs();
        println!(
            "Length[{}]: orig={}, rot={}, diff={}, symprec={}, diff>symprec? {}",
            i,
            length_orig[i],
            length_rot[i],
            diff,
            symprec,
            diff > symprec
        );
    }

    // Check angles
    let elem_sets = [[0, 1], [0, 2], [1, 2]];
    for i in 0..3 {
        let j = elem_sets[i][0];
        let k = elem_sets[i][1];
        if angle_symprec > 0.0 {
            // Not used
        } else {
            let cos1 = metric_orig[j][k] / length_orig[j] / length_orig[k];
            let cos2 = metric_rot[j][k] / length_rot[j] / length_rot[k];
            let x = cos1 * cos2 + (1.0 - cos1 * cos1).sqrt() * (1.0 - cos2 * cos2).sqrt();
            let sin_dtheta2 = 1.0 - x * x;
            let length_ave2 =
                ((length_orig[j] + length_rot[j]) * (length_orig[k] + length_rot[k])) / 4.0;
            println!(
                "Angle check {}{}: cos1={}, cos2={}, x={}, sin_dtheta2={}, length_ave2={}",
                j, k, cos1, cos2, x, sin_dtheta2, length_ave2
            );
            println!(
                "  sin_dtheta2 > 1e-12? {}, sin_dtheta2 * length_ave2 > 1e-10? {}",
                sin_dtheta2 > 1e-12,
                sin_dtheta2 * length_ave2 > 1e-10
            );
        }
    }
}
