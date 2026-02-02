fn main() {
    let a = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    let b = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

    println!("Testing mat_get_similar_matrix_d3 with identity matrices...");
    println!("a = {:?}", a);
    println!("b = {:?}", b);

    // We need to use the actual function from the library
    // For now, let's manually compute
    // Similarity transform: b^-1 * a * b
    // For identity matrices, this should be identity

    let det = a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
        - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
        + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);
    println!("det(a) = {}", det);

    if det.abs() < 0.0 {
        println!("det.abs() < 0.0 is {}", det.abs() < 0.0);
    }
}

