pub mod arithmetic;
pub mod cell;
pub mod debug;
pub mod delaunay;
pub mod hall_symbol;
pub mod kgrid;
pub mod kpoint;
pub mod mathfunc;
pub mod niggli;
pub mod overlap;
pub mod pointgroup;
pub mod primitive;
pub mod sitesym_database;
pub mod spacegroup;
pub mod spg_database;
pub mod spin;
pub mod symmetry;
pub fn add(left: u64, right: u64) -> u64 {
    left + right
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_works() {
        let result = add(2, 2);
        assert_eq!(result, 4);
    }
}
