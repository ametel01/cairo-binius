use alexandria_math::BitShift;

mod binary_field;
pub mod utils;
mod binary_ntt;

pub trait BitLenTrait<T> {
    fn bit_len(self: T) -> T;
}

impl BitLenImpl of BitLenTrait<u128> {
    fn bit_len(self: u128) -> u128 {
        let mut lenght = 0;
        if self == 0 {
            return 0;
        } else {
            let mut n = self;
            while (n > 0) {
                lenght += 1;
                n = BitShift::shr(n, 1);
            }
        }
        lenght
    }
}
