use core::num::traits::zero::Zero;
use core::traits::{BitAnd, BitXor, BitOr};
use core::integer::u128_wide_mul;
use alexandria_math::{BitShift, pow};
use binius::BitLenTrait;
use binius::utils::max;

#[derive(Copy, Default, Debug, Drop)]
pub struct BinaryFieldElement {
    pub value: u128
}

pub fn bnf(value: u128) -> BinaryFieldElement {
    BinaryFieldElement { value }
}

#[generate_trait]
pub impl BnFOpsImpl of BnFOpsTrait {
    fn pow(self: BinaryFieldElement, n: u128) -> BinaryFieldElement {
        if n == 0 {
            return bnf(1);
        } else if n == 1 {
            return self;
        } else if n == 2 {
            return self * self;
        } else {
            return self.pow(n % 2) * self.pow(n / 2).pow(2);
        }
    }

    fn inv(self: BinaryFieldElement) -> BinaryFieldElement {
        let L = BitShift::shl(1, (self.value.bit_len() - 1).bit_len());
        return self.pow(pow(2, L) - 2);
    }
}

impl BnFAdd of Add<BinaryFieldElement> {
    fn add(lhs: BinaryFieldElement, rhs: BinaryFieldElement) -> BinaryFieldElement {
        bnf(BitXor::bitxor(lhs.value, rhs.value))
    }
}

impl BnFSub of Sub<BinaryFieldElement> {
    fn sub(lhs: BinaryFieldElement, rhs: BinaryFieldElement) -> BinaryFieldElement {
        lhs + rhs
    }
}

impl BnFMul of Mul<BinaryFieldElement> {
    fn mul(lhs: BinaryFieldElement, rhs: BinaryFieldElement) -> BinaryFieldElement {
        bnf(binmul(lhs.value, rhs.value, Option::None))
    }
}

impl BnFDiv of Div<BinaryFieldElement> {
    fn div(lhs: BinaryFieldElement, rhs: BinaryFieldElement) -> BinaryFieldElement {
        lhs * rhs.inv()
    }
}

pub impl U128BitShift of BitShift<BinaryFieldElement> {
    fn shl(x: BinaryFieldElement, n: BinaryFieldElement) -> BinaryFieldElement {
        let (_, bottom_word) = u128_wide_mul(x.value, pow(2, n.value));
        bnf(bottom_word)
    }

    fn shr(x: BinaryFieldElement, n: BinaryFieldElement) -> BinaryFieldElement {
        bnf(x.value / pow(2, n.value))
    }
}

impl BnFBitLenImpl of BitLenTrait<BinaryFieldElement> {
    fn bit_len(self: BinaryFieldElement) -> BinaryFieldElement {
        bnf(BitShift::shl(1, (self.value.bit_len() - 1).bit_len()))
    }
}

impl PartialEqBnf of PartialEq<BinaryFieldElement> {
    fn eq(lhs: @BinaryFieldElement, rhs: @BinaryFieldElement) -> bool {
        lhs.value == rhs.value
    }
    fn ne(lhs: @BinaryFieldElement, rhs: @BinaryFieldElement) -> bool {
        lhs.value != rhs.value
    }
}


pub fn binmul(lhs: u128, rhs: u128, len: Option<u128>) -> u128 {
    if lhs == 0 || rhs == 0 {
        return 0;
    }
    if lhs < 2 || rhs < 2 {
        return lhs * rhs;
    }
    let mut length = 0;
    match len {
        Option::None => { length = BitShift::shl(1, (max(lhs, rhs).bit_len() - 1).bit_len()) },
        Option::Some(len) => { length = len; }
    }
    let halflen = length / 2;
    let quarterlen = length / 4;
    let halfmask = BitShift::shl(1, halflen) - 1;

    let L1 = BitAnd::bitand(lhs, halfmask);
    let R1 = BitShift::shr(lhs, halflen);
    let L2 = BitAnd::bitand(rhs, halfmask);
    let R2 = BitShift::shr(rhs, halflen);

    if L1 == 0 && R1 == 1 {
        let outR = BitXor::bitxor(
            binmul(BitShift::shl(1, quarterlen), R2, Option::Some(halflen)), L2
        );
        return BitXor::bitxor(R2, BitShift::shl(outR, halflen));
    }

    let L1L2 = binmul(L1, L2, Option::Some(halflen));
    let R1R2 = binmul(R1, R2, Option::Some(halflen));
    let R1R2_high = binmul(BitShift::shl(1, quarterlen), R1R2, Option::Some(halflen));
    let Z3 = binmul(BitXor::bitxor(L1, R1), BitXor::bitxor(L2, R2), Option::Some(halflen));

    let t0 = BitXor::bitxor(Z3, L1L2);
    let t1 = BitXor::bitxor(t0, R1R2);
    let t2 = BitXor::bitxor(t1, R1R2_high);
    let t3 = BitShift::shl(t2, halflen);
    let t4 = BitXor::bitxor(L1L2, R1R2);

    return BitXor::bitxor(t4, t3);
}

#[cfg(test)]
mod test {
    use super::{bnf, BinaryFieldElement, BitLenTrait, binmul, BnFOpsTrait};

    #[test]
    fn add_test() {
        let res = bnf(16) + bnf(251);
        assert!(res.value == 235, "16 + 251 = {}", res.value)
    }

    #[test]
    fn bit_len_test() {
        let res = bnf(1000).bit_len();
        assert!(res.value == 16, "bit_len(1000) = {}", res.value);
    }

    #[test]
    fn binmul_test() {
        let res = binmul(4, 15, Option::None);
        assert!(res == 11, "4 * 15 = {}", res);
    }

    #[test]
    fn bnf_mul_test() {
        let res = bnf(23) * bnf(17);
        assert!(res.value == 38, "23 * 17 = {}", res.value)
    }

    #[test]
    fn bnf_pow_test() {
        let res = bnf(16).pow(251);
        assert!(res.value == 174, "16 ** 251 = {}", res.value);
    }

    #[test]
    fn bnf_inv_test() {
        let res = bnf(23).inv();
        assert!(res.value == 126, "23 ** -1 = {}", res.value);
    }

    #[test]
    fn bnf_div_test() {
        let res = bnf(23) / bnf(11);
        assert!(res.value == 207, "23 / 11 = {}", res.value);
    }

    #[test]
    fn test_biniry_operations() {
        let t0 = (bnf(3) + bnf(14)) * bnf(15);
        let t1 = bnf(3) * bnf(15) + bnf(14) * bnf(15);
        assert!(t0.value == t1.value, "t0: {:?}, t1: {:?}", t0, t1);

        let res = bnf(13).pow(255);
        assert!(res.value == 1, "13 ** 255 = {}", res.value);

        let res = bnf(420).pow(255);
        assert!(res.value != 1, "420 ** 255 = {}", res.value);

        let res = bnf(420).pow(65535);
        assert!(res.value == 1, "420 ** 65535 = {}", res.value);
    }
}
