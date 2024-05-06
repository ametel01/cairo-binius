use core::array::ArrayTrait;
use core::traits::Into;
use core::option::OptionTrait;
use core::traits::TryInto;
use alexandria_data_structures::array_ext::SpanTraitExt;
use core::array::SpanTrait;
use binius::binary_field::{BinaryFieldElement as B, BnFOpsTrait, bnf};
use binius::utils::{log2, PolyOps};
use alexandria_math::BitShift;
use alexandria_data_structures::vec::{VecTrait, Felt252Vec};

fn get_Wi(i: u128) -> Array<B> {
    if i == 0 {
        return array![bnf(0), bnf(1)];
    } else {
        let mut prev = get_Wi(i - 1);
        let mut b = array![bnf(1)];
        let mut o: Array<B> = PolyOps::mul_polys(prev.clone(), PolyOps::add_polys(ref prev, ref b));
        let inv_quot = PolyOps::eval_poly_at(o.span(), bnf(BitShift::shl(1, i))).inv();
        let mut out: Array<B> = array![];
        let len = o.len();
        let mut j = 0;
        while j < len {
            out.append(*o.at(j) * inv_quot);
            j += 1;
        };
        out
    }
}

fn get_Wi_eval(dim: usize, pt: u128, ref cache: Felt252Vec<u128>) -> B {
    // let mut Wi_eval_cache = VecTrait::<Felt252Vec, u128>::new();
    let mut i = 0;
    while i <= dim {
        cache.push(0);
        i += 1;
    };
    if dim == 0 {
        return bnf(pt);
    }
    let mut new_value = Default::default();
    if pt != cache.at(dim) {
        let prev = get_Wi_eval(dim - 1, pt, ref cache);
        let prev_quot = get_Wi_eval(dim - 1, BitShift::shl(1, dim.into()), ref cache);
        new_value = prev * (prev + bnf(1)) / (prev_quot * (prev_quot + bnf(1)));
        cache.set(dim, new_value.value);
    }
    new_value
}

fn additive_ntt(vals: Span<B>, start: u128, ref cache: Felt252Vec<u128>) -> Array<B> {
    let vals_len = vals.len();
    if vals_len == 1 {
        return vals.dedup();
    }
    let halflen = vals_len / 2;
    let mut L = array![];
    let mut R = array![];
    let mut i = 0;
    while i < halflen {
        L.append(*vals.at(i));
        R.append(*vals.at(i + halflen));
        i += 1;
    };
    let coeff1 = get_Wi_eval(log2(halflen.into()).try_into().unwrap(), start, ref cache);
    let mut sub_input1 = array![];
    let mut sub_input2 = array![];
    let mut j = 0;
    while j < halflen {
        sub_input1.append(*L.at(j) + *R.at(j) * coeff1);
        sub_input2.append(*sub_input1.at(j) - *R.at(j));
        j += 1;
    };
    sub_input1 = additive_ntt(sub_input1.span(), start, ref cache);
    sub_input2 = additive_ntt(sub_input2.span(), start + halflen.into(), ref cache);
    let mut o = array![];
    let mut k = 0;
    let len = sub_input1.len();
    while k < len {
        o.append(*sub_input1.at(k));
        k += 1;
    };
    let mut l = 0;
    let len = sub_input2.len();
    while l < len {
        o.append(*sub_input2.at(l));
        l += 1;
    };
    o
}

fn inv_additive_ntt(vals: Span<B>, start: u128, ref cache: Felt252Vec<u128>) -> Array<B> {
    let vals = vals.dedup();
    let vals_len = vals.len();
    if vals_len == 1 {
        return vals;
    }
    let halflen = vals_len / 2;
    let mut L = array![];
    let mut R = array![];
    let mut i = 0;
    while i < halflen {
        L.append(*vals.at(i));
        R.append(*vals.at(i + halflen));
        i += 1;
    };
    L = inv_additive_ntt(L.span(), start, ref cache);
    R = inv_additive_ntt(R.span(), start + halflen.into(), ref cache);
    let coeff1 = get_Wi_eval(log2(halflen.into()).try_into().unwrap(), start, ref cache);
    let coeff2 = coeff1 + bnf(1);
    let mut o = array![];
    let mut j = 0;
    let L_cp = L.clone();
    let R_cp = R.clone();
    while j < halflen {
        o.append(*L.at(j) * coeff2 + *R.at(j) * coeff1);
        j += 1;
    };
    let mut k = 0;
    while k < halflen {
        o.append(*L_cp.at(k) + *R_cp.at(k));
        k += 1;
    };
    o
}

fn extend(data: Span<B>, expansion_factor: usize) -> Array<B> {
    let mut cache = VecTrait::<Felt252Vec, u128>::new();
    let mut t0 = inv_additive_ntt(data, 0, ref cache);
    // let mut t1 = array![];
    let len_data = data.len();
    let mut i = 0;
    while i < len_data * (expansion_factor - 1) {
        t0.append(bnf(0));
        i += 1;
    };
    additive_ntt(t0.span(), 0, ref cache)
}

#[cfg(test)]
mod test {
    use super::{get_Wi, bnf, VecTrait, Felt252Vec, additive_ntt, inv_additive_ntt, extend};

    #[test]
    fn get_Wi_test() {
        let res = get_Wi(0);
        assert!(res == array![bnf(0), bnf(1)], "got {:?}", res);

        let res = get_Wi(1);
        assert!(res == array![bnf(0), bnf(1), bnf(1)], "got {:?}", res);

        let res = get_Wi(2);
        assert!(res == array![bnf(0), bnf(3), bnf(0), bnf(0), bnf(3)], "got {:?}", res);
    }

    #[test]
    fn additive_ntt_test() {
        let vals = array![bnf(1), bnf(2), bnf(3), bnf(4)];
        let mut cache = VecTrait::<Felt252Vec, u128>::new();
        let res = super::additive_ntt(vals.span(), 0, ref cache);
        assert!(res == array![bnf(1), bnf(3), bnf(9), bnf(15)], "got {:?}", res);
    }

    #[test]
    fn inv_additive_ntt_test() {
        let vals = array![bnf(1), bnf(3), bnf(9), bnf(15)];
        let mut cache = VecTrait::<Felt252Vec, u128>::new();
        let res = super::inv_additive_ntt(vals.span(), 0, ref cache);
        assert!(res == array![bnf(1), bnf(2), bnf(3), bnf(4)], "got {:?}", res);
    }

    #[test]
    fn extend_test() {
        let vals = array![bnf(1), bnf(3), bnf(9), bnf(15)];
        let res = extend(vals.span(), 2);
        assert!(
            res == array![bnf(1), bnf(3), bnf(9), bnf(15), bnf(14), bnf(15), bnf(14), bnf(11)],
            "got {:?}",
            res
        );
    }
}
