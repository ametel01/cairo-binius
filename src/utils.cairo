use core::clone::Clone;
use core::option::OptionTrait;
use core::array::ArrayTrait;
use alexandria_math::pow;
use binius::BitLenTrait;
use binius::binary_field::{BinaryFieldElement as B, bnf};

fn multilinear_poly_eval(ref evals: Array<u128>, ref pt: Array<u128>) -> u128 {
    let evals_len = evals.len();
    let pt_len = pt.len();
    assert!(evals_len == pow(2, pt_len), "evals.len() must be 2^pt.len()");
    if pt_len == 0 {
        return evals.pop_front().unwrap();
    }
    let mut top_evals = array![];
    let mut i = 0;
    while i < evals_len / 2 {
        top_evals.append(evals.pop_front().unwrap());
        i += 1;
    };
    let mut bottom_evals = array![];
    i = 0;
    while i < evals_len / 2 {
        bottom_evals.append(evals.pop_front().unwrap());
        i += 1;
    };
    let mut top_point = array![];
    let mut i = 0;
    let point_cp = pt.clone();
    while i < pt_len - 1 {
        top_point.append(*point_cp.at(i));
        i += 1;
    };
    let mut top = multilinear_poly_eval(ref top_evals, ref top_point);
    let mut bottom = multilinear_poly_eval(ref bottom_evals, ref top_point);

    let pt_last = *pt.at(pt_len - 1);
    return (bottom - top) * pt_last + top;
}

pub fn log2(x: u128) -> u128 {
    assert!(BitAnd::bitand(x, x - 1) == 0, "x must be a power of 2");
    x.bit_len() - 1
}

pub trait PolyOps<T> {
    fn add_polys(ref poly1: Array<T>, ref poly2: Array<T>) -> Array<T>;
    fn mul_polys(poly1: Array<T>, poly2: Array<T>) -> Array<T>;
    fn eval_poly_at(poly: Span<T>, pt: T) -> T;
}

impl PolyOptU128 of PolyOps<u128> {
    fn add_polys(ref poly1: Array<u128>, ref poly2: Array<u128>) -> Array<u128> {
        let len_a = poly1.len();
        let len_b = poly2.len();
        let minl = if len_a < len_b {
            len_a
        } else {
            len_b
        };
        let maxl = if len_a > len_b {
            len_a
        } else {
            len_b
        };
        let mut i = 0;
        while i < maxl - minl {
            poly1.append(0);
            poly2.append(0);
            i += 1;
        };
        let mut o = array![];
        let mut i = 0;
        while i < maxl {
            o.append(poly1.pop_front().unwrap() + poly2.pop_front().unwrap());
            i += 1;
        };
        o
    }

    fn mul_polys(poly1: Array<u128>, poly2: Array<u128>) -> Array<u128> {
        let len_a = poly1.len();
        let len_b = poly2.len();
        let mut o = array![];
        let mut i = 0;
        while i < len_a + len_b - 1 {
            o.append(0);
            i += 1;
        };
        let mut i = 0;
        while i < len_a {
            let aval = *poly1.at(i);
            if aval != 0 {
                let mut temp: Array<u128> = array![];
                let mut j = 0;
                while j < i {
                    temp.append(0);
                    j += 1;
                };
                let mut j = 0;
                let poly2_cp = poly2.clone();
                while j < len_b {
                    temp.append(aval * *poly2_cp.at(j));
                    j += 1;
                };
                o = PolyOps::add_polys(ref o, ref temp);
            }
            i += 1;
        };
        o
    }

    fn eval_poly_at(poly: Span<u128>, pt: u128) -> u128 {
        let mut o = 0;
        let mut power = 1;
        let mut i = 0;
        let len = poly.len();
        while i < len {
            o += *poly.at(i) * power;
            power *= pt;
            i += 1;
        };
        o
    }
}

impl PolyOptBinaryField of PolyOps<B> {
    fn add_polys(ref poly1: Array<B>, ref poly2: Array<B>) -> Array<B> {
        let len_a = poly1.len();
        let len_b = poly2.len();
        let minl = if len_a < len_b {
            len_a
        } else {
            len_b
        };
        let maxl = if len_a > len_b {
            len_a
        } else {
            len_b
        };
        let mut i = 0;
        while i < maxl - minl {
            poly1.append(bnf(0));
            poly2.append(bnf(0));
            i += 1;
        };
        let mut o = array![];
        let mut i = 0;
        while i < maxl {
            o.append(poly1.pop_front().unwrap() + poly2.pop_front().unwrap());
            i += 1;
        };
        o
    }

    fn mul_polys(poly1: Array<B>, poly2: Array<B>) -> Array<B> {
        let len_a = poly1.len();
        let len_b = poly2.len();
        let mut o = array![];
        let mut i = 0;
        while i < len_a + len_b - 1 {
            o.append(bnf(0));
            i += 1;
        };
        let mut i = 0;
        while i < len_a {
            let aval = *poly1.at(i);
            if aval.value != 0 {
                let mut temp: Array<B> = array![];
                let mut j = 0;
                while j < i {
                    temp.append(bnf(0));
                    j += 1;
                };
                let mut j = 0;
                let poly2_cp = poly2.clone();
                while j < len_b {
                    temp.append(aval * *poly2_cp.at(j));
                    j += 1;
                };
                o = PolyOps::add_polys(ref o, ref temp);
            }
            i += 1;
        };
        o
    }

    fn eval_poly_at(poly: Span<B>, pt: B) -> B {
        let mut o = bnf(0);
        let mut power = bnf(1);
        let mut i = 0;
        let len = poly.len();
        while i < len {
            o = o + *poly.at(i) * power;
            power = power * pt;
            i += 1;
        };
        o
    }
}

fn evaluation_tensor_product(pt: Span<B>) -> Array<B> {
    let mut o = array![bnf(1)];
    let len = pt.len();
    let mut i = 0;
    while i < len {
        let o_cp = o.clone();
        let mut new_o = array![];
        let len_o = o.len();
        let mut j = 0;
        while j < len_o {
            new_o.append((bnf(1) - *pt.at(i)) * *o.at(j));
            j += 1;
        };
        let mut k = 0;
        while k < len_o {
            new_o.append(*pt.at(i) * *o_cp.at(k));
            k += 1;
        };
        o = new_o;
        i += 1;
    };
    o
}

pub fn max(a: u128, b: u128) -> u128 {
    if a > b {
        return a;
    }
    b
}

#[cfg(test)]
mod test {
    use super::{multilinear_poly_eval, log2, PolyOps, evaluation_tensor_product, bnf};

    #[test]
    fn test_multilinear_poly_eval() {
        let mut evals = array![3, 14, 15, 92];
        let mut pt = array![0, 0];
        let res = multilinear_poly_eval(ref evals, ref pt);
        assert!(res == 3, "expected 3, got {}", res);

        let mut evals = array![3, 14, 15, 92];
        let mut pt = array![1, 0];
        let res = multilinear_poly_eval(ref evals, ref pt);
        assert!(res == 14, "expected 14, got {}", res);

        let mut evals = array![3, 14, 15, 92];
        let mut pt = array![2, 5];
        let res = multilinear_poly_eval(ref evals, ref pt);
        assert!(res == 745, "expected 745, got {}", res);
    }

    #[test]
    fn test_log2() {
        let res = log2(1);
        assert!(res == 0, "expected 0, got {}", res);

        let res = log2(2);
        assert!(res == 1, "expected 1, got {}", res);

        let res = log2(4);
        assert!(res == 2, "expected 2, got {}", res);

        let res = log2(8);
        assert!(res == 3, "expected 3, got {}", res);

        let res = log2(16);
        assert!(res == 4, "expected 4, got {}", res);
    }

    #[test]
    fn test_eval_poly_at() {
        let mut poly = array![3_u128, 1, 4, 1, 5];
        let res = PolyOps::eval_poly_at(poly.span(), 10);
        assert!(res == 51413, "expected 51413, got {}", res);
    }

    #[test]
    fn test_add_polys() {
        let mut poly1 = array![3_u128, 1, 4, 1, 5];
        let mut poly2 = array![9, 2, 6, 5, 3];
        let res = PolyOps::add_polys(ref poly1, ref poly2);
        assert!(res == array![12, 3, 10, 6, 8], "expected [12, 3, 10, 6, 8], got {:?}", res);

        let mut poly1 = array![3_u128, 1, 4, 1, 5];
        let mut poly2 = array![9, 2, 6];
        let res = PolyOps::add_polys(ref poly1, ref poly2);
        assert!(res == array![12, 3, 10, 1, 5], "expected [12, 3, 10, 6, 8], got {:?}", res);
    }

    #[test]
    fn test_mul_polys() {
        let mut poly1 = array![3_u128, 1, 4, 1, 5];
        let mut poly2 = array![9, 2, 6, 5, 3];
        let res = PolyOps::mul_polys(poly1, poly2);
        assert!(
            res == array![27, 15, 56, 38, 85, 39, 47, 28, 15],
            "expected [27, 15, 56, 38, 85, 39, 47, 28, 15], got {:?}",
            res
        );

        let mut poly1 = array![3_u128, 1, 4, 1, 5];
        let mut poly2 = array![9, 2, 6];
        let res = PolyOps::mul_polys(poly1, poly2);
        assert!(
            res == array![27, 15, 56, 23, 71, 16, 30],
            "expected [27, 15, 56, 23, 71, 16, 30], got {:?}",
            res
        );
    }

    #[test]
    fn test_evaluation_tensor_product() {
        let pt = array![bnf(12), bnf(5)];
        let res = evaluation_tensor_product(pt.span());
        assert!(
            res == array![bnf(3), bnf(7), bnf(14), bnf(11)],
            "expected [12, 8, 15, 10], got {:?}",
            res
        );
    }
}
