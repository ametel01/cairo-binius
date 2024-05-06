# Binary Field Arithmetic and Polynomial Operations

This repository contains Cairo code for performing arithmetic operations on binary field elements and polynomial operations. The code is divided into three main sections:

1. Binary Field Arithmetic
2. Additive Number Theoretic Transform (NTT)
3. Polynomial Operations and Utilities

## Sources

The code in this repository is inspired by the following research paper:

- "Succinct Arguments over Towers of Binary Fields" by Benjamin E. Diamond and Jim Posen . Available at: [https://eprint.iacr.org/2023/1784.pdf](https://eprint.iacr.org/2023/1784.pdf)

Additionally, a Python implementation of the concepts discussed in the paper can be found in the following GitHub repository:

- [ethereum/research/binius](https://github.com/ethereum/research/tree/master/binius) by Vitalik Buterin

## Binary Field Arithmetic

The `BinaryFieldElement` struct represents an element in a binary field. It provides implementations for various arithmetic operations, including addition, subtraction, multiplication, division, exponentiation, and inversion. The code also includes a custom `binmul` function for efficient binary field multiplication.

## Additive Number Theoretic Transform (NTT

The code includes functions for performing additive NTT and its inverse. The `get_Wi` function generates the necessary coefficients for the NTT, while `get_Wi_eval` evaluates the coefficients at a specific point. The `additive_ntt` and `inv_additive_ntt` functions perform the forward and inverse NTT, respectively. The `extend` function extends a sequence of binary field elements by a specified expansion factor using the NTT.

## Polynomial Operations and Utilities

The code provides utility functions for polynomial operations, including:

- `multilinear_poly_eval`: Evaluates a multilinear polynomial at a given point.
- `log2`: Computes the logarithm base 2 of a power of 2.
- `PolyOps` trait: Defines operations for adding, multiplying, and evaluating polynomials.
- `evaluation_tensor_product`: Computes the evaluation tensor product of a sequence of binary field elements.
- `max`: Returns the maximum of two `u128` values.


## Testing

The code includes test modules for each section, which contain unit tests to verify the correctness of the implemented functions. To run the tests, use the `snforge test` command.