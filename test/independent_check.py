"""
Independent python implementation of popgen estimators

All equation numbers reference Version 2023-04-12 of
the "Pool-Sequencing corrections for population genetics statistics"
document.

This script is intended as an independent reimplementation of those statistics
by Jeffrey P. Spence (jspence@stanford.edu). As such the code is concise and
tries to follow the notation in the equations document. It is not intended to
be beautiful or easy to read or efficient.
"""

from typing import List
import tqdm
from pandas import DataFrame
from numpy import sqrt
from scipy.stats import binom

# Speed up the repetitive part
from functools import lru_cache

@lru_cache(maxsize=None)
def harmonic(n: float, power: float) -> float:
    """Returns the nth harmonic number of type power"""
    int_n = round(n)
    return sum([1/x**power for x in range(1, int_n+1)])


@lru_cache(maxsize=None)
def binom_pmf(r, n, p):
    return binom.pmf(r, n, p)


@lru_cache(maxsize=None)
def pi_within(
    n1: int, c1: int, k1: int, n2: int, c2: int, k2: int
) -> float:
    """
    Equation (32)

    n1: number of individuals in pool 1
    c1: number of reads in pool 1 (coverage)
    k1: number of reads with "1" allele in pool 1
    n2: number of individuals in pool 2
    c2: number of reads in pool 2 (coverage)
    k2: number of reads with "1" allele in pool 2
    """
    f11 = k1/c1
    f10 = 1 - k1/c1
    pi_one = n1 / (n1 - 1) * c1 / (c1 - 1) * (1 - f11**2 - f10**2)

    f21 = k2/c2
    f20 = 1 - k2/c2
    pi_two = n2 / (n2 - 1) * c2 / (c2 - 1) * (1 - f21**2 - f20**2)

    return (pi_one + pi_two) / 2.


@lru_cache(maxsize=None)
def pi_between(
    n1: int, c1: int, k1: int, n2: int, c2: int, k2: int
) -> float:
    """
    Equation (33)

    n1: number of individuals in pool 1
    c1: number of reads in pool 1 (coverage)
    k1: number of reads with "1" allele in pool 1
    n2: number of individuals in pool 2
    c2: number of reads in pool 2 (coverage)
    k2: number of reads with "1" allele in pool 2
    """
    f11 = k1/c1
    f10 = 1 - k1/c1
    f21 = k2/c2
    f20 = 1 - k2/c2
    return 1 - f11*f21 - f10*f20


@lru_cache(maxsize=None)
def pi_total(
    n1: int, c1: int, k1: int, n2: int, c2: int, k2: int
) -> float:
    """
    Equation (34)

    n1: number of individuals in pool 1
    c1: number of reads in pool 1 (coverage)
    k1: number of reads with "1" allele in pool 1
    n2: number of individuals in pool 2
    c2: number of reads in pool 2 (coverage)
    k2: number of reads with "1" allele in pool 2
    """
    pi_w = pi_within(n1, c1, k1, n2, c2, k2)
    pi_b = pi_between(n1, c1, k1, n2, c2, k2)
    return (pi_w + pi_b) / 2.


@lru_cache(maxsize=None)
def theta_pi_denom( n, c_ell, b ):
    den = 0.
    for m in range(b, c_ell - b + 1):
        pre_factor = 2 * m * (c_ell - m) / c_ell / (c_ell - 1)
        for k_val in range(1, n):
            den += pre_factor * binom_pmf(m, c_ell, k_val/n) / k_val
    return den

def theta_pi(
    n: int,
    c: List[int],
    k: List[int],
    b: int
) -> float:
    """
    Equation (13)

    n: number of individuals in the pool
    c: list of coverages at each position
    k: number of derived alleles at each position
    b: minimum number of minor alleles to consider
    """
    assert len(c) == len(k)
    to_return = 0.
    for c_ell, k_ell in zip(c, k):
        if k_ell < b or (c_ell - k_ell) < b:
            continue
        f1 = k_ell / c_ell
        f0 = 1 - f1
        num = c_ell / (c_ell - 1) * (1 - f1**2 - f0**2)
        den = theta_pi_denom( n, c_ell, b )
        to_return += num / den
    return to_return / len(c)


@lru_cache(maxsize=None)
def theta_w_denom( n, c_ell, b ):
    den = 0.
    for m in range(b, c_ell - b + 1):
        for k_val in range(1, n):
            den += binom_pmf(m, c_ell, k_val/n) / k_val
    return den

def theta_w(
    n: int,
    c: List[int],
    k: List[int],
    b: int
) -> float:
    """
    Equation (16)

    n: number of individuals in the pool
    c: list of coverages at each position
    k: number of derived alleles at each position
    b: minimum number of minor alleles to consider
    """
    assert len(c) == len(k)
    to_return = 0.
    for c_ell, k_ell in zip(c, k):
        if k_ell < b or (c_ell - k_ell) < b:
            continue
        to_return += 1.0 / theta_w_denom( n, c_ell, b )

        # den = 0.
        # for m in range(b, c_ell - b + 1):
        #     for k_val in range(1, n):
        #         den += binom_pmf(m, c_ell, k_val/n) / k_val
        # to_return += 1/den

    # Up until the equations document version 2023-08-02, we had an algebraic
    # mistake in the equation here, where we had an additional harmonic a_n
    # in the final estimator. It was just a mistake in the document, but made
    # it into the code here (line below). We fixed this now, but keep this line
    # here in order to back-check later in case that's needed.
    # See also the `test-watterson` directory here for the simulation test
    # that we did to make sure that this was just an algebraic mistake.
    a1 = 1.0 # harmonic(n-1, 1)
    return to_return / a1 / len(c)


@lru_cache(maxsize=None)
def f_star(n: int) -> float:
    """Equation (20)"""
    a1 = harmonic(n-1, 1)
    return (n-3) / (a1 * (n-1) - n)


def n_tilde(n: int, c: int) -> float:
    """Equation (25)"""
    return n * (1 - ((n-1)/n)**c)


@lru_cache(maxsize=None)
def alpha_star(n: int) -> float:
    """Equation (21)"""
    f = f_star(n)
    a1 = harmonic(n-1, 1)
    term_1 = f**2 * (a1 - n / (n-1))
    term_2 = f * (a1*4*(n+1)/(n-1)**2 - 2*(n+3)/(n-1))
    term_3 = -a1*8*(n+1)/n/(n-1)**2
    term_4 = (n**2 + n + 60) / (3 * n * (n-1))
    return term_1 + term_2 + term_3 + term_4


@lru_cache(maxsize=None)
def beta_star(n: int) -> float:
    """Equation (22)"""
    f = f_star(n)
    a2 = harmonic(n-1, 2)
    a1 = harmonic(n-1, 1)
    term_1 = f**2 * (a2 - (2*n-1)/(n-1)**2)
    term_2 = f * (a1*8/(n-1)
                  - a1*4/n/(n-1)
                  - (n**3 + 12*n**2 - 35*n + 18) / n / (n-1)**2)
    term_3 = -a1*16/n/(n-1)
    term_4 = a1*8/n**2/(n-1)
    term_5 = 2 * (n**4 + 110*n**2 - 255*n + 126) / (9 * n**2 * (n-1)**2)
    return term_1 + term_2 + term_3 + term_4 + term_5


def achaz_var_d(
    n: int,
    c: List[int],
    k: List[int],
    b: int
) -> float:
    """
    Equation (26)

    n: number of individuals in the pool
    c: list of coverages at each position
    k: number of derived alleles at each position
    b: minimum number of minor alleles to consider
    """
    assert b == 2
    assert len(c) == len(k)
    theta_hat_w = theta_w(n, c, k, b)
    effective_n = n_tilde(n, b)
    a_star = alpha_star(effective_n)
    b_star = beta_star(effective_n)
    return a_star / len(c) * theta_hat_w + b_star * theta_hat_w**2


def tajimas_d(
    n: int,
    c: List[int],
    k: List[int],
    b: int
) -> float:
    """
    Equation (26)

    n: number of individuals in the pool
    c: list of coverages at each position
    k: number of derived alleles at each position
    b: minimum number of minor alleles to consider
    """
    d_pool = theta_pi(n, c, k, b) - theta_w(n, c, k, b)
    var_pool = achaz_var_d(n, c, k, b)
    return d_pool / sqrt(var_pool)


def fst_nei(
    n1: int,
    c1: List[int],
    k1: List[int],
    n2: int,
    c2: List[int],
    k2: List[int]
) -> float:
    """
    Equation (37)

    n1: number of individuals in pool 1
    c1: list of coverages at each position in the window in pool 1
    k1: list of number of derived alleles at each position in pool 1
    n2: number of individuals in pool 2
    c2: list of coverages at each position in the window in pool 2
    k2: list of number of derived alleles at each position in pool 2
    """
    assert len(c1) == len(k1)
    assert len(c1) == len(c2)
    assert len(c1) == len(k2)
    pi_ws = [pi_within(n1, c1ell, k1ell, n2, c2ell, k2ell)
             for c1ell, k1ell, c2ell, k2ell in zip(c1, k1, c2, k2)]
    pi_ts = [pi_total(n1, c1ell, k1ell, n2, c2ell, k2ell)
             for c1ell, k1ell, c2ell, k2ell in zip(c1, k1, c2, k2)]
    return 1 - sum(pi_ws) / sum(pi_ts)


def fst_hudson(
    n1: int,
    c1: List[int],
    k1: List[int],
    n2: int,
    c2: List[int],
    k2: List[int]
) -> float:
    """
    Equation (38)

    n1: number of individuals in pool 1
    c1: list of coverages at each position in the window in pool 1
    k1: list of number of derived alleles at each position in pool 1
    n2: number of individuals in pool 2
    c2: list of coverages at each position in the window in pool 2
    k2: list of number of derived alleles at each position in pool 2
    """
    assert len(c1) == len(k1)
    assert len(c1) == len(c2)
    assert len(c1) == len(k2)
    pi_ws = [pi_within(n1, c1ell, k1ell, n2, c2ell, k2ell)
             for c1ell, k1ell, c2ell, k2ell in zip(c1, k1, c2, k2)]
    pi_bs = [pi_between(n1, c1ell, k1ell, n2, c2ell, k2ell)
             for c1ell, k1ell, c2ell, k2ell in zip(c1, k1, c2, k2)]
    return 1 - sum(pi_ws) / sum(pi_bs)


def karlsson_h(n: int, c: int, k: int) -> float:
    """
    Unlabeled equation between (50) and (51)

    n: number of individuals in pool (intentionally unused)
    c: coverage
    k: number of derived alleles
    """
    u = k/c
    v = 1 - k/c
    return c / (c-1) * u * v


def karlsson_n(
    n1: int, c1: int, k1: int, n2: int, c2: int, k2: int
) -> float:
    """
    Equation (49)

    n1: number of individuals in pool 1
    c1: list of coverages at each position in the window in pool 1
    k1: list of number of derived alleles at each position in pool 1
    n2: number of individuals in pool 2
    c2: list of coverages at each position in the window in pool 2
    k2: list of number of derived alleles at each position in pool 2
    """
    u1 = k1/c1
    u2 = k2/c2
    h1 = karlsson_h(n1, c1, k1)
    h2 = karlsson_h(n2, c2, k2)
    return (u1 - u2)**2 - (h1/c1 + h2/c2)


def karlsson_d(
    n1: int, c1: int, k1: int, n2: int, c2: int, k2: int
) -> float:
    """
    Equation (50)

    n1: number of individuals in pool 1
    c1: list of coverages at each position in the window in pool 1
    k1: list of number of derived alleles at each position in pool 1
    n2: number of individuals in pool 2
    c2: list of coverages at each position in the window in pool 2
    k2: list of number of derived alleles at each position in pool 2
    """
    h1 = karlsson_h(n1, c1, k1)
    h2 = karlsson_h(n2, c2, k2)
    return karlsson_n(n1, c1, k1, n2, c2, k2) + h1 + h2


def fst_karlsson(
    n1: int,
    c1: List[int],
    k1: List[int],
    n2: int,
    c2: List[int],
    k2: List[int]
) -> float:
    """
    Equation (51)

    n1: number of individuals in pool 1
    c1: list of coverages at each position in the window in pool 1
    k1: list of number of derived alleles at each position in pool 1
    n2: number of individuals in pool 2
    c2: list of coverages at each position in the window in pool 2
    k2: list of number of derived alleles at each position in pool 2
    """
    assert len(c1) == len(k1)
    assert len(c1) == len(c2)
    assert len(c1) == len(k2)
    ns = [karlsson_n(n1, c1ell, k1ell, n2, c2ell, k2ell)
          for c1ell, k1ell, c2ell, k2ell in zip(c1, k1, c2, k2)]
    ds = [karlsson_d(n1, c1ell, k1ell, n2, c2ell, k2ell)
          for c1ell, k1ell, c2ell, k2ell in zip(c1, k1, c2, k2)]
    return sum(ns) / sum(ds)


def pi_kofler(n: int, c: int, k: float) -> float:
    """
    Equations (39), (40), and (41)

    n: number of individuals in pool (intentionally unused)
    c: coverage
    k: number of derived alleles
    """
    return c / (c - 1) * (1 - (k/c)**2 - (1 - k/c)**2)


def fst_kofler(
    n1: int,
    c1: List[int],
    k1: List[int],
    n2: int,
    c2: List[int],
    k2: List[int]
) -> float:
    """
    Equations (42), (43), (44), and (45)

    n1: number of individuals in pool 1
    c1: list of coverages at each position in the window in pool 1
    k1: list of number of derived alleles at each position in pool 1
    n2: number of individuals in pool 2
    c2: list of coverages at each position in the window in pool 2
    k2: list of number of derived alleles at each position in pool 2
    """
    pi_w1 = n1 / (n1 - 1) * sum(
        [pi_kofler(n1, c1ell, k1ell) for c1ell, k1ell in zip(c1, k1)]
    )
    pi_w2 = n2 / (n2 - 1) * sum(
        [pi_kofler(n1, c2ell, k2ell) for c2ell, k2ell in zip(c2, k2)]
    )
    min_n = min([n1, n2])
    min_c = [min([c1ell, c2ell]) for c1ell, c2ell in zip(c1, c2)]
    fs = [(k1ell/c1ell + k2ell/c2ell) / 2
          for k1ell, c1ell, k2ell, c2ell in zip(k1, c1, k2, c2)]
    pi_t = min_n / (min_n - 1) * sum(
        [pi_kofler(min_n, c_ell, c_ell*f) for c_ell, f in zip(min_c, fs)]
    )

    return (pi_t - (pi_w1 + pi_w2)/2) / pi_t


def main():
    n1_list = []
    n2_list = []
    c1_list = []
    c2_list = []
    k1_list = []
    k2_list = []
    w_list = []
    theta_pi_one_list = []
    theta_pi_two_list = []
    theta_w_one_list = []
    theta_w_two_list = []
    achaz_var_d_one_list = []
    achaz_var_d_two_list = []
    tajimas_d_one_list = []
    tajimas_d_two_list = []
    fst_nei_list = []
    fst_hudson_list = []
    fst_karlsson_list = []
    fst_kofler_list = []

    run_vars = []
    for n1 in [10, 100]:
        for n2 in [5, 50]:
            for c1 in [10, 100]:
                for c2 in [5, 50]:
                    for f1 in [0.1, 0.25, 0.5, 0.75, 0.9]:
                        for f2 in [0.2, 0.4, 0.6, 0.8]:
                            for w in [1, 10, 100]:
                                k1 = round(f1*c1)
                                k2 = round(f2*c2)
                                run_vars.append((n1, n2, c1, c2, k1, k2, w))

    for run_var in tqdm.tqdm(run_vars):
        n1, n2, c1, c2, k1, k2, w = run_var
        n1_list.append(n1)
        n2_list.append(n2)
        c1_list.append(c1)
        c2_list.append(c2)
        k1_list.append(k1)
        k2_list.append(k2)
        w_list.append(w)
        c1s = [c1] * w
        c2s = [c2] * w
        k1s = [k1] + [0] * (w-1)
        k2s = [k2] + [0] * (w-1)
        theta_pi_one_list.append(
            theta_pi(n1, c1s, k1s, 2)
        )
        theta_pi_two_list.append(
            theta_pi(n2, c2s, k2s, 2)
        )
        theta_w_one_list.append(
            theta_w(n1, c1s, k1s, 2)
        )
        theta_w_two_list.append(
            theta_w(n2, c2s, k2s, 2)
        )
        achaz_var_d_one_list.append(
            achaz_var_d(n1, c1s, k1s, 2)
        )
        achaz_var_d_two_list.append(
            achaz_var_d(n2, c2s, k2s, 2)
        )
        tajimas_d_one_list.append(
            tajimas_d(n1, c1s, k1s, 2)
        )
        tajimas_d_two_list.append(
            tajimas_d(n2, c2s, k2s, 2)
        )
        fst_nei_list.append(
            fst_nei(n1, c1s, k1s, n2, c2s, k2s)
        )
        fst_hudson_list.append(
            fst_hudson(n1, c1s, k1s, n2, c2s, k2s)
        )
        fst_karlsson_list.append(
            fst_karlsson(n1, c1s, k1s, n2, c2s, k2s)
        )
        fst_kofler_list.append(
            fst_kofler(n1, c1s, k1s, n2, c2s, k2s)
        )
    to_save = DataFrame()
    to_save['pool_size_1'] = n1_list
    to_save['pool_size_2'] = n2_list
    to_save['coverage_1'] = c1_list
    to_save['coverage_2'] = c2_list
    to_save['derived_count_1'] = k1_list
    to_save['derived_count_2'] = k2_list
    to_save['window_size'] = w_list
    to_save['theta_pi_1'] = theta_pi_one_list
    to_save['theta_pi_2'] = theta_pi_two_list
    to_save['theta_w_1'] = theta_w_one_list
    to_save['theta_w_2'] = theta_w_two_list
    to_save['achaz_var_1'] = achaz_var_d_one_list
    to_save['achaz_var_2'] = achaz_var_d_two_list
    to_save['tajimas_d_1'] = tajimas_d_one_list
    to_save['tajimas_d_2'] = tajimas_d_two_list
    to_save['fst_nei'] = fst_nei_list
    to_save['fst_hudson'] = fst_hudson_list
    to_save['fst_karlsson'] = fst_karlsson_list
    to_save['fst_kofler'] = fst_kofler_list
    to_save.to_csv('independent_check_statistics.tsv', index=False, sep='\t')


if __name__ == '__main__':
    main()
