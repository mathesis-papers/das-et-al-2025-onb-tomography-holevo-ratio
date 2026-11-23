# Fixed-purity Holevo minima and the full-ONB conjecture of Das et al. (2025)

`note.tex` proves the full-ONB conjecture immediately after Theorem 3 in
Section 3.2.1 of Das et al., Quantum 9, 1867 (2025), arXiv:2405.09622v2.
It extends their analysis by proving the Holevo minimality that they already
identified as sufficient in Section 5.4, immediately after Figure 5.
The spectral evaluation specializes established D-invariant Holevo/RLD theory
(Suzuki, Entropy 21, 703 (2019), Section 4.1.2, Proposition 2).
The exact fixed-purity minimum and the individual-measurement extension are
the quantitative conclusions; the reduction, feasible NH blocks, and
weighted-two-design framework are credited to their earlier sources.
Scott (2006), Section VI, Equations (101)--(103), already gives the exact
estimator variance used here: Hilbert--Schmidt squared reconstruction error
is the sum of squared ONB-coordinate errors. The individual-measurement
ratio bound combines that prior estimator formula with the Holevo minimum.
For dimension `d >= 2`, identity weight, and purity `P`, it establishes

```text
C_H = d - P + sum_{a<b} |lambda_a - lambda_b|
C_H >= (d - 1)(1 + P)
C_NH <= d^2 + d - 1 - P
m_d(P) = d - P + sqrt((d - 1)(d*P - 1))
C_H >= m_d(P) >= (d - 1)(1 + P)
C_NH / C_H <= (d^2 + d - 1 - P) / m_d(P) <= d + 1.
```

The last inequality is strict away from the maximally mixed state.
The same upper bounds hold for `C_MI / C_H`, where `C_MI` is the optimal
individual-measurement cost for full-rank states. A finite weighted
projective two-design supplies a globally unbiased estimator of total
variance `d^2 + d - 1 - P`, without a SIC-existence assumption. This also
improves Theorem 5 of the source from `d+2` to the sharp `d+1`.
The Gill--Massar/Holevo ratio inherits the same purity envelope and
also attains `d+1` exactly at the maximally mixed state.
The function `m_d(P)` is the exact minimum Holevo cost at purity `P`, with
all minimizing spectra described in the note. The ratio bound is strictly
decreasing with purity, without asserting the same for the actual maximal
NH/H ratio at each purity. The note also proves the exact global Holevo
range and gives two full-rank qutrit states showing that purity does not
order individual Holevo costs.
The sharp maximally mixed NH value is credited to Das et al., Lemma 3.
The model uses a full traceless Hermitian Hilbert--Schmidt ONB normalized
by `tr(B_j B_k) = delta_jk`. The main theorem concerns the NHCRB/HCRB ratio;
the MI/H ratio is a written corollary with an explicit POVM argument. Arbitrary weights,
arbitrary lower-dimensional models, and fixed-purity NH extremizers are
outside the theorem's coverage; the PDF supplies exact source locators.

## Build the PDF

Use Tectonic 0.17.0:

```sh
tectonic note.tex
```

Or use a standard TeX distribution:

```sh
pdflatex -interaction=nonstopmode -halt-on-error note.tex
pdflatex -interaction=nonstopmode -halt-on-error note.tex
```

The bibliography is embedded in the TeX source; BibTeX is unnecessary.

## Exact finite checks

Use Python 3.10 or newer and SymPy 1.14.0:

```sh
python3 -m venv .venv
.venv/bin/python -m pip install sympy==1.14.0
.venv/bin/python verify.py > verification.json
```

The checker verifies exact matrix identities in dimensions 2 and 3 and
all rational probability vectors with denominator 12 in dimensions 2--5.
It checks the squared form of the fixed-purity spectral inequality and
its equality classification, 182 attaining spectra in dimensions 2--8,
the global cost range, the exact full-rank counterexample, and the
quantitative stability inequality.
It also verifies explicit six- and twelve-outcome two-design POVMs in
dimensions two and three, their global unbiasedness and score norms,
and exact Fisher-matrix costs at four full-rank states.
`verification.json` records the finite ranges. These checks complement
the written all-dimensional arguments and do not replace them.

The declaration
`D5.S3.Quantum.Measurements.OnbTomographyHolevoRatio.result` in the
[pinned Lean file](https://github.com/the-omega-institute/trureturing/blob/3d88df0f99ef3c8bfb034fa52be669b13746a844/D5/S3/Quantum/Measurements/OnbTomographyHolevoRatio.lean)
proves `C_NH <= (d + 1) C_H` for the full identity-weight model and
positive semidefinite states. The spectral identity, both purity-dependent
ratio bounds, the exact fixed-purity minimum, the MI-cost bound, the global range,
the counterexample, and the stability estimate have written proofs in the note.

Creators: Haobo Ma (ChronoAI Pte Ltd, Singapore), then Wenlin Zhang
(National University of Singapore; ChronoAI Pte Ltd, Singapore).
License: Creative Commons Attribution 4.0 International (CC BY 4.0).
