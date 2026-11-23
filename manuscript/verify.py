import json
from fractions import Fraction
from itertools import combinations

import sympy as sp


def zero_matrix(matrix):
    return all(sp.simplify(entry) == 0 for entry in matrix)


def basis(dimension):
    matrices = []
    for first, second in combinations(range(dimension), 2):
        symmetric = sp.zeros(dimension)
        symmetric[first, second] = symmetric[second, first] = 1 / sp.sqrt(2)
        antisymmetric = sp.zeros(dimension)
        antisymmetric[first, second] = -sp.I / sp.sqrt(2)
        antisymmetric[second, first] = sp.I / sp.sqrt(2)
        matrices.extend([symmetric, antisymmetric])
    for level in range(1, dimension):
        diagonal = [1] * level + [-level] + [0] * (dimension - level - 1)
        matrices.append(sp.diag(*diagonal) / sp.sqrt(level * (level + 1)))
    return matrices


def block_matrix(blocks):
    return sp.Matrix.vstack(*(sp.Matrix.hstack(*row) for row in blocks))


def matrix_checks(dimension):
    matrices = basis(dimension)
    count = len(matrices)
    identity = sp.eye(dimension)
    for first, left in enumerate(matrices):
        assert left.H == left and sp.trace(left) == 0
        for second, right in enumerate(matrices):
            assert sp.simplify(sp.trace(left * right)) == int(first == second)
    for first in range(dimension):
        for second in range(dimension):
            unit = sp.zeros(dimension)
            unit[first, second] = 1
            contraction = sum((matrix * unit * matrix for matrix in matrices), sp.zeros(dimension))
            assert zero_matrix(contraction - sp.trace(unit) * identity + unit / dimension)
    casimir = sum((matrix * matrix for matrix in matrices), sp.zeros(dimension))
    assert zero_matrix(casimir - sp.Rational(count, dimension) * identity)
    scale = sp.Rational(dimension + 1, dimension + 2)
    ratio = sp.Rational(1, dimension + 1)
    involution = block_matrix([
        [right * left - ratio * left * right for right in matrices]
        for left in matrices
    ])
    assert zero_matrix(involution - involution.H)
    assert zero_matrix(involution * involution - sp.eye(count * dimension))
    primal = [
        [scale * (left * right + right * left + int(first == second) * identity)
         for second, right in enumerate(matrices)]
        for first, left in enumerate(matrices)
    ]
    constraint = block_matrix([
        [primal[first][second] - left * right for second, right in enumerate(matrices)]
        for first, left in enumerate(matrices)
    ])
    assert zero_matrix(constraint - scale * (sp.eye(count * dimension) + involution))
    assert zero_matrix(sum((primal[index][index] for index in range(count)), sp.zeros(dimension))
                       - sp.Rational(count * (dimension + 1), dimension) * identity)
    symbols = sp.symbols('lambda_0:' + str(dimension), real=True)
    state = sp.diag(*symbols)
    covariance = sp.Matrix(count, count, lambda first, second: sp.trace(state * matrices[first] * matrices[second]))
    imaginary = covariance.applyfunc(sp.im)
    expected = sp.zeros(count)
    for pair_index, (first, second) in enumerate(combinations(range(dimension), 2)):
        expected[2 * pair_index, 2 * pair_index + 1] = (symbols[first] - symbols[second]) / 2
        expected[2 * pair_index + 1, 2 * pair_index] = -(symbols[first] - symbols[second]) / 2
    assert zero_matrix(imaginary - expected)
    return {'dimension': dimension, 'basis_size': count, 'block_size': count * dimension,
            'identities': ['orthonormality', 'completeness_on_matrix_units', 'casimir',
                           'hermitian_involution', 'primal_constraint_factorization',
                           'primal_diagonal_sum', 'symbolic_imaginary_covariance']}


def compositions(total, count):
    if count == 1:
        yield (total,)
    else:
        for first in range(total + 1):
            for remaining in compositions(total - first, count - 1):
                yield (first,) + remaining


def spectral_checks(dimension, denominator):
    cases = 0
    for numerator in compositions(denominator, dimension):
        spectrum = [Fraction(value, denominator) for value in numerator]
        purity = sum(value * value for value in spectrum)
        absolute_gaps = sum(abs(spectrum[first] - spectrum[second])
                            for first, second in combinations(range(dimension), 2))
        squared_gaps = sum((spectrum[first] - spectrum[second]) ** 2
                           for first, second in combinations(range(dimension), 2))
        assert squared_gaps == dimension * purity - 1
        assert absolute_gaps >= squared_gaps
        assert absolute_gaps ** 2 >= (dimension - 1) * squared_gaps
        holevo = dimension - purity + absolute_gaps
        lower = (dimension - 1) * (1 + purity)
        upper = dimension ** 2 + dimension - 1 - purity
        assert holevo >= lower > 0
        assert upper / holevo <= upper / lower <= dimension + 1
        assert upper / lower < dimension + 1 or purity == Fraction(1, dimension)
        assert (dimension + 1) * lower - upper == dimension * (dimension * purity - 1)
        assert Fraction(dimension ** 2 - 1, dimension) <= holevo <= 2 * (dimension - 1)
        assert (holevo == 2 * (dimension - 1)) == (purity == 1)
        ordered = sorted(spectrum, reverse=True)
        nonzero_gaps = [index for index in range(1, dimension)
                        if ordered[index - 1] != ordered[index]]
        equality = absolute_gaps ** 2 == (dimension - 1) * squared_gaps
        assert equality == (not nonzero_gaps or
                            len(nonzero_gaps) == 1 and nonzero_gaps[0] in (1, dimension - 1))
        stability = Fraction(dimension + 1) - upper / lower
        assert stability == dimension * (dimension * purity - 1) / lower
        assert stability >= dimension ** 2 * (purity - Fraction(1, dimension)) / (2 * (dimension - 1))
        cases += 1
    return {'dimension': dimension, 'denominator': denominator, 'spectra': cases}


def attaining_family_checks():
    cases = 0
    for dimension in range(2, 9):
        for numerator in range(13):
            depolarization = Fraction(numerator, 12)
            for gap, multiplicity in [(depolarization, 1),
                                      (depolarization / (dimension - 1), dimension - 1)]:
                higher = Fraction(1, dimension) + Fraction(dimension - multiplicity, dimension) * gap
                lower = Fraction(1, dimension) - Fraction(multiplicity, dimension) * gap
                spectrum = [higher] * multiplicity + [lower] * (dimension - multiplicity)
                assert min(spectrum) >= 0 and sum(spectrum) == 1
                purity = sum(value ** 2 for value in spectrum)
                gaps = sum(abs(spectrum[first] - spectrum[second])
                           for first, second in combinations(range(dimension), 2))
                assert gaps ** 2 == (dimension - 1) * (dimension * purity - 1)
                if multiplicity == dimension - 1:
                    assert purity <= Fraction(1, dimension - 1)
                cases += 1
    return {'dimensions': [2, 3, 4, 5, 6, 7, 8], 'parameter_denominator': 12,
            'spectra': cases, 'families': ['one larger eigenvalue', 'one smaller eigenvalue']}


def purity_ordering_counterexample():
    first = [Fraction(3, 4), Fraction(1, 5), Fraction(1, 20)]
    second = [Fraction(4, 5), Fraction(1, 10), Fraction(1, 10)]
    values = []
    for spectrum in [first, second]:
        purity = sum(value ** 2 for value in spectrum)
        gaps = sum(abs(spectrum[left] - spectrum[right]) for left, right in combinations(range(3), 2))
        values.append((purity, 3 - purity + gaps))
    assert values == [(Fraction(121, 200), Fraction(759, 200)),
                      (Fraction(33, 50), Fraction(187, 50))]
    assert values[0][0] < values[1][0] and values[0][1] > values[1][1]
    return {'dimension': 3, 'full_rank': True,
            'spectra': [[str(value) for value in first], [str(value) for value in second]],
            'purities': [str(pair[0]) for pair in values],
            'holevo_costs': [str(pair[1]) for pair in values]}


def design_checks(dimension):
    projectors = []
    for index in range(dimension):
        vector = sp.eye(dimension)[:, index]
        projectors.append(vector * vector.H)
    if dimension == 2:
        for phase in [sp.Integer(1), -sp.Integer(1), sp.I, -sp.I]:
            vector = sp.Matrix([1, phase]) / sp.sqrt(2)
            projectors.append(vector * vector.H)
    else:
        omega = -sp.Rational(1, 2) + sp.I * sp.sqrt(3) / 2
        for coefficient in range(dimension):
            for linear in range(dimension):
                vector = sp.Matrix([sp.expand(omega ** ((coefficient * index ** 2 + linear * index) % dimension))
                                    for index in range(dimension)]) / sp.sqrt(dimension)
                projectors.append((vector * vector.H).applyfunc(sp.simplify))
    weight = sp.Rational(1, len(projectors))
    identity = sp.eye(dimension)
    swap = sp.zeros(dimension ** 2)
    for first in range(dimension):
        for second in range(dimension):
            swap[first * dimension + second, second * dimension + first] = 1
    moment = sum((weight * sp.kronecker_product(projector, projector) for projector in projectors),
                 sp.zeros(dimension ** 2))
    assert zero_matrix(moment - (sp.eye(dimension ** 2) + swap) / (dimension * (dimension + 1)))
    effects = [dimension * weight * projector for projector in projectors]
    assert zero_matrix(sum(effects, sp.zeros(dimension)) - identity)
    matrices = basis(dimension)
    scores = [[sp.simplify((dimension + 1) * sp.trace(matrix * projector)) for matrix in matrices]
              for projector in projectors]
    for index, matrix in enumerate(matrices):
        reconstruction = sum((score[index] * effect for score, effect in zip(scores, effects)),
                             sp.zeros(dimension))
        assert zero_matrix(reconstruction - matrix)
    for score in scores:
        assert sp.simplify(sum(value ** 2 for value in score)) == sp.Rational((dimension ** 2 - 1) * (dimension + 1), dimension)
    fisher_values = []
    spectra = [[sp.Rational(1, dimension)] * dimension,
               [sp.Rational(dimension - index, dimension * (dimension + 1) // 2) for index in range(dimension)]]
    for spectrum in spectra:
        state = sp.diag(*spectrum)
        probabilities = [sp.simplify(sp.trace(state * effect)) for effect in effects]
        derivatives = [[sp.simplify(sp.trace(matrix * effect)) for matrix in matrices] for effect in effects]
        fisher = sp.zeros(len(matrices))
        for probability, derivative in zip(probabilities, derivatives):
            vector = sp.Matrix(derivative)
            fisher += vector * vector.T / probability
        fisher = fisher.applyfunc(sp.simplify)
        cost = sp.simplify(sp.trace(fisher.inv()))
        purity = sum(value ** 2 for value in spectrum)
        upper = dimension ** 2 + dimension - 1 - purity
        assert sp.simplify(upper - cost) >= 0
        if spectrum == spectra[0]:
            assert zero_matrix(fisher - sp.Rational(dimension, dimension + 1) * sp.eye(len(matrices)))
            assert cost == upper
        fisher_values.append({'purity': str(purity), 'fisher_inverse_trace': str(cost),
                              'estimator_total_variance': str(upper)})
    return {'dimension': dimension, 'outcomes': len(projectors),
            'identities': ['weighted_two_design_moment', 'POVM_normalization',
                           'global_unbiasedness', 'constant_score_norm'],
            'full_rank_fisher_examples': fisher_values}


def main():
    result = {
        'coverage': 'Exact finite matrix identities and rational-spectrum inequalities; all-dimensional conclusions have written proofs in note.pdf.',
        'arithmetic': 'SymPy exact algebra and Python Fraction',
        'matrix_checks': [matrix_checks(dimension) for dimension in (2, 3)],
        'spectral_checks': [spectral_checks(dimension, 12) for dimension in range(2, 6)],
        'spectral_inequalities': ['sum_absolute_gaps_squared >= (d-1)(d*P-1)',
                                 'equality classification by ordered spectral gaps',
                                 'global Holevo range and pure-state equality',
                                 'quantitative gap below d+1'],
        'attaining_families': attaining_family_checks(),
        'purity_ordering_counterexample': purity_ordering_counterexample(),
        'design_checks': [design_checks(dimension) for dimension in (2, 3)],
        'status': 'pass',
    }
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
