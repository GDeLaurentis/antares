import pytest

from lips import Particles
from syngular import Field

from antares.terms.term import Term


def test_term_analytic_division_monomial():
    oPs = Particles(6, field=Field("finite field", 2 ** 31 - 1, 1))

    a = Term("""+(1/2⟨1|2⟩⁴[1|2][2|3]⟨3|1+2|5]⁴)/(⟨1|3⟩⁴[4|5][5|6]⟨1|2+3|4]⟨3|1+2|6]s_123)""")
    b = Term("""+(1/2[5|6]⁴⟨5|6⟩⟨4|5⟩⟨2|(1+3)|4]⁴)/(⟨1|2⟩⟨2|3⟩[4|6]⁴⟨1|(2+3)|4]⟨3|(1+2)|6]s_456)""")

    assert (a / b)(oPs) == a(oPs) / b(oPs)


def test_term_analytic_multiplication_monomial():
    oPs = Particles(6, field=Field("finite field", 2 ** 31 - 1, 1))

    a = Term("""+(1/2⟨1|2⟩⁴[1|2][2|3]⟨3|1+2|5]⁴)/(⟨1|3⟩⁴[4|5][5|6]⟨1|2+3|4]⟨3|1+2|6]s_123)""")
    b = Term("""+(1/2[5|6]⁴⟨5|6⟩⟨4|5⟩⟨2|(1+3)|4]⁴)/(⟨1|2⟩⟨2|3⟩[4|6]⁴⟨1|(2+3)|4]⟨3|(1+2)|6]s_456)""")

    assert (a * b)(oPs) == a(oPs) * b(oPs)


def test_term_analytic_division_polynomial():
    oPs = Particles(6, field=Field("finite field", 2 ** 31 - 1, 1))

    a = Term("""+(⟨1|2⟩³[2|3]⟨3|1+2|5]⁴(1/2⟨1|2⟩[1|2]+7/13⟨1|3⟩[1|3]))/(⟨1|3⟩⁴[4|5][5|6]⟨1|2+3|4]⟨3|1+2|6]s_123)""")
    b = Term("""+(1/2[5|6]⁴⟨5|6⟩⟨4|5⟩⟨2|(1+3)|4]⁴)/(⟨1|2⟩⟨2|3⟩[4|6]⁴⟨1|(2+3)|4]⟨3|(1+2)|6]s_456)""")

    assert (a / b)(oPs) == a(oPs) / b(oPs)


def test_term_analytic_multiplication_polynomial():
    oPs = Particles(6, field=Field("finite field", 2 ** 31 - 1, 1))

    a = Term("""+(⟨1|2⟩³[2|3]⟨3|1+2|5]⁴(1/2⟨1|2⟩[1|2]+7/13⟨1|3⟩[1|3]))/(⟨1|3⟩⁴[4|5][5|6]⟨1|2+3|4]⟨3|1+2|6]s_123)""")
    b = Term("""+(1/2[5|6]⁴⟨5|6⟩⟨4|5⟩⟨2|(1+3)|4]⁴)/(⟨1|2⟩⟨2|3⟩[4|6]⁴⟨1|(2+3)|4]⟨3|(1+2)|6]s_456)""")

    assert (a * b)(oPs) == a(oPs) * b(oPs)


def test_terms_string_with_typo_raises_assertion():
    bad_term = """+(1⟨6|2+4|3]⟨6|2+3|4][2|4])/(⟨5|6⟩⟨1|2+3|4]s_234[3|4]²"""
    with pytest.raises(AssertionError):
        Term(bad_term)
