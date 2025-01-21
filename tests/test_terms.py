import pytest

from lips import Particles
from lips.fields import Field

from antares.core.settings import settings
from antares.terms.terms import Term, Terms


def test_str_and_rstr_big_powers():
    oPs = Particles(6, field=Field("finite field", 2 ** 31 - 1, 1))
    string = "(1⟨2|3⟩⟨4|5⟩^10[4|5]⁷[3|5])/(⟨1|2⟩⟨1|4⟩³⟨1|5⟩⟨2|3⟩⟨2|4⟩³[1|3]²[1|4][2|3][2|5]⟨2|1+3|2]³)"
    assert Terms(string)(oPs) == Terms(str(Terms(string)))(oPs)


def test_terms_hash():
    oTerms1 = Terms("+(+1⟨4|(2+3)|1]³)/([1|6]⟨2|3⟩⟨3|4⟩[5|6]⟨2|(1+6)|5]s_234)")
    oTerms2 = Terms("+(+1⟨6|(1+2)|3]³)/(⟨1|2⟩⟨1|6⟩[3|4][4|5]⟨2|(1+6)|5]s_345)")
    assert hash(oTerms1) == hash(oTerms1)
    assert hash(oTerms1) != hash(oTerms2)
    assert oTerms1 == oTerms1
    assert oTerms1 != oTerms2


def test_terms_phase_weights_and_mass_dimensions():
    settings.field = Field("finite field", 2 ** 31 - 1, 1)
    oTerms = Terms("""+(1)/([1|6]⟨2|3⟩⟨3|4⟩[5|6]⟨2|(1+6)|5]s_234)
    +(1)/(⟨1|2⟩⟨1|6⟩[3|4][4|5]⟨2|(1+6)|5]s_345)
    +(1)/(⟨1|2⟩⟨1|6⟩[1|6]⟨2|3⟩⟨3|4⟩[3|4][4|5][5|6]⟨2|(1+6)|5]s_234s_345)""")
    oTerms.multiplicity = 6
    with pytest.raises(ValueError):
        oTerms.mass_dimension
    with pytest.raises(ValueError):
        oTerms.phase_weights
    assert oTerms.mass_dimensions == [-8.0, -8.0, -14.0]
    assert oTerms.lphase_weights == [[1, -2, -2, -1, 2, 2], [-2, -2, 1, 2, 2, -1], [-1, -3, -1, 1, 3, 1]]
    settings.field = Field("mpc", 0, 300)


def test_terms_rstr_ansatz_template():
    assert Terms(
        str(Terms("""
        +(1)/(⟨1|2⟩⟨1|3⟩⟨2|3⟩[1|4][3|5]⟨1|2+4|1])
        +(1)/(⟨1|2⟩⟨1|3⟩⟨2|3⟩[2|5][3|4]⟨2|1+5|2])
        """))) == Terms("""
        +(1)/(⟨1|2⟩⟨1|3⟩⟨2|3⟩[1|4][3|5]⟨1|2+4|1])
        +(1)/(⟨1|2⟩⟨1|3⟩⟨2|3⟩[2|5][3|4]⟨2|1+5|2])
        """)


def test_terms_image_vs_particles_image():
    oPs = Particles(5, field=Field("finite field", 2 ** 31 - 1, 1))
    # this is a basis function from a permutation of +-+++ 2L Nc1 Nf0
    oTerms = Terms("""
        +(-4[1|4]⟨3|4⟩²⟨3|5⟩[3|5])/(⟨1|4⟩⟨2|4⟩⟨2|5⟩⟨4|5⟩²[4|5])
        +('12354', False, '-')
        +(-4[1|3]⟨3|4⟩²⟨3|5⟩²)/(⟨1|3⟩⟨2|4⟩⟨2|5⟩⟨4|5⟩³)""")
    assert oTerms(oPs.image(('13452', False))) == oTerms.Image(("15234", False, ))(oPs)  # 13452 and 15234 are each other's inverse


def test_terms_explicit_representation():
    oPs = Particles(5, field=Field("finite field", 2 ** 31 - 1, 1))
    # this is the rational part of 1L 5g ++-+- in D=4
    oTerms = Terms("""+(-1/3[1|4]²[2|4]³⟨4|5⟩)/([2|3][3|4][4|5]⟨4|(1+5)|4]²)
        +(-1/3[1|4]²[2|4]²⟨4|5⟩)/(⟨2|4⟩[2|3][3|4][4|5]⟨4|(1+5)|4])
        +(-2/3[1|2][2|4]³⟨2|5⟩²)/(⟨1|2⟩[2|3][3|4]⟨2|(1+5)|2]²)
        +(+1/3[1|4]²[1|2]⟨1|3⟩⟨1|5⟩[2|5])/(⟨1|2⟩⟨1|4⟩[1|5][2|3][4|5]⟨1|(2+3)|1])
        +(+1[2|4]²⟨4|5⟩⟨2|5⟩²)/(⟨1|2⟩⟨1|5⟩⟨2|4⟩²[2|3][3|4])
        +(-1/3⟨3|4⟩[2|4]²[1|3][1|4])/(⟨1|4⟩⟨2|4⟩[1|5][2|3][3|4][4|5])
        +(+1/3⟨1|3⟩[2|4]²[1|3][1|4])/(⟨1|2⟩⟨1|4⟩[1|5][2|3][3|4][4|5])
        +('21543', False, '-')
        +(-2/9⟨3|5⟩⁴)/(⟨1|2⟩⟨1|5⟩⟨2|3⟩⟨3|4⟩⟨4|5⟩)""")
    oTerms.multiplicity = 5
    assert oTerms(oPs) == oTerms.explicit_representation()(oPs)
    assert oTerms(oPs) == (oTerms[:7] - oTerms[:7].Image(oTerms[7].tSym) + oTerms[8:])(oPs)


def test_terms_with_trace():
    coeff_v1 = Terms("""
    +(-8⟨1|3⟩³Δ_123|45|67)/(⟨1|2⟩⟨2|3⟩⟨1|4+5|6+7|1⟩²)
    +(-4⟨1|3⟩⟨3|4+5|6+7|3⟩)/(⟨1|2⟩⟨2|3⟩⟨1|4+5|6+7|1⟩)
    +(+2⟨1|3⟩⟨3|4+5-6-7|2])/(⟨1|2⟩⟨1|4+5|6+7|1⟩)
    +('2134567', False, '-')
    +(⟨3|4+5-6-7|3]-⟨1|4+5-6-7|1]-⟨2|4+5-6-7|2])(1/2⟨1|3⟩⟨2|3⟩(⟨3|4+5-6-7|3]-⟨1|4+5-6-7|1]-⟨2|4+5-6-7|2])-1/2⟨1|2⟩⟨3|4+5-6-7|1-2|3⟩)/(⟨1|2⟩⟨1|4+5|6+7|1⟩⟨2|4+5|6+7|2⟩)
    """)
    coeff_v2 = Terms("""
    +(-8⟨1|3⟩³Δ_123|45|67)/(⟨1|2⟩⟨2|3⟩⟨1|4+5|6+7|1⟩²)
    +(-4⟨1|3⟩⟨3|4+5|6+7|3⟩)/(⟨1|2⟩⟨2|3⟩⟨1|4+5|6+7|1⟩)
    +(+2⟨1|3⟩⟨3|4+5-6-7|2])/(⟨1|2⟩⟨1|4+5|6+7|1⟩)
    +('2134567', False, '-')
    +tr(3-1-2|4+5-6-7)(1/2⟨1|3⟩⟨2|3⟩tr(3-1-2|4+5-6-7)-1/2⟨1|2⟩⟨3|4+5-6-7|1-2|3⟩)/(⟨1|2⟩⟨1|4+5|6+7|1⟩⟨2|4+5|6+7|2⟩)
    """)
    oPs = Particles(7, field=Field("finite field", 2 ** 31 - 1, 1))
    assert coeff_v1(oPs) - coeff_v2(oPs) == 0


def test_term_with_mass():
    string = """+(-7/3⟨2|4⟩⟨1|4|2]⟨3|4|1]+7/3⟨3|4⟩mt2²-14/3⟨3|4⟩⟨1|3|1]mt2+7/3⟨3|4⟩⟨1|3|1]²+7/3⟨2|4⟩⟨3|4|2]⟨1|4|1]-7/3⟨4|3|1]⟨1|3⟩mt2)/((s_123-mt2)⟨1|3|1]⟨1|3|2])"""
    assert len(str(Term(string))) == len(string)


def test_coeffs_normalization():
    oPs = Particles(6, field=Field("finite field", 2 ** 31 - 1, 1))
    oTerms = Terms("""
    +(+96/127[35]⟨4|2+3|1]⟨16⟩)/(⟨56⟩[56]⟨1|2+4|3]⟨2|1+4|3])
    +("432165", True, "+")
    +(-192/127[35]⟨26⟩⟨4|2+3|1])/(⟨56⟩[56]⟨2|1+4|3]^2)
    """)
    assert oTerms.with_normalized_coeffs(oPs) * oTerms.common_coeff_factor == oTerms(oPs)


def test_coeffs_normalization_big_and_outliers():
    from big_terms import oTerms
    assert oTerms.with_normalized_coeffs.max_denominator == 64
