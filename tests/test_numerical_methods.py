import pytest
import numpy

from lips.fields.field import Field

from antares.core.settings import settings
from antares.core.numerical_methods import num_func, tensor_function


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def func(oPs):
    """Scalar function with crazy powers, to stress test."""
    return oPs('1/(⟨1|2⟩^10⟨1|3⟩³⟨2|4⟩²⟨1|(3+4)|1]⟨2|(3+4)|2]s_134s_234tr5_1234²[2|4]^20)')


def tensor_func(oPs):
    """Tensor function with crazy powers and constants (including zero)."""
    return numpy.array([[oPs.field(0), oPs('1/(⟨1|2⟩^10⟨1|3⟩³⟨2|4⟩²⟨1|(3+4)|1]⟨2|(3+4)|2]s_134s_234tr5_1234²[2|4]^20)')], [oPs("s12"), oPs.field(1)]])


oTF = tensor_function(tensor_func)
oTF.multiplicity = 6
func.multiplicity = 6
oF = num_func(func)


@pytest.mark.parametrize("field", [
    Field("mpc", 0, 300),
    Field("padic", 2 ** 31 - 19, 5),
    Field("finite field", 2 ** 31 - 1, 1),
])
class TestParametrizedOverFields:

    @staticmethod
    def test_scalar_mass_dimension(field):
        settings.field = field
        assert oF.mass_dimension == -51

    @staticmethod
    def test_tensor_mass_dimension(field):
        settings.field = field
        assert numpy.all(oTF.mass_dimension == numpy.array([[0, -51], [2, 0]]))

    @staticmethod
    def test_scalar_phase_weights(field):
        settings.field = field
        assert oF.phase_weights == [-13, 8, -3, 18, 0, 0]

    @staticmethod
    def test_tensor_phase_weights(field):
        settings.field = field
        assert numpy.all(oTF.phase_weights == numpy.array([[[0, 0, 0, 0, 0, 0], [-13, 8, -3, 18, 0, 0]], [[0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0]]]))


def test_tensor_function_len_and_shape():
    oTF = tensor_function(lambda x: numpy.array([[x, x ** 2], [2 * x, -x]]))
    with pytest.raises(AttributeError):
        len(oTF)
    with pytest.raises(AttributeError):
        oTF.shape
    oTF(1)  # evaluating it once should set the properties
    assert len(oTF) == 2
    assert oTF.shape == (2, 2)
