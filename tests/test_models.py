import pytest

from mutalyzer_crossmapper.models import Point, CodingPoint, GenomicPoint, NonCodingPoint, Point, ProteinPoint

# Point tests
def test_point_valid_creation():
    p = Point(position=-10)
    assert p.position == -10
    assert p.offset == 0
    assert p.region == ''


def test_point_custom_values():
    p = Point(position=0, offset=4, region='*')
    assert p.position == 0
    assert p.offset == 4
    assert p.region == '*'


def test_genomic_point_valid_creation():
    p = GenomicPoint(position=10)
    assert p.position == 10


def test_genomic_point_invalid_creation():
    with pytest.raises(ValueError):
        GenomicPoint(position=-1)

    with pytest.raises(ValueError):
        GenomicPoint(position=1.0)

    with pytest.raises(ValueError):
        GenomicPoint(position='chr1')

    with pytest.raises(TypeError):
        GenomicPoint(**{})


def test_genomic_point_invalid_keyword():
    with pytest.raises(TypeError):
        GenomicPoint(position=1, region='')

    with pytest.raises(TypeError):
        GenomicPoint(**{'123': "abc"})


def test_genomic_point_to_string():
    assert str(GenomicPoint(123456789)) == "123456789"


# NonCodingPoint tests
def test_noncoding_point_default_values():
    p = NonCodingPoint(position=10)
    assert p.position == 10
    assert p.offset == 0
    assert p.region == ''


def test_noncoding_point_valid_creation():
    p = NonCodingPoint(position=10, offset=1, region='u')
    assert p.position == 10
    assert p.offset == 1
    assert p.region == 'u'


def test_noncoding_point_invalid_creation():
    with pytest.raises(ValueError):
        NonCodingPoint(position=0)

    with pytest.raises(ValueError):
        NonCodingPoint(position=10, region='-')

    with pytest.raises(ValueError):
        NonCodingPoint(position=10, region=123)

    with pytest.raises(TypeError):
        NonCodingPoint(position=1, offset='+1')

    with pytest.raises(TypeError):
        NonCodingPoint(offset=0)


def test_noncoding_point_invalid_keyword():
    with pytest.raises(TypeError):
        NonCodingPoint(position=10, other='test')


def test_noncoding_point_to_string():
    assert str(NonCodingPoint(position=123, offset=0)) == '123'
    assert str(NonCodingPoint(position=123, offset=11)) == '123+11'
    assert str(NonCodingPoint(position=123, region='u')) == 'u123'
    assert str(NonCodingPoint(position=123, offset=-10, region='')) == '123-10'
    assert str(NonCodingPoint(position=123, offset=-11, region='d')) == 'd123-11'


# CodingPoint tests
def test_coding_point_default_creation():
    p = CodingPoint(position=11)
    assert p.position == 11
    assert p.offset == 0
    assert p.region == ''


def test_coding_point_valid_creation():
    p = CodingPoint(position=987654321, offset=-1111, region='-')
    assert p.position == 987654321
    assert p.offset == -1111
    assert p.region == '-'

def test_coding_point_invalid_creation():
    with pytest.raises(ValueError):
        CodingPoint(position=0, offset=-1, region='')


def test_coding_point_to_string():
    assert str(CodingPoint(position=123, offset=0)) == '123'
    assert str(CodingPoint(position=123, offset=11)) == '123+11'
    assert str(CodingPoint(position=123, offset=11, region='-')) == '-123+11'
    assert str(CodingPoint(position=123, offset=11, region='*')) == '*123+11'
    assert str(CodingPoint(position=123, offset=-11, region='-')) == '-123-11'
    assert str(CodingPoint(position=123, region='u')) == 'u123'
    assert str(CodingPoint(position=123, offset=-10, region='')) == '123-10'
    assert str(CodingPoint(position=123, offset=-11, region='d')) == 'd123-11'


# ProteinPoint tests
def test_protein_point_default_position_in_codon():
    p = ProteinPoint(position=10)
    assert p.position == 10
    assert p.region == ''
    assert p.offset == 0
    assert p.position_in_codon == 1


def test_protein_point_valid_creation():
    p = ProteinPoint(position=10, position_in_codon=2)
    assert p.position == 10
    assert p.region == ''
    assert p.offset == 0
    assert p.position_in_codon == 2


def test_protein_point_invalid_creation():
    with pytest.raises(ValueError):
        ProteinPoint(position=10, position_in_codon=0)


def test_protein_point_to_string():
    assert str(ProteinPoint(position=11)) == '11'
    assert str(ProteinPoint(position=11, offset=0, region='', position_in_codon=2)) == '11'
    assert str(ProteinPoint(position=11, offset=1, region='*')) == '*11+1'
    assert str(ProteinPoint(position=11, offset=-1, region='-', position_in_codon=3)) == '-11-1'
    assert str(ProteinPoint(position=11, offset=0, region='d')) == 'd11'
    assert str(ProteinPoint(position=11, region='u', offset=10)) == 'u11+10'
