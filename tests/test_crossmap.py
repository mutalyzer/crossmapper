from crossmapper import Crossmap


class TestCrossmap(object):
    def setup(self):
        self._crossmap = Crossmap({
            'features': [
                {
                    'start': {
                        'position': 9},
                    'end': {
                        'position': 19},
                    'type': 'exon'},
                {
                    'start': {
                        'position': 29},
                    'end': {
                        'position': 49},
                    'type': 'exon'},
                {
                    'start': {
                        'position': 59},
                    'end': {
                        'position': 69},
                    'type': 'exon'},
                {
                    'start': {
                        'position': 79},
                    'end': {
                        'position': 89},
                    'type': 'exon'},
                {
                    'start': {
                        'position': 39},
                    'end': {
                        'position': 84},
                    'type': 'cds'}],
            'start': {
                'position': 9},
            'end': {
                'position': 89}})

    def _test_position(self, f, f_inv, coordinate, position):
        """Test coordinate to position conversion and its inverse."""
        assert f(coordinate) == position
        assert f_inv(position) == coordinate

    def test_variables(self):
        assert not self._crossmap._inverted
        assert self._crossmap._locus == (9, 89)
        assert self._crossmap._exons == [(9, 19), (29, 49), (59, 69), (79, 89)]
        assert self._crossmap._cds == (39, 84)

    def test_genomic(self):
        self._test_position(
            self._crossmap.coordinate_to_genomic,
            self._crossmap.genomic_to_coordinate,
            {'position': 49}, {'position': 50})

    def test_locus(self):
        self._test_position(
            self._crossmap.coordinate_to_locus,
            self._crossmap.locus_to_coordinate,
            {'position': 49}, {
                'position': 10,
                'offset': {
                    'position': 40}})
        self._test_position(
            self._crossmap.coordinate_to_locus,
            self._crossmap.locus_to_coordinate,
            {'position': 4}, {
                'position': 10,
                'offset': {
                    'position': -5}})
        self._test_position(
            self._crossmap.coordinate_to_locus,
            self._crossmap.locus_to_coordinate,
            {'position': 99}, {
                'position': 10,
                'offset': {
                    'position': 90}})

    def test_locus_inverted(self):
        self._crossmap._inverted = True
        self._test_position(
            self._crossmap.coordinate_to_locus,
            self._crossmap.locus_to_coordinate,
            {'position': 49}, {
                'position': 90,
                'offset': {
                    'position': 40}})
        self._test_position(
            self._crossmap.coordinate_to_locus,
            self._crossmap.locus_to_coordinate,
            {'position': 4}, {
                'position': 90,
                'offset': {
                    'position': 85}})
        self._test_position(
            self._crossmap.coordinate_to_locus,
            self._crossmap.locus_to_coordinate,
            {'position': 99}, {
                'position': 90,
                'offset': {
                    'position': -10}})

    def test_noncoding(self):
        pass
