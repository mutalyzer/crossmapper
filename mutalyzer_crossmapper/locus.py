class Locus(object):
    """Locus object."""
    def __init__(self, location, inverted=False):
        """
        :arg tuple location: Locus location.
        :arg bool inverted: Orientation.
        """
        self._location = location
        self._inverted = inverted

    def _to_position(self, coordinate):
        """Convert a coordinate to a proper position.

        :arg int coordinate: Coordinate.

        :returns tuple: Position.
        """
        if self._inverted:
            if coordinate >= self._location[1]:
                return 1, self._location[1] - coordinate - 1
            if coordinate < self._location[0]:
                return (
                    self._location[1] - self._location[0],
                    self._location[0] - coordinate)
            return self._location[1] - coordinate, 0

        if coordinate < self._location[0]:
            return 1, coordinate - self._location[0]
        if coordinate >= self._location[1]:
            return (
                self._location[1] - self._location[0],
                coordinate - self._location[1] + 1)
        return coordinate - self._location[0] + 1, 0

    def to_position(self, coordinate, degenerate=False):
        """Convert a coordinate to a position.

        :arg int coordinate: Coordinate.
        :arg bool degenerate: Return a degenerate position.

        :returns tuple: Position.
        """
        position = self._to_position(coordinate)

        if degenerate:
            if position[0] == 1:
                return position[1], 0
            return position[0] + position[1], 0

        return position

    def to_coordinate(self, position):
        """Convert a position to a coordinate.

        :arg int position: Position.

        :returns int: Coordinate.
        """
        if self._inverted:
            if position[0] < 1:  # Degenerate position.
                return self._location[1] - position[0] - position[1] - 1
            return self._location[1] - position[0] - position[1]

        if position[0] < 1:      # Degenerate position.
            return self._location[0] + position[0] + position[1]
        return self._location[0] + position[0] + position[1] - 1
