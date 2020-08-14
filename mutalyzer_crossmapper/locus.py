class Locus(object):
    """Locus object."""
    def __init__(self, location, inverted=False):
        """
        :arg tuple location: Locus location.
        :arg bool inverted: Orientation.
        """
        self._inverted = inverted

        self.boundary = location[0], location[1] - 1
        self._end = self.boundary[1] - self.boundary[0]

    def to_position(self, coordinate):
        """Convert a coordinate to a proper position.

        :arg int coordinate: Coordinate.

        :returns tuple: Position.
        """
        if self._inverted:
            if coordinate > self.boundary[1]:
                return 0, self.boundary[1] - coordinate
            if coordinate < self.boundary[0]:
                return self._end, self.boundary[0] - coordinate
            return self.boundary[1] - coordinate, 0

        if coordinate < self.boundary[0]:
            return 0, coordinate - self.boundary[0]
        if coordinate > self.boundary[1]:
            return self._end, coordinate - self.boundary[1]
        return coordinate - self.boundary[0], 0

    #def to_position(self, coordinate, degenerate=False):
    #    """Convert a coordinate to a position.

    #    :arg int coordinate: Coordinate.
    #    :arg bool degenerate: Return a degenerate position.

    #    :returns tuple: Position.
    #    """
    #    position = self._to_position(coordinate)

    #    if degenerate:
    #        if position[0] == 1:
    #            return position[1], 0
    #        return position[0] + position[1], 0

    #    return position

    def to_coordinate(self, position):
        """Convert a position to a coordinate.

        :arg int position: Position.

        :returns int: Coordinate.
        """
        if self._inverted:
            #if position[0] < 1:  # Degenerate position.
            #    return self._location[1] - position[0] - position[1] - 1
            return self.boundary[1] - position[0] - position[1]

        #if position[0] < 1:      # Degenerate position.
        #    return self._location[0] + position[0] + position[1]
        return self.boundary[0] + position[0] + position[1]
