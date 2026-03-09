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
        """Convert a coordinate to a proper position model.

        :arg int coordinate: Coordinate.

        :returns dict: Position model with 'position' and 'offset' keys.
        """
        if self._inverted:
            if coordinate > self.boundary[1]:
                return {"position": 0, "offset": self.boundary[1] - coordinate}
            if coordinate < self.boundary[0]:
                return {"position": self._end, "offset": self.boundary[0] - coordinate}
            return {"position": self.boundary[1] - coordinate, "offset": 0}

        if coordinate < self.boundary[0]:
            return {"position": 0, "offset": coordinate - self.boundary[0]}
        if coordinate > self.boundary[1]:
            return {"position": self._end, "offset": coordinate - self.boundary[1]}
        return {"position": coordinate - self.boundary[0], "offset": 0}

    def to_coordinate(self, pos_m):
        """Convert a position model to a coordinate.

        :arg dict position: Position model with 'position' and 'offset' keys.

        :returns int: Coordinate.
        """
        if self._inverted:
            return self.boundary[1] - pos_m["position"] - pos_m["offset"]
        return self.boundary[0] + pos_m["position"] + pos_m["offset"]
