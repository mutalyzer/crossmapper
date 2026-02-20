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

        if coordinate < self.boundary[0]:           # upstream of an exon, re
            return {"position": 0, "offset": coordinate - self.boundary[0]}
        if coordinate > self.boundary[1]:           # downstream of an exon
            return {"position": self._end, "offset": coordinate - self.boundary[1]}
        return {"position": coordinate - self.boundary[0], "offset": 0}

    def to_coordinate(self, position_m):
        """Convert a position model to a coordinate.

        :arg dict position: Position model with 'position' and 'offset' keys.

        :returns int: Coordinate.
        """
        if self._inverted:
            if position_m["region"] == "u":
                return self.boundary[1] + position_m["position"]
            elif position_m["region"] == "d":
                return self.boundary[0] - position_m["position"]
            return self.boundary[1] - position_m["position"] - position_m["offset"]
        if position_m["region"] == "u":
            return self.boundary[0] - position_m["position"]
        elif position_m["region"] == "d":
            return self.boundary[1] + position_m["position"]
        else:
            return self.boundary[0] + position_m["position"] + position_m["offset"]
