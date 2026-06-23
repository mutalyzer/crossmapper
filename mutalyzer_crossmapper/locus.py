from .models import Point


class Locus(object):
    """Locus object."""
    def __init__(self, location: tuple[int, int], inverted: bool = False) -> None:
        """
        :arg tuple location: Locus location.
        :arg bool inverted: Orientation.
        """
        self._inverted = inverted

        self.boundary = location[0], location[1] - 1
        self._end = self.boundary[1] - self.boundary[0]

    def to_position(self, coordinate: int) -> Point:
        """Convert a coordinate to a proper point model.

        :arg int coordinate: Coordinate.

        :returns Point: Position point model.
        """
        if self._inverted:
            if coordinate > self.boundary[1]:
                return Point(position=0, offset=self.boundary[1] - coordinate)
            if coordinate < self.boundary[0]:
                return Point(position=self._end, offset=self.boundary[0] - coordinate)
            return Point(position=self.boundary[1] - coordinate, offset=0)

        if coordinate < self.boundary[0]:
            return Point(position=0, offset=coordinate - self.boundary[0])
        if coordinate > self.boundary[1]:
            return Point(position=self._end, offset=coordinate - self.boundary[1])
        return Point(position=coordinate - self.boundary[0], offset=0)

    def to_coordinate(self, point: Point) -> int:
        """Convert a point model to a coordinate.

        :arg Point point: Point model.

        :returns int: Coordinate.
        """
        if self._inverted:
            return self.boundary[1] - point.position - point.offset
        return self.boundary[0] + point.position + point.offset
