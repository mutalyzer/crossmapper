from bisect import bisect_right
from itertools import accumulate

from .location import nearest_location
from .locus import Locus


def _offsets(locations, orientation):
    """For each location, calculate the length of the preceding locations.

    :arg list locations: List of locations.
    :arg int orientation: Direction of {locations}.

    :returns list: List of cumulative location lengths.
    """
    return [0] + list(accumulate(map(
        lambda x: x[1] - x[0], locations[::orientation][:-1])))


class MultiLocus(object):
    """MultiLocus object."""
    def __init__(self, locations:list, inverted=False):
        """
        :arg list locations: List of locus locations.
        :arg bool inverted: Orientation.
        """
        self._locations = locations
        self._inverted = inverted

        self._loci = [Locus(location, inverted) for location in locations]
        self._orientation = -1 if inverted else 1
        self._offsets = _offsets(locations, self._orientation)

    def _direction(self, index):
        if self._inverted:
            return len(self._offsets) - index - 1
        return index

    def outside(self, coordinate:int):
        """Calculate the offset relative to this MultiLocus.

        :arg int coordinate: Coordinate.

        :returns int: Negative: upstream, 0: inside, positive: downstream.
        """
        if coordinate < self._loci[0].boundary[0]:
            return coordinate - self._loci[0].boundary[0]
        if coordinate > self._loci[-1].boundary[1]:
            return coordinate - self._loci[-1].boundary[1]
        return 0

    def to_position(self, coordinate:int):
        """Convert a coordinate to a position.

        :arg int coordinate: Coordinate.

        :returns dict: Position model 'position', 'offset' and 'region' keys.
        """
        index = nearest_location(self._locations, coordinate, self._inverted)
        outside = self._orientation * self.outside(coordinate)
        region = "u" if outside < 0 else "d" if outside > 0 else ""
        location = self._loci[index].to_position(coordinate)

        if outside:
            return {
                "position": abs(location["offset"]),
                "offset": 0,
                "region": region
            }

        if location["offset"] == 0:
            return {
                "position": location["position"] + self._offsets[self._direction(index)],
                "offset": 0,
                "region": ""
            }

        if location["offset"] < 0:
            return {
                "position": self._offsets[self._direction(index)],
                "offset": location["offset"],
                "region": ""
            }

        return{
            "position": location["position"] + self._offsets[self._direction(index)],
            "offset": location["offset"],
            "region": ""
        }

    def to_coordinate(self, pos_m:dict):
        """Convert a position model to a coordinate.

        :arg dict position: Position model with 'position','offset' and 'region' keys.

        :returns int: Coordinate.
        """
        region = pos_m["region"]

        if region == "u":
            if self._inverted:
                return abs(pos_m["position"]) + self._locations[-1][1] + pos_m["offset"] - 1
            return self._locations[0][0] - abs(pos_m["position"]) + pos_m["offset"]

        if region == "d":
            if self._inverted:
                return self._locations[0][0] - abs(pos_m["position"]) + pos_m["offset"]
            return abs(pos_m["position"]) + self._locations[-1][1] + pos_m["offset"] - 1

        index = min(
            len(self._offsets),
            max(0, bisect_right(self._offsets, pos_m["position"]) - 1)
        )
        pos_m["position"] = pos_m["position"] - self._offsets[index]
        return self._loci[self._direction(index)].to_coordinate(pos_m)
