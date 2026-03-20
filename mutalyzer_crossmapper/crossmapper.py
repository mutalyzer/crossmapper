from .multi_locus import MultiLocus


class Genomic(object):
    """Genomic crossmap object."""
    def coordinate_to_genomic(self, coordinate):
        """Convert a coordinate to a genomic position (g./m./o.).

        :arg int coordinate: Coordinate.

        :returns dict: Genomic position model.
        """
        return {"position": coordinate + 1}

    def genomic_to_coordinate(self, pos_m):
        """Convert a genomic position (g./m./o.) to a coordinate.

        :arg dict pos_m: Genomic position model.

        :returns int: Coordinate.
        """
        return pos_m["position"] - 1


class NonCoding(Genomic):
    """NonCoding crossmap object."""
    def __init__(self, locations, inverted=False):
        """
        :arg list locations: List of locus locations.
        :arg bool inverted: Orientation.
        """
        self._inverted = inverted

        self._noncoding = MultiLocus(locations, inverted)

    def coordinate_to_noncoding(self, coordinate):
        """Convert a coordinate to a noncoding position (n./r.).

        :arg int coordinate: Coordinate.

        :returns dict: Noncoding position model.
        """
        pos_m = self._noncoding.to_position(coordinate)
        if pos_m["region"] == "":
            pos_m["position"] = pos_m["position"] + 1
        return pos_m

    def noncoding_to_coordinate(self, pos_m):
        """Convert a noncoding position (n./r.) to a coordinate.

        :arg dict pos_m: Noncoding position model.

        :returns int: Coordinate.
        """
        multilocus_pos_m = {**pos_m}
        if pos_m["region"] == "":
            multilocus_pos_m["position"] = pos_m["position"] - 1
        return self._noncoding.to_coordinate(multilocus_pos_m)


class Coding(NonCoding):
    """Coding crossmap object."""
    def __init__(self, locations, cds, inverted=False):
        """
        :arg list locations: List of locus locations.
        :arg tuple cds: Locus location.
        :arg bool inverted: Orientation.
        """
        NonCoding.__init__(self, locations, inverted)

        b0 = self._noncoding.to_position(cds[0])
        b1 = self._noncoding.to_position(cds[1]-1)
        e0 = self._noncoding.to_position(locations[0][0])
        e1 = self._noncoding.to_position(locations[-1][1]-1)

        if self._inverted:
            self._coding = (b1["position"] + b1["offset"], b0["position"] + b0["offset"] + 1)
            self._exons = (e1["position"], e0["position"])
        else:
            self._coding = (b0["position"] + b0["offset"], b1["position"] + b1["offset"] +1)
            self._exons = (e0["position"], e1["position"])

    def _degenerate_position(self, pos_m):
        """Degenerate a coding position model (c./r.).

        :arg dict pos_m: Coding position model.

        :returns dict: a generate coding position model.
        """
        region = pos_m["region"]
        position = pos_m["position"]

        degenerated_pos_m = {"offset": pos_m["offset"]}

        if region == "u":
            if self._inverted:
                degenerated_pos_m["position"] = position + self._exons[1] - self._coding[1] + 1
            else:
                degenerated_pos_m["position"] = position + self._coding[0]
            degenerated_pos_m["region"] = "-"
        if region == "d":
            if self._inverted:
                degenerated_pos_m["position"] = position + self._coding[0]
            else:
                degenerated_pos_m["position"] = position + self._exons[1]- self._coding[1] + 1
            degenerated_pos_m["region"] = "*"
        return degenerated_pos_m

    def _normalize_position(self, pos_m):
        """Normalize a coding position model (c./r.).

        :arg dict pos_m: Coding position model.

        :returns dict: a normalized coding postion model.
        """
        initial_pos = {**pos_m, "offset": 0}
        coordinate = self._coding_to_coordinate(initial_pos)
        if self._inverted:
            coordinate = coordinate - pos_m["offset"]
        else:
            coordinate = coordinate + pos_m["offset"]
        return self.coordinate_to_coding(coordinate)

    def _coordinate_to_coding(self, coordinate):
        """Convert a coordinate to a coding position (c./r.).

        :arg int coordinate: Coordinate.

        :returns dict: Coding position model (c./r.).
        """
        noncoding_pos_m = self._noncoding.to_position(coordinate)

        if noncoding_pos_m["region"] in ["u", "d"]:
            return noncoding_pos_m

        location = noncoding_pos_m["position"]
        offset = noncoding_pos_m["offset"]
        if location < self._coding[0]:
            return {
                "position": self._coding[0] - location,
                "offset": offset,
                "region": "-"
            }
        if location >= self._coding[1]:
            return {
                "position": location - self._coding[1] + 1,
                "offset": offset,
                "region": "*"
            }
        return {
            "position": location - self._coding[0] + 1,
            "offset": offset,
            "region": ""
        }

    def coordinate_to_coding(self, coordinate, degenerate=False):
        """Convert a coordinate to a coding position (c./r.).

        :arg int coordinate: Coordinate.
        :arg bool degenerate: Return a degenerate position.

        :returns dict: Coding position model (c./r.).
        """
        pos_m = self._coordinate_to_coding(coordinate)

        if degenerate and pos_m["region"] in ("u", "d"):
            pos_m = self._degenerate_position(pos_m)

        return pos_m

    def _coding_to_coordinate(self, pos_m):
        """Convert a coding position (c./r.) to a coordinate.

        :arg dict pos_m: Coding position model (c./r.).

        :returns int: Coordinate.
        """
        position = pos_m["position"]
        region = pos_m["region"]

        if region in ["u", "d"]:
            return self._noncoding.to_coordinate(pos_m)

        noncoding_pos_m = {"offset": pos_m["offset"], "region": ""}
        if region == "":
            noncoding_pos_m["position"] = position + self._coding[0] - 1
        elif region == "-":
            noncoding_pos_m["position"] = self._coding[0] - position
        else:
            noncoding_pos_m["position"] = self._coding[1] + position - 1

        return self._noncoding.to_coordinate(noncoding_pos_m)

    def coding_to_coordinate(self, pos_m):
        """Convert a coding position (c./r.) to a coordinate.

        :arg dict pos_m: Coding position model (c./r.).

        :returns int: Coordinate.
        """
        normalized_pos_m = self._normalize_position(pos_m)

        return self._coding_to_coordinate(normalized_pos_m)

    def coordinate_to_protein(self, coordinate):
        """Convert a coordinate to a protein position (p.).

        :arg int coordinate: Coordinate.

        :returns dict: Protein position model(p.).
        """
        pos = self.coordinate_to_coding(coordinate)

        if pos["region"] == "u":
            pos = self.coordinate_to_coding(coordinate + pos["position"])
        elif pos["region"] == "d":
            pos = self.coordinate_to_coding(coordinate - pos["position"])

        position = pos["position"]
        if pos["region"] == "-":
            return {
                "position": abs(-position // 3),
                "position_in_codon": -position % 3 + 1,
                **{k: v for k, v in pos.items() if k != "position"}}
        return {
                "position": (position + 2) // 3,
                "position_in_codon": (position + 2) % 3 + 1,
                **{k: v for k, v in pos.items() if k != "position"}}

    def protein_to_coordinate(self, pos_m):
        """Convert a protein position (p.) to a coordinate.

        :arg dict position: Protein position model(p.).

        :returns int: Coordinate.
        """
        if pos_m["region"] == "-":
            return self.coding_to_coordinate(
                {"position": 3 * pos_m["position"] - pos_m["position_in_codon"] + 1,
                 "offset": pos_m["offset"],
                 "region": pos_m["region"]})

        return self.coding_to_coordinate(
                {"position": 3 * pos_m["position"] + pos_m["position_in_codon"] - 3,
                 "offset": pos_m["offset"],
                 "region": pos_m["region"]})
