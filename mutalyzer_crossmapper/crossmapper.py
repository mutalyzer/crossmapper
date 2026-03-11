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
        if pos_m["region"] == "":
            pos_m["position"] = pos_m["position"] - 1
        return self._noncoding.to_coordinate(pos_m)


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
            self._cds_len = abs((b0["position"] + b0["offset"]) - (b1["position"] + b1["offset"]))
            self._exons_end = e1["position"]
            self._exons_start = e0["position"]
        else:
            self._coding = (b0["position"] + b0["offset"], b1["position"] + b1["offset"] +1)
            self._cds_len = abs((b1["position"] + b1["offset"]) - (b0["position"] + b0["offset"]))
            self._exons_end = e1["position"]
            self._exons_start = e0["position"]

    def _coordinate_to_coding(self, coordinate):
        """Convert a coordinate to a coding position (c./r.).

        :arg int coordinate: Coordinate.

        :returns dict: Coding position model (c./r.).
        """
        noncoding_pos_m = self._noncoding.to_position(coordinate)
        location = noncoding_pos_m["position"]
        if noncoding_pos_m["region"] == "":
            if location < self._coding[0]:
                return {
                    "position": self._coding[0] - location,
                    "offset": noncoding_pos_m["offset"],
                    "region": "-"
                }
            elif location >= self._coding[1]:
                return {
                    "position": location - self._coding[1] + 1,
                    "offset": noncoding_pos_m["offset"],
                    "region": "*"
                }
            else:
                return {
                    "position": location - self._coding[0] + 1,
                    "offset": noncoding_pos_m["offset"],
                    "region": ""
                }
        else:
            return noncoding_pos_m

    def coordinate_to_coding(self, coordinate, degenerate=False):
        """Convert a coordinate to a coding position (c./r.).

        :arg int coordinate: Coordinate.
        :arg bool degenerate: Return a degenerate position.

        :returns dict: Coding position model (c./r.).
        """
        pos_m = self._coordinate_to_coding(coordinate)
        if degenerate:
            if pos_m["region"] == "u":
                if self._inverted:
                    pos_m["position"] = pos_m["position"] + self._exons_start - self._coding[1] + 1
                else:
                    pos_m["position"] = pos_m["position"] + self._coding[0]
                pos_m["region"] = "-"
            if pos_m["region"] == "d":
                if self._inverted:
                    pos_m["position"] = pos_m["position"] + self._coding[0]
                else:
                    pos_m["position"] = pos_m["position"] + self._exons_end - self._coding[1] + 1
                pos_m["region"] = "*"
        return pos_m

    def coding_to_coordinate(self, pos_m):
        """Convert a coding position (c./r.) to a coordinate.

        :arg dict pos_m: Coding position model (c./r.).

        :returns int: Coordinate.
        """
        region = pos_m["region"]
        position = pos_m["position"]
        offset = pos_m["offset"]

        if region == "u":
            noncoding_pos_m = {
                "position": position - offset,
                "offset": 0,
                "region": "u"
            }
        elif region == "d":
            noncoding_pos_m = {
                "position": position + offset,
                "offset": 0,
                "region": "d"
            }
        elif region == "":
            noncoding_pos_m = {
                "position": position + self._coding[0] -1,
                "offset": offset,
                "region": ""
            }
        # add checks for degenerate results?
        elif region == "-":
            if position > self._coding[0]:
                noncoding_pos_m = {
                    "position": position - self._coding[0],
                    "offset": offset,
                    "region": "u"
                }
            else:
                noncoding_pos_m = {
                    "position": self._coding[0] - position,
                    "offset": offset,
                    "region": ""
                }
        else: # *
            if position > self._coding[1]:
                noncoding_pos_m = {
                    "position": position - self._coding[1],
                    "offset": offset,
                    "region": "d"
                }
            else:
                noncoding_pos_m = {
                    "position": self._coding[1] + position - 1,
                    "offset": offset,
                    "region": ""
                }
        return self._noncoding.to_coordinate(noncoding_pos_m)


    def coordinate_to_protein(self, coordinate):
        """Convert a coordinate to a protein position (p.).

        :arg int coordinate: Coordinate.

        :returns tuple: Protein position (p.).
        """
        pos = self.coordinate_to_coding(coordinate)

        if pos["region"] in ["-", "*"]:
            return {
                "position": pos["position"] // 3 + 1,
                "position_in_codon": pos["position"] % 3,
                **{k: v for k, v in pos.items() if k != "position"}}
        return {
                "position": (pos["position"]+2) // 3,
                "position_in_codon": (pos["position"]+2) % 3 + 1,
                **{k: v for k, v in pos.items() if k != "position"}}

    def protein_to_coordinate(self, position):
        """Convert a protein position (p.) to a coordinate.

        :arg tuple position: Protein position (p.).

        :returns int: Coordinate.
        """
        if position["region"] in ["-", "*"]:
            return self.coding_to_coordinate(
                {"position": 3 * position["position"] + position["position_in_codon"] - 3,
                 "offset": position["offset"],
                 "region": position["region"]})

        return self.coding_to_coordinate(
                {"position": 3 * position["position"] + position["position_in_codon"] - 3,
                 "offset": position["offset"],
                 "region": position["region"]})
