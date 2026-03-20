"""a script to check degenerate option"""
from mutalyzer_crossmapper import Coding, Genomic, NonCoding


def serialize(pos_m: dict):
    if pos_m["offset"] > 0:
        return f"{pos_m['region']}{pos_m['position']}+{pos_m['offset']}"
    elif pos_m["offset"] < 0:
        return f"{pos_m['region']}{pos_m['position']}{pos_m['offset']}"
    else:
        return f"{pos_m['region']}{pos_m['position']}"

_exons = [(5, 8), (14, 20), (30, 35), (40, 44), (50, 52), (70, 72)]
_cds = (32, 43)


nc_crossmap = NonCoding(_exons, True)
c_crossmap = Coding(_exons, _cds)
test = Coding([(10, 11)], (10, 11))


for i in range(0, 80):
    # nc = nc_crossmap.coordinate_to_noncoding(i)
    c = c_crossmap.coordinate_to_protein(i)

    # print(i, c, c_crossmap.coding_to_coordinate(c))
    # c_de = nc_crossmap.coordinate_to_coding(i, True)
    # nc_de = nc_crossmap.coordinate_to_coding(i, True)
    print(i, c)
    # print(f'"{i}", "{c["position"]}", "{c["position_in_codon"]}", "{c["offset"]}", "{c["region"]}"')


# crossmap = Coding(_exons, _cds)
# for i in range(0, 80):
#     print(i,  crossmap.coordinate_to_coding(i), crossmap.coordinate_to_coding(i, degenerate=True))

# nc_crossmap = NonCoding(_exons)
# for i in range(0, 80):
#     print(i, nc_crossmap.coordinate_to_noncoding(i))






# degereate option
"""
With this option, it keeps counting c_pos outside the exons range
e.g., (-16, 0, -1, -5) means
c_pos=16, offset_to_c_pos=0, before_CDS, offset_to_exons_range = -5
"""
