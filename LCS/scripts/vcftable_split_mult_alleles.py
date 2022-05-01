
import sys

h = open(sys.argv[1],"rt") if len(sys.argv)>1 else sys.stdin

headers = h.readline().rstrip().split("\t")

def filter_ann_fields(ann):
    ann_fs = ann.split('|')
    fs = map(lambda i: ann_fs[i], [1,2,3,6,9,10])
    return "|".join(fs)

def normalize_ref_alt(ref, alt):
    if ref[0] != alt[0] and ref[1:] == alt[1:]:
        return ref[0], alt[0]
    else:
        return ref, alt

# remove sample name from AD and AF
headers = list(map(lambda i: "AD" if "AD" in i else "AF" if "AF" in i else i, headers))
print("\t".join(headers))
for l in h:
    l = l.rstrip()
    rows = l.split("\t")
    if "," in rows[headers.index("ALT")]:
        ref = rows[headers.index("REF")]
        alts = rows[headers.index("ALT")].split(",")
        anns = rows[headers.index("ANN")].split(",")
        ads = rows[headers.index("AD")].split(",")
        afs = rows[headers.index("AF")].split(",")
        for i in range(0,len(alts)):
            alt = alts[i]
            rows[headers.index("ALT")] = alt
            if len(ref) == len(alt) and len(ref)>1:
                ref2, alt2 = normalize_ref_alt(ref, alt)
                rows[headers.index("REF")] = ref2
                rows[headers.index("ALT")] = alt2
            x1 = filter(lambda s: s.startswith(alts[i]), anns)
            x2 = map(filter_ann_fields, x1)
            rows[headers.index("ANN")] = ",".join(x2)
            rows[headers.index("AD")] = ads[0]+","+ads[i+1]
            rows[headers.index("AF")] = afs[i]
            print("\t".join(rows))
    else:
        anns = rows[headers.index("ANN")].split(",")
        x2 = map(filter_ann_fields, anns)
        rows[headers.index("ANN")] = ",".join(x2)
        print("\t".join(rows))
