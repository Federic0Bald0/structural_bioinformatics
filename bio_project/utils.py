# coding: utf-8


def parse_rdb(filename):
    """
    Parse RepeadDB protein file 
    """
    obj = {}
    with open(filename) as f:
        for line in f:
            line = line.strip().split()
            if line[0] == "SOURCE":
                obj['source'] = line[1]
            elif line[0] == "PDB":
                obj['id'] = line[1]
            elif line[0] == "CHAIN":
                obj['chain'] = line[1]
            elif line[0] == "REG":
                obj.setdefault('regions', [])
                obj['regions'].append((line[1], line[2], 
                                       line[3] + '.' + line[4]))
            elif line[0] == "UNIT":
                start, end = line[1:3]
                try:
                    start_id = (' ', int(start), ' ')
                except:
                    start_id = (' ', int(start[:-1]), start[-1])
                try:
                    end_id = (' ', int(end), ' ')
                except:
                    end_id = (' ', int(end[:-1]), end[-1])

                obj.setdefault('units', [])
                obj['units'].append((start_id, end_id))

            elif line[0] == "INS":
                obj.setdefault('insertions', [])
                obj['insertions'].append((line[1], line[2]))

    for i, unit in enumerate(obj['units']):
        print "Unit {:<5} {}".format(i, unit)

    # for i, unit in enumerate(obj['insertions']):
    #     print "Insertion {:<5} {}".format(i, unit)

    return obj


# if __name__ == "__main__":
#     db_f = 'data/rdb/2z7xb.db'
#     print parse_rdb(db_f)