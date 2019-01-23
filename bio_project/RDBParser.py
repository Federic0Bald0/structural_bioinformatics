class RDBParser:
    '''
    Class representing a RepeatsDB file parsed
    '''

    source = ""
    pdb = ""
    chain = ""
    regions = []
    units = []
    insertions = []

    def parse_rdb(self, filename):
        '''
        Parse a RepeatsDB file
        :param filename: the file source
        :return: RDBParser object
        '''
        with open(filename) as myfile:
            for line in myfile:
                line = line.strip().split()
                if line[0] == "SOURCE":
                    self.source = line[1]
                elif line[0] == "PDB":
                    self.pdb = line[1]
                elif line[0] == "CHAIN":
                    self.chain = line[1]
                elif line[0] == "REG":
                    self.regions.append(line[1])
                    self.regions.append(line[2])
                    self.regions.append(line[3] + '.' + line[4])
                elif line[0] == "UNIT":
                    start_r = int(line[1])
                    end_r = int(line[2])
                    self.units.append((start_r, end_r))
                elif line[0] == "INS":
                    self.insertions.append((line[1], line[2]))
        return self

