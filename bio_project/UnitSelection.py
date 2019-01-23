from Bio.PDB import Select


class UnitSelection(Select):
    '''
    Class used to select a unit starting at start_r and finishing at end_r
    '''
    def __init__(self, start_r, end_r):
        '''
        Constructor
        :param start_r: from where the unit starts
        :param end_r: the unit end
        '''
        self.start_r = start_r
        self.end_r = end_r

    def accept_residue(self, residue):
        '''
        Select only the residue in the interval
        :param residue: the residue
        :return: True if the residue is in the interval, false otherwise
        '''
        hetatm_flag, resseq, icode = residue.get_id()
        if self.start_r <= resseq <= self.end_r:
            return True
        return False