import os


class CDFfile(object):

    def __init__(self, fname=None, read=True):

        self.fname = fname

        if read and self.fname:
            self.read_nc()

    def read_nc(self, fname=None):

        if not os.path.exists(fname):
            raise OSError('File not found : {}'.format(fname))
