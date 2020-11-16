
class DataFile:
    def __init__(self):
        self.pfad = None
        self.versuch = None
        self.df = None
        self.directory = None
        self.cnc = None
        self.slp_out = None

    def __getattr__(self, item):
        if item == "pfad":
            return self.pfad
        elif item == "versuch":
            return self.versuch
        elif item == "df":
            if self.df is not None:
                return self.df
            else:
                raise ValueError("Dataframe nicht eingestellt")
        elif item == "cnc":
            if self.cnc is not None:
                return self.cnc
            else:
                raise ValueError("cnc nicht eingestellt")
        elif item == "slp_out":
            if self.slp_out is not None:
                return self.slp_out
            else:
                raise ValueError("slp_out nicht eingestellt")
        else:
            raise ValueError("Wert nicht vorhanden")

    def set_data_frame(self, df):
        self.df = df

    def set_pfad(self, pfad):
        self.pfad = pfad

    def set_versuch(self, versuch):
        self.versuch = versuch

    def set_directory(self, directory):
        self.directory = directory

    def set_cnc(self, cnc):
        self.cnc = cnc

    def set_slp_out(self, slp_out):
        self.slp_out = slp_out
