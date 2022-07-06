class Annotation:

    def __init__(self,id,inchikey,adduct,ion_mode,english_name,is_istd,is_known):
        self.id=id
        self.inchikey=inchikey
        self.adduct=adduct
        self.ion_mode=ion_mode
        self.english_name=english_name
        self.is_istd=is_istd
        self.is_known=is_known

    def set_comment(self,comment):
        self.comment=comment

    # def set_bin(self,bin):
    #     self.bin=bin

