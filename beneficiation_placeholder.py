# -*- coding: utf-8 -*-
"""
Created on Fri Aug 12 16:41:41 2022


DL Test comment 1
"""


class Benef_class():
    
    def __init__(self,B_in_regolith,pre_benef_ilmenite_grade):
    
        self.B_in_regolith = B_in_regolith
        self.pre_benef_ilmenite_grade=pre_benef_ilmenite_grade
        
        self.enrichment_factor = 6
        self.benef_ilmenite_recovery= 0.51
    
        self.B_in_ilmenite = self.B_in_regolith * self.pre_benef_ilmenite_grade
        self.B_out_ilmenite = self.B_in_ilmenite * self.benef_ilmenite_recovery     
        self.B_out_gangue = (self.B_out_ilmenite-self.enrichment_factor*self.B_out_ilmenite*self.pre_benef_ilmenite_grade)/(self.enrichment_factor *self.pre_benef_ilmenite_grade)
        self.B_out_regolith = self.B_out_ilmenite + self.B_out_gangue

#output ilmenite 
#output gangue
#output regolith




