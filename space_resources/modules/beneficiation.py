# -*- coding: utf-8 -*-
"""
Created on Fri Aug 12 16:41:41 2022

Author: Dorian Leger
"""

#transforms regolith input values into the corresponding output values given the beneficiation parameters and saves those in a class object
class Benef_class():

    def __init__(self, B_in_regolith, pre_benef_ilmenite_grade, enrichment_factor = 6, benef_ilmenite_recovery = 0.51):

        self.B_in_regolith = B_in_regolith
        self.pre_benef_ilmenite_grade = pre_benef_ilmenite_grade

        self.enrichment_factor = enrichment_factor
        self.benef_ilmenite_recovery = benef_ilmenite_recovery

        self.post_benef_ilmenite_grade = self.pre_benef_ilmenite_grade * enrichment_factor
        if(self.post_benef_ilmenite_grade >= 0.99):
            self.post_benef_ilmenite_grade = 0.99

        self.B_in_ilmenite = self.B_in_regolith * self.pre_benef_ilmenite_grade
        self.B_out_ilmenite = self.B_in_ilmenite * self.benef_ilmenite_recovery
        
        self.B_out_regolith =  self.B_out_ilmenite / self.post_benef_ilmenite_grade

        #self.B_out_gangue = (self.B_out_ilmenite-self.enrichment_factor*self.B_out_ilmenite *
                             #self.pre_benef_ilmenite_grade)/(self.enrichment_factor * self.pre_benef_ilmenite_grade)

        #self.B_out_regolith = self.B_out_ilmenite + self.B_out_gangue
        self.B_out_gangue = self.B_out_regolith - self.B_out_ilmenite
