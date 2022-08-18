# -*- coding: utf-8 -*-
"""
Created on Fri Aug 12 16:41:41 2022

"""
  

benef_gangue_recovery = 0.045
enrichment_factor = 6
benef_ilmenite_recovery= 0.42


head_grade_min = 0.01
head_grade_max = 0.15
step = 0.01

def benef(head_grade):

    regolith_in_pre_benef = 1000   # grams
    post_benef_grade = head_grade*enrichment_factor   # %
    gangue_in_pre_benef = regolith_in_pre_benef*(1-head_grade) # mass 
    ilmenite_in_pre_benef = regolith_in_pre_benef * head_grade # mass 

    #post_benef_ilmenite = ilmenite_in_pre_benef * benef_ilmenite_recovery
    post_benef_gangue = gangue_in_pre_benef * benef_gangue_recovery
    post_benef_ilmenite = ilmenite_in_pre_benef * benef_ilmenite_recovery
    #regolith_recovered =  post_benef_gangue+post_benef_ilmenite
    #post_benef_ilmenite_grade = post_benef_grade/(post_benef_gangue+post_benef_ilmenite)


    #print(post_benef_grade)
    print("post benef ilmenite g", post_benef_ilmenite)
    print("post benef gangue g", post_benef_gangue)
    print("post benef ilm %", round(post_benef_grade*100), "% \n")


for i in range(1,16):
    head_grade = i*step
    benef(head_grade)
    
