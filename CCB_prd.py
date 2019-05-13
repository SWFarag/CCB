#!/usr/bin/env python
# coding: utf-8

from itertools import combinations,combinations_with_replacement,tee, chain
import pandas as pd
import numpy as np
import random
import argparse
import sys
import os
from os import listdir
from os.path import isfile, join
from IPython.display import display

class CCB:
    def __init__(self, in_df_path=None, out_df_path=None, peptide_length=None, genus_type=None, replacement=None):
        self.in_df_path = in_df_path
        self.out_df_path = out_df_path
        self.peptide_length = peptide_length
        self.genus_type = genus_type
        self.replacement = replacement


    def unziping(self, vip_seqs):
        res_list = []
        for seq in vip_seqs:
            res =[]
            merged = list(chain(*seq))
            res = sorted(set(merged), key=merged.index)
            res_list.append(res)
        return res_list

    def unziping2(self, vip_seqs):
        res_list = []
        for t in vip_seqs:
            temp = [t[0][0], t[0][1]]
            for i in range(1,len(t)):
                a = t[i]
                temp.append(a[1])
            res_list.append(temp)
        return res_list


    def pairwise(self, iterable):
        "s -> (s0,s1), (s1,s2), (s2, s3), ..."
        a, b = tee(iterable)
        next(b, None)
        return zip(a, b)

    def lookForVipSeq(self, result, unique_pairs_edited):
        vip_seqs = []
        for pot_seq in result:
            pairs_ = []
            potential_pairs = self.pairwise(pot_seq)
            for j in potential_pairs:
                pairs_.append(j)
            if(all((i in unique_pairs_edited) for i in pairs_)):
                vip_seqs.append(pairs_)
        return vip_seqs


    def rSubset(self, arr, r):

        # return list of all subsets of length r
        # to deal with duplicate subsets use
        # set(list(combinations(arr, r)))
        return list(combinations(arr, r))

    def rSubset_w_rep(self, arr, r):

        # return list of all subsets of length r
        # to deal with duplicate subsets use
        # set(list(combinations(arr, r)))
        return list(combinations_with_replacement(arr, r))

    def perUniquepairsFormat(self, unique_pairs):
        unique_pairs_edited = []
        for p in unique_pairs:
            temp = p.split("_")
            unique_pairs_edited.append((temp[0].lower(),temp[1].lower()))
        return unique_pairs_edited

    def getVipSeq(self, in_df_path, out_df_path, peptide_length, genus_type, replacement):
        suffix = peptide_length
        all_linkers = pd.read_csv(in_df_path)
        print("all_linkers_shape: ", all_linkers.shape)
        print("all_linkers_pairs: ", len(all_linkers['Pairs']))

        if(genus_type):
            all_linkers_genus =  all_linkers[all_linkers['Genus']==genus_type]
            unique_pairs = list(set(all_linkers_genus['Pairs']))
            print("genus_all_unique_pairs: ", len(unique_pairs))

            a1 = list(set(all_linkers_genus["A1"]))
            a2 = list(set(all_linkers_genus["A2"]))
            unique_mono = list(set(a1+a2))
            unique_mono.sort()
            print("all_unique_mono: ", len(unique_mono))

        else:
            unique_pairs = list(set(all_linkers['Pairs']))
            print("all_unique_pairs: ", len(unique_pairs))


            a1 = list(set(all_linkers["A1"]))
            a2 = list(set(all_linkers["A2"]))

            unique_mono = list(set(a1+a2))
            unique_mono.sort()
            print("all_unique_mono: ", len(unique_mono))

        if(replacement):
            result = self.rSubset_w_rep(unique_mono, peptide_length)
            print("Number of Possible combination with replacement: ", len(result))
        else:
            result = self.rSubset(unique_mono, peptide_length)
            print("Number of Possible combination without replacement: ", len(result))


        unique_pairs_edited = self.perUniquepairsFormat(unique_pairs)
        print("unique_pairs_edited is done", len(unique_pairs_edited))

        vip_seqs = self.lookForVipSeq(result, unique_pairs_edited)
        print("vip_seqs is done", len(vip_seqs))

        res_list = self.unziping2(vip_seqs)
        print("res_list is done", len(res_list))

        self.writeRes(res_list, genus_type, out_df_path, suffix)

    def writeRes(self, res_list, genus_type, out_df_path, suffix):
        fileName = ""
        if(genus_type):
            fileName = out_df_path + "/" + "potential_sequences_" +  str(suffix) + "_" + genus_type +".csv"
        else:
            fileName = out_df_path + "/" + "potential_sequences_" + str(suffix )+ ".csv"
        dict_out = {}
        dict_out["Peptides"] = res_list
        df = pd.DataFrame.from_dict(dict_out)
        df.to_csv(fileName)

    def test_run(self):
        folderName = ""
        if(self.replacement is None): self.replacement = False
        if ((not self.replacement) and (self.genus_type is None) ): folderName = "no_rep_no_genus/"
        elif((not self.replacement) and  (self.genus_type is not None)):
            folderName = "no_rep_genus/"
            if(not self.checkGenus(self.genus_type)):
                raise Exception('genus type does not exist. The genus type was: {}'.format(self.genus_type))
        elif( (self.replacement) and (self.genus_type is None)): folderName = "rep_no_genus/"
        else:
            folderName = "rep_genus/"
            if(not self.checkGenus(self.genus_type)):
                raise Exception('genus type does not exist. The genus type was: {}'.format(self.genus_type))


        # define the name of the directory to be created
        new_output_path = self.out_df_path + folderName
        try:
            os.mkdir(new_output_path)
        except OSError:
            print ("Creation of the directory %s failed" % new_output_path)
        else:
            print ("Successfully created the directory %s " % new_output_path)

        self.getVipSeq(self.in_df_path, new_output_path, self.peptide_length, self.genus_type, self.replacement)

    def checkGenus(self, genus):
        isGenus = False
        print("ich bin da")
        all_linkers = pd.read_csv(self.in_df_path)
        try:
            genus_df = all_linkers[all_linkers['Genus']==self.genus_type]
            isGenus = True
        except e:
            print("eshta" + e)
        return isGenus


if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-in','--in_df_path',
                        help='Path of input file path.',  required=True, default=None)
    parser.add_argument('-o', '--out_df_path',
                        help='Path to new .csv file for saving the potential peptide dataframe.', required=True, default=None)
    parser.add_argument('-l', '--peptide_length', type=int,
                        help='desired peptide length.', required=True, default=None)
    parser.add_argument('-g', '--genus_type',
                        help='Genus type of the Bacteria.', default=None)
    parser.add_argument('-r', '--replacement',
                        help='Boolean for combinatorial with replacement.', default=False)

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()
    ccb = CCB(in_df_path=args.in_df_path,
                           out_df_path=args.out_df_path,
                           peptide_length=args.peptide_length,
                           genus_type=args.genus_type,
                           replacement=args.replacement)
    ccb.test_run()
