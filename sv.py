#!/usr/bin/env python

import sys
from utilities import calculate_overlap, calculate_stddev, get_most_common_item

SV_TYPES = ['UNK', 'BND', 'DEL', 'INS', 'INV', 'NOT_SV']

def make_info_dictionary(info):
    spl_info = info.split(";")
    info_dict = {}

    for key_val in spl_info:
        if "=" in key_val:
            key, val = key_val.split("=")
            info_dict[key] = val
        else:
            info_dict[key_val] = ""

    return info_dict


class SV(object):
    """
    Constructs a structural variant (SV) candidate
    """
    def __init__(self,
                 vcf_record,
                 check_type = True,
                 join_mode = False,
                 output_ids = False,
                 ignore_bnd = False,
                 ignore_inv = False):
        self.is_outputting_ids = output_ids
        self.join_mode = join_mode
        self.is_sv = True

        spl_line = vcf_record.rstrip("\n").split("\t")[:8]
        self.check_type = check_type
        self.chromosome = spl_line[0]
        self.begin = int(spl_line[1])
        self.end = self.begin

        if self.is_outputting_ids:
            self.ids = [spl_line[2]]

        # If no INFO, add SVTYPE=UNK
        if spl_line[7] == ".":
            spl_line[7] = ""

        info_dict = make_info_dictionary(spl_line[7])

        if "SVTYPE" not in info_dict:
            ref = spl_line[3]
            alt = spl_line[4]

            # Only biallelic
            alt_spl = [x for x in alt.split(",") if x != "*"]

            if len(alt_spl) == 1:
                if len(ref) >= len(alt) + 50:
                    spl_line[4] = alt_spl[0]

                    if len(spl_line[7]) > 0:
                        spl_line[7] += ";"

                    spl_line[7] = "%sSVTYPE=DEL;SVLEN=%d" % (spl_line[7], len(alt) - len(ref))
                    info_dict["SVTYPE"] = "DEL"
                    info_dict["SVSIZE"] = "%d" % (len(ref) - len(alt))
                    info_dict["SVLEN"] = "%d" % (len(alt) - len(ref))
                elif len(alt) >= len(ref) + 50:
                    spl_line[4] = alt_spl[0]

                    if len(spl_line[7]) > 0:
                        spl_line[7] += ";"

                    spl_line[7] = "%sSVTYPE=INS;SVLEN=%d" % (spl_line[7], len(alt) - len(ref))
                    info_dict["SVTYPE"] = "INS"
                    info_dict["SVSIZE"] = "%d" % (len(alt) - len(ref))
                    info_dict["SVLEN"] = "%d" % (len(alt) - len(ref))
                else:
                  self.is_sv = False
                  return None
            else:
              self.is_sv = False
              return None

            if "SVTYPE" not in info_dict:
                info_dict["SVTYPE"] = "NOT_SV"
                self.is_sv = False
                return None
        else:
            if (ignore_bnd and (info_dict["SVTYPE"] == "BND" or info_dict["SVTYPE"] == "TRA")) or \
               (ignore_inv and info_dict["SVTYPE"] == "INV"):
                self.is_sv = False
                return None

            # Join related SV types
            if info_dict["SVTYPE"] == "DEL_ALU" or info_dict["SVTYPE"] == "DEL_LINE1":
                info_dict["SVTYPE"] = "DEL"
            elif info_dict["SVTYPE"] == "ALU" or info_dict["SVTYPE"] == "LINE1" or info_dict["SVTYPE"] == "SVA" or \
                 info_dict["SVTYPE"] == "DUP" or info_dict["SVTYPE"] == "CNV" or info_dict["SVTYPE"] == "INVDUP" or \
                 info_dict["SVTYPE"] == "INV":
                info_dict["SVTYPE"] = "INS"
            elif info_dict["SVTYPE"] == "TRA":
                info_dict["SVTYPE"] = "BND"



        # Remove old values
        self.old_num_merged_svs = -1

        if not join_mode and "NUM_MERGED_SVS" in info_dict:
          #print(int(info_dict["NUM_MERGED_SVS"]))
          self.old_num_merged_svs = int(info_dict["NUM_MERGED_SVS"])
          spl_line[7] = spl_line[7].replace("NUM_MERGED_SVS=%s;" % info_dict["NUM_MERGED_SVS"], "")

        if "STDDEV_POS" in info_dict:
          spl_line[7] = spl_line[7].replace("STDDEV_POS=%s;" % info_dict["STDDEV_POS"], "")

        if "SVTYPE" in info_dict and info_dict["SVTYPE"] == "DEL":
            if "END" in info_dict:
                self.end = int(info_dict["END"])
            else:
                if "SVSIZE" in info_dict:
                    self.end = self.begin + abs(int(info_dict["SVSIZE"]))
                elif "SVLEN" in info_dict:
                    self.end = self.begin + abs(int(info_dict["SVLEN"]))

            if self.end - 50 < self.begin:
                self.is_sv = False
                return None
        elif "SVTYPE" in info_dict and info_dict["SVTYPE"] == "INS":
            svlen = int(info_dict["SVLEN"]) if "SVLEN" in info_dict else -1

            if svlen == -1:
                svlen = int(info_dict["SVSIZE"]) if "SVSIZE" in info_dict else -1

            if svlen < 50 and "SVINSSEQ" not in info_dict and "LEFT_SVINSSEQ" not in info_dict and "RIGHT_SVINSSEQ" not in info_dict:
                self.is_sv = False
                return None

        if check_type:
            if info_dict["SVTYPE"] in SV_TYPES:
                self.type = SV_TYPES.index(info_dict["SVTYPE"])
            else:
                self.type = 0  # Unknown type
                assert False
        else:
            self.type = 0

        # Get all begins and ends
        self.min_begin = self.begin
        self.max_begin = self.begin
        self.begins = [self.begin]
        self.ends = [self.end]
        self.infos = [spl_line[7]]
        self.unique_begins_and_ends = set()
        self.unique_begins_and_ends.add((int(self.begin), int(self.end)))
        self.refs = [spl_line[3]]
        self.alts = [spl_line[4]]


    """
    Converts the SV to a string
    """
    def __str__(self):
        info = self.infos[0]
        if len(info) > 0:
            info += ";"

        num_text = "JOINED" if self.join_mode else "MERGED"

        num_svs = len(self.begins)

        if self.old_num_merged_svs > 0:
            num_svs += self.old_num_merged_svs - 1

        if self.is_outputting_ids:
            info += "MERGED_IDS=%s;NUM_%s_SVS=%d;STDDEV_POS=%.2f,%.2f" % (",".join(self.ids),
                                                                          num_text,
                                                                          num_svs,
                                                                          calculate_stddev(self.begins),
                                                                          calculate_stddev(self.ends)
                                                                          )
        else:
            info += "NUM_%s_SVS=%d;STDDEV_POS=%.2f,%.2f" % (num_text,
                                                            num_svs,
                                                            calculate_stddev(self.begins),
                                                            calculate_stddev(self.ends)
                                                            )

        return "%s\t%d\t%s\t%s\t%s\t%d\t%s\t%s\n" % (
            self.chromosome,
            self.begin,
            ".",  # variant ID
            self.refs[0],
            self.alts[0],
            0,  # quality
            ".",  # filter
            info
            )


    """
    Finalize SV before printing
    """
    def finalize(self):
        # Get the most common begin and end combo
        begins_and_ends = list(zip(self.begins, self.ends))
        self.begin, self.end = get_most_common_item(begins_and_ends)

        for i, (begin, end) in enumerate(begins_and_ends):
            if begin == self.begin and end == self.end:
                if i > 0:
                    self.infos[0] = self.infos[i]
                    self.refs[0] = self.refs[i]
                    self.alts[0] = self.alts[i]
                return

    """
    Defines how to represent the SV (in printing)
    """
    def __repr__(self):
        return self.__str__()

    """
    Check if some SV should merge with this one. Merging should happen when both SVs are the same type with and overlap
    highly with each other.

    :parameter otherSV The other SV to check.
    :returns True if the two SVs should merge.
    """
    def should_merge(self, other_sv, max_sv_distance, max_size_difference):
        assert isinstance(other_sv, SV)
        assert len(self.begins) == len(self.ends)

        if other_sv.type != self.type:
            return False  # SVs of different types or chromosome should not merge

        # For sanity
        if abs(other_sv.max_begin - self.min_begin) > 10000 or\
                abs(other_sv.min_begin - self.max_begin) > 10000:
            return False

        for begin, end in self.unique_begins_and_ends:
            # For other SVs, calculate their overlap
            if abs(begin - other_sv.begin) <= max_sv_distance and\
                    abs(end - other_sv.end) <= max_sv_distance and\
                    abs((other_sv.end - other_sv.begin) - (end - begin)) <= max_size_difference:
                return True  # Overlap is enough to merge

        # We could not find any interval in this SV that did not enough, we should therefore merge these two SVs
        return False

    """
    Merges two SVs. Assumes that they should be merged, which can be checked with the 'should merge' function

    :parameter otherSV The other SV to merge with.
    :returns Nothing
    """
    def merge(self, other_sv):
        assert isinstance(other_sv, SV)
        assert self.type == other_sv.type

        self.min_begin = min(self.min_begin, other_sv.min_begin)
        self.max_begin = max(self.max_begin, other_sv.max_begin)
        self.begins += other_sv.begins
        self.ends += other_sv.ends
        self.infos += other_sv.infos
        self.refs += other_sv.refs
        self.alts += other_sv.alts
        self.unique_begins_and_ends = self.unique_begins_and_ends.union(other_sv.unique_begins_and_ends)

        if other_sv.old_num_merged_svs > 0:
            self.old_num_merged_svs += other_sv.old_num_merged_svs - 1

        if self.is_outputting_ids:
            self.ids += other_sv.ids
