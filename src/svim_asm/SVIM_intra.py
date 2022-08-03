from __future__ import print_function

import sys
from bisect import bisect_left

from svim_asm.SVCandidate import CandidateDeletion, CandidateInsertion


def analyze_cigar_indel(tuples, min_length):
    """Parses CIGAR tuples (op, len) and returns Indels with a length > minLength"""
    pos_ref = 0
    pos_read = 0
    indels = []
    for operation, length in tuples:
        if operation == 0:                     # alignment match
            pos_ref += length
            pos_read += length
        elif operation == 1:                   # insertion
            if length >= min_length:
                indels.append((pos_ref, pos_read, length, "INS"))
            pos_read += length
        elif operation == 2:                   # deletion
            if length >= min_length:
                indels.append((pos_ref, pos_read, length, "DEL"))
            pos_ref += length
        elif operation == 4:                   # soft clip
            pos_read += length
        elif operation == 7 or operation == 8:        # match or mismatch
            pos_ref += length
            pos_read += length
    return indels


def combine_indels(vntr_indels):
    if len(vntr_indels) == 1:
        return vntr_indels[0]

    combined_size = 0
    for indel in vntr_indels:
        if indel[3] == "INS":
            combined_size += indel[2]
        else:
            combined_size -= indel[2]
    if combined_size >= 0:
        return vntr_indels[0][0], vntr_indels[0][1], combined_size, "INS"
    else:
        return vntr_indels[0][0], vntr_indels[0][1], -combined_size, "DEL"


def group_vntr_indels(aln_indels, tandem_annotations, min_length, ref_start):
    vntr_starts = [x[0] for x in tandem_annotations]
    if not vntr_starts:
        return aln_indels

    new_indels = []
    prev_vntr_id = 0
    vntr_cluster = []
    for ref_pos, qry_pos, indel_size, indel_type in aln_indels:
        vntr_id = 0
        idx = bisect_left(vntr_starts, ref_pos + ref_start)
        if idx > 0 and ref_pos + ref_start < tandem_annotations[idx - 1][1]:
            vntr_id = idx - 1

        if vntr_id == 0:
            new_indels.append((ref_pos, qry_pos, indel_size, indel_type))
            continue

        if vntr_id != prev_vntr_id:
            if vntr_cluster:
                new_indels.append(combine_indels(vntr_cluster))
            prev_vntr_id = vntr_id
            vntr_cluster = []
        vntr_cluster.append((ref_pos, qry_pos, indel_size, indel_type))

    if vntr_cluster:
        new_indels.append(combine_indels(vntr_cluster))

    new_indels = [x for x in new_indels if x[2] >= min_length]

    return new_indels


def analyze_alignment_indel(alignment, bam, query_name, options, tandem_annotations):
    sv_candidates = []
    ref_chr = bam.getrname(alignment.reference_id)
    ref_start = alignment.reference_start

    cigar_min_size = options.min_sv_size if not tandem_annotations else 10
    indels = analyze_cigar_indel(alignment.cigartuples, cigar_min_size)
    if tandem_annotations:
        indels = group_vntr_indels(indels, tandem_annotations[ref_chr], options.min_sv_size, ref_start)

    for pos_ref, pos_read, length, typ in indels:
        if typ == "DEL":
            sv_candidates.append(CandidateDeletion(ref_chr, ref_start + pos_ref, ref_start + pos_ref + length, [query_name], bam))
        elif typ == "INS":
            insertion_seq = alignment.query_sequence[pos_read:pos_read+length]
            sv_candidates.append(CandidateInsertion(ref_chr, ref_start + pos_ref, ref_start + pos_ref + length, [query_name], insertion_seq, bam))
    return sv_candidates


