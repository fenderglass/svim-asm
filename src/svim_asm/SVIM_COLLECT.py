import logging
import pysam
from copy import deepcopy

from collections import defaultdict

from svim_asm.SVIM_intra import analyze_alignment_indel
from svim_asm.SVIM_inter import analyze_read_segments


def retrieve_other_alignments(main_alignment, bam):
    """Reconstruct other alignments of the same read for a given alignment from the SA tag"""
    #reconstructing other alignments from SA tag does not work if sequence of main_alignment is hard-clipped
    if main_alignment.get_cigar_stats()[0][5] > 0:
        return []
    try:
        sa_tag = main_alignment.get_tag("SA").split(";")
    except KeyError:
        return []
    other_alignments = []
    # For each other alignment encoded in the SA tag
    for element in sa_tag:
        # Read information from the tag
        fields = element.split(",")
        if len(fields) != 6:
            continue
        rname = fields[0]
        pos = int(fields[1])
        strand = fields[2]
        # CIGAR string encoded in SA tag is shortened
        cigar = fields[3]
        mapq = int(fields[4])
        nm = int(fields[5])

        # Generate an aligned segment from the information
        a = pysam.AlignedSegment()
        a.query_name = main_alignment.query_name
        a.query_sequence = ''
        if strand == "+":
            a.flag = 2048
        else:
            a.flag = 2064
        a.reference_id = bam.get_tid(rname)
        a.reference_start = pos - 1
        try:
            a.mapping_quality = mapq
        except OverflowError:
            a.mapping_quality = 0
        try:
            a.cigarstring = cigar
        except OverflowError:
            logging.error("OverflowError while retrieving supplementary CIGAR string. Read name: {0}, Position: {1}, CIGAR: {2}".format(rname, pos, cigar))
            continue
        a.next_reference_id = -1
        a.next_reference_start = -1
        a.template_length = 0
        a.query_qualities = ''
        a.set_tags([("NM", nm, "i")])

        other_alignments.append(a)
    return other_alignments


def filter_contained_alignments(alignments):
    CONTAINMENT_RATE = 0.9
    aln_intervals = defaultdict(list)
    for aln in alignments:
        ref_id, ref_start, ref_end = aln.reference_name, aln.reference_start, aln.reference_end
        if ref_id:
            aln_intervals[ref_id].append((ref_start, ref_end))

    for chr_id in aln_intervals:
        aln_intervals[chr_id].sort(key=lambda p: p[1] - p[0], reverse=True)

    filtered_alignments = []
    for aln in alignments:
        if not aln.reference_name:
            continue
        ref_id, ref_start, ref_end = aln.reference_name, aln.reference_start, aln.reference_end

        contained = False
        for other_aln in aln_intervals[ref_id]:
            overlap = min(other_aln[1], ref_end) - max(other_aln[0], ref_start)
            if overlap > CONTAINMENT_RATE * (ref_end - ref_start) and (other_aln[1] - other_aln[0]) > (ref_end - ref_start):
                contained = True
                break

        if not contained:
            filtered_alignments.append(aln)

    return filtered_alignments


def parse_tandem_bed(filename):
    tandem_by_chr = defaultdict(list)
    for line in open(filename, "r"):
        if line.startswith("#"):
            continue
        fields = line.strip().split()
        tandem_by_chr[fields[0]].append((int(fields[1]), int(fields[2])))
    return tandem_by_chr


def analyze_alignment_file_coordsorted(bam, options, bed_file):
    tandem_annotations = None
    if options.tandem:
        tandem_annotations = parse_tandem_bed(options.tandem)

    all_alignments = []
    for aln in bam.fetch():
        if aln.is_unmapped or aln.is_secondary or aln.mapping_quality < options.min_mapq:
            continue
        all_alignments.append(aln)

    if options.filter_contained:
        all_alignments = filter_contained_alignments(all_alignments)

    with open(bed_file, "w") as f:
        for aln in all_alignments:
            f.write("{0}\t{1}\t{2}\t{3}\n".format(aln.reference_name, aln.reference_start,
                                                  aln.reference_end, aln.query_name))

    supplementary_aln_by_read = defaultdict(list)
    for aln in all_alignments:
        if aln.is_supplementary:
            supplementary_aln_by_read[aln.query_name].append(aln)

    sv_candidates = []
    for current_alignment in all_alignments:
        sv_candidates.extend(analyze_alignment_indel(current_alignment, bam, current_alignment.query_name,
                                                     options, tandem_annotations))
        if not current_alignment.is_supplementary:
            good_suppl_alns = supplementary_aln_by_read[current_alignment.query_name]
            #_supplementary_alignments = retrieve_other_alignments(current_alignment, bam)
            #_good_suppl_alns = [aln for aln in _supplementary_alignments if not aln.is_unmapped and aln.mapping_quality >= options.min_mapq]

            sv_candidates.extend(analyze_read_segments(current_alignment, good_suppl_alns, bam, options))

    return sv_candidates
