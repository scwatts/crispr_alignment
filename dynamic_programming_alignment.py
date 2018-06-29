#!/usr/bin/env python3
from itertools import product


import skbio.alignment._pairwise
import skbio.sequence
import skbio.util
import numpy


class SpacerSequence(skbio.sequence.GrammaredSequence):

    __validation_mask = None
    _number_of_extended_ascii_codes = 256

    @skbio.util.classproperty
    def definite_chars(cls):
        return set(i for i in numpy.array(range(1,256), dtype=numpy.uint8))

    @skbio.util.classproperty
    def degenerate_map(cls):
        return dict()

    @skbio.util.classproperty
    def default_gap_char(cls):
        return numpy.uint8(0)

    @skbio.util.classproperty
    def gap_chars(cls):
        return {numpy.uint8(0)}

    @skbio.util.classproperty
    def _validation_mask(cls):
        # TODO These masks could be defined (as literals) on each concrete
        # object. For now, memoize!
        if cls.__validation_mask is None:
            cls.__validation_mask = numpy.invert(numpy.bincount(
                numpy.array(list(cls.alphabet), dtype=numpy.uint8),
                minlength=cls._number_of_extended_ascii_codes).astype(bool))
        return cls.__validation_mask


def compute_substitution_score(aln1_chars, aln2_chars, substitution_matrix,
                                gap_substitution_score, gap_chars):
    substitution_score = 0
    for aln1_char, aln2_char in product(aln1_chars, aln2_chars):
        if aln1_char in gap_chars or aln2_char in gap_chars:
                substitution_score += gap_substitution_score
        elif aln1_char == aln2_char:
                substitution_score += 1
    substitution_score /= (len(aln1_chars) * len(aln2_chars))
    return substitution_score


def traceback(traceback_matrix, score_matrix, aln1, aln2, start_row,
               start_col):
    # cache some values for simpler reference
    aend = skbio.alignment._pairwise._traceback_encoding['alignment-end']
    match = skbio.alignment._pairwise._traceback_encoding['match']
    vgap = skbio.alignment._pairwise._traceback_encoding['vertical-gap']
    hgap = skbio.alignment._pairwise._traceback_encoding['horizontal-gap']
    gap_character = aln1.dtype.default_gap_char

    # initialize the result alignments
    aln1_sequence_count = aln1.shape.sequence
    aligned_seqs1 = [[] for e in range(aln1_sequence_count)]

    aln2_sequence_count = aln2.shape.sequence
    aligned_seqs2 = [[] for e in range(aln2_sequence_count)]

    current_row = start_row
    current_col = start_col

    best_score = score_matrix[current_row, current_col]
    current_value = None

    while current_value != aend:
        current_value = traceback_matrix[current_row, current_col]

        if current_value == match:
            for aligned_seq, input_seq in zip(aligned_seqs1, aln1):
                aligned_seq.append(str(input_seq[current_col-1]))
            for aligned_seq, input_seq in zip(aligned_seqs2, aln2):
                aligned_seq.append(str(input_seq[current_row-1]))
            current_row -= 1
            current_col -= 1
        elif current_value == vgap:
            for aligned_seq in aligned_seqs1:
                aligned_seq.append(gap_character)
            for aligned_seq, input_seq in zip(aligned_seqs2, aln2):
                aligned_seq.append(str(input_seq[current_row-1]))
            current_row -= 1
        elif current_value == hgap:
            for aligned_seq, input_seq in zip(aligned_seqs1, aln1):
                aligned_seq.append(str(input_seq[current_col-1]))
            for aligned_seq in aligned_seqs2:
                aligned_seq.append(gap_character)
            current_col -= 1
        elif current_value == aend:
            continue
        else:
            raise ValueError(
                "Invalid value in traceback matrix: %s" % current_value)

    #for i, (aligned_seq, original) in enumerate(zip(aligned_seqs1, aln1)):
    #    aligned_seq = ''.join(str(i) for i in aligned_seq)[::-1]
    #    constructor = aln1.dtype
    #    metadata = None
    #    if original.has_metadata():
    #        metadata = original.metadata
    #    aligned_seqs1[i] = constructor(aligned_seq, metadata=metadata,
    #                                   validate=False)

    #for i, (aligned_seq, original) in enumerate(zip(aligned_seqs2, aln2)):
    #    aligned_seq = ''.join(str(i) for i in aligned_seq)[::-1]
    #    constructor = aln2.dtype
    #    metadata = None
    #    if original.has_metadata():
    #        metadata = original.metadata
    #    aligned_seqs2[i] = constructor(aligned_seq, metadata=metadata,
    #                                   validate=False)

    return aligned_seqs1, aligned_seqs2, best_score, current_col, current_row


def global_pairwise_align(seq1, seq2, gap_open_penalty, gap_extend_penalty,
                        substitution_matrix, penalize_terminal_gaps=False):
    for seq in seq1, seq2:
        # We don't need to check the case where `seq` is a `TabularMSA` with a
        # dtype that isn't a subclass of `GrammaredSequence`, this is
        # guaranteed by `TabularMSA`.
        if not isinstance(seq, (skbio.sequence.GrammaredSequence, skbio.alignment.TabularMSA)):
            raise TypeError(
                "`seq1` and `seq2` must be GrammaredSequence subclasses or "
                "TabularMSA, not type %r" % type(seq).__name__)

    seq1 = skbio.alignment._pairwise._coerce_alignment_input_type(seq1)
    seq2 = skbio.alignment._pairwise._coerce_alignment_input_type(seq2)

    if seq1.dtype is not seq2.dtype:
        raise TypeError(
            "`seq1` and `seq2` must have the same dtype: %r != %r"
            % (seq1.dtype.__name__, seq2.dtype.__name__))

    if penalize_terminal_gaps:
        init_matrices_f = skbio.alignment._pairwise._init_matrices_nw
    else:
        init_matrices_f = skbio.alignment._pairwise._init_matrices_nw_no_terminal_gap_penalty

    score_matrix, traceback_matrix = \
        skbio.alignment._pairwise._compute_score_and_traceback_matrices(
            seq1, seq2, gap_open_penalty, gap_extend_penalty,
            substitution_matrix, new_alignment_score=-numpy.inf,
            init_matrices_f=init_matrices_f,
            penalize_terminal_gaps=penalize_terminal_gaps)

    end_row_position = traceback_matrix.shape[0] - 1
    end_col_position = traceback_matrix.shape[1] - 1

    aligned1, aligned2, score, seq1_start_position, seq2_start_position = \
        skbio.alignment._pairwise._traceback(traceback_matrix, score_matrix, seq1, seq2,
                   end_row_position, end_col_position)
    start_end_positions = [(seq1_start_position, end_col_position-1),
                           (seq2_start_position, end_row_position-1)]

    result = [convert_alignments(aligned1[0]),
              convert_alignments(aligned2[0])]
    return result, score, start_end_positions


def convert_alignments(a):
    n = list()
    for d in a:
        if d == 0:
            n.append('-')
            continue
        try:
            n.append(ord(d))
        except TypeError:
            n.append(int(d))
    return n[::-1]


skbio.alignment._pairwise._compute_substitution_score = compute_substitution_score
skbio.alignment._pairwise._traceback = traceback
skbio.alignment._pairwise.global_pairwise_align = global_pairwise_align


def main():
    #seq1_elements = [1, 2, 3, 4, 5, 10, 13, 14, 15]
    #seq2_elements = [1, 2, 3, 4, 5, 6, 7, 8, 9, 13, 14, 15, 127]

    seq1_elements = [1, 2, 3, 4, 6, 6, 8, 9]
    seq2_elements = [1, 2, 3, 4, 6, 5, 7, 8, 9]

    seq1 = SpacerSequence(numpy.array(seq1_elements, dtype=numpy.uint8))
    seq2 = SpacerSequence(numpy.array(seq2_elements, dtype=numpy.uint8))

    alphabet = [numpy.uint8(i) for i in SpacerSequence.alphabet]
    identity_mat = {c1: {c2: 0 for c2 in alphabet} for c1 in alphabet}
    for c in SpacerSequence.alphabet:
        identity_mat[c][c] = 1

    alignment, score, pos = skbio.alignment._pairwise.global_pairwise_align(seq1, seq2, 0.5, 0, identity_mat, False)
    seq1_aln, seq2_aln = alignment
    print(*seq1_aln, sep='\t')
    print(*seq2_aln, sep='\t')


if __name__ == '__main__':
    main()
