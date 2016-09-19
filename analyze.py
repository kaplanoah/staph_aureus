from collections import defaultdict
import operator
from math import ceil


# settings

sequence_length = 12
ignore_center_length = 4

min_print_representation = 5
max_print_count = 50


# get genome

with open ('staph_aureus_genome_clean.txt', 'r') as staph_aureus_file:
    staph_aureus_genome = staph_aureus_file.read().upper()


# calculate base frequencies

base_counts = defaultdict(int)

for base in staph_aureus_genome:
    base_counts[base] += 1

base_frequencies = {}

for base, count in base_counts.iteritems():
    base_frequencies[base] = float(count) / len(staph_aureus_genome)


# count sequence occurrences

sequence_occurrence_counts = defaultdict(int)

for start_index in xrange(len(staph_aureus_genome) - sequence_length + 1):

    sequence_end_length = int(ceil((float(sequence_length) - ignore_center_length) / 2))

    sequence = staph_aureus_genome[start_index:start_index + sequence_length]

    sequence_ignoring_center = '%s%s%s' % (
            sequence[:sequence_end_length],
            '-' * ignore_center_length,
            sequence[sequence_end_length + ignore_center_length:]
        )

    sequence_occurrence_counts[sequence_ignoring_center] += 1


# calculate sequence probabilities

sequence_probabilities = {}

for sequence, occurrence_count in sequence_occurrence_counts.iteritems():

    sequence_probability = 1

    for base in sequence:
        sequence_probability *= base_frequencies.get(base, 1)

    sequence_probabilities[sequence] = sequence_probability


# compare actual and expected counts

sequence_representations = {}

for sequence, probability in sequence_probabilities.iteritems():

    expected_occurrence_count = probability * len(staph_aureus_genome)
    actual_occurrence_count   = sequence_occurrence_counts[sequence]

    sequence_representations[sequence] = actual_occurrence_count / expected_occurrence_count


# sort by overrepresented to underrepresented

sorted_representations = sorted(
        sequence_representations.items(),
        key=operator.itemgetter(1),
        reverse=True,
    )


# print overrepresented sequences

for print_count, (sequence, representation) in enumerate(sorted_representations):

    if representation < min_print_representation or print_count == max_print_count:
        break

    expected_occurrence_count = sequence_probabilities[sequence] * len(staph_aureus_genome)
    actual_occurrence_count   = sequence_occurrence_counts[sequence]

    print sequence, '%6.2f' % representation,
    print '%6.2f' % expected_occurrence_count, '%4i' % actual_occurrence_count
