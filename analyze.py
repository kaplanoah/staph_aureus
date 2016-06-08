from collections import defaultdict
import operator


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


# count sequence occurences

sequence_occurence_counts = defaultdict(int)

sequence_length = 8

for start_index in xrange(len(staph_aureus_genome) - sequence_length + 1):

    sequence = staph_aureus_genome[start_index:start_index + sequence_length]

    sequence_occurence_counts[sequence] += 1


# calculate sequence probabilities

sequence_probabilities = {}

for sequence, occurence_count in sequence_occurence_counts.iteritems():

    sequence_probability = 1

    for base in sequence:
        sequence_probability *= base_frequencies[base]

    sequence_probabilities[sequence] = sequence_probability


# compare actual and expected counts

sequence_representations = {}

for sequence, probability in sequence_probabilities.iteritems():

    expected_occurence_count = probability * len(staph_aureus_genome)
    actual_occurence_count   = sequence_occurence_counts[sequence]

    sequence_representations[sequence] = actual_occurence_count / expected_occurence_count


# sort by overrepresented to underrepresented

sorted_representations = sorted(
        sequence_representations.items(),
        key=operator.itemgetter(1),
        reverse=True,
    )


# print overrepresented sequences

for sequence, representation in sorted_representations:

    if representation < 10:
        break

    print sequence, '  ',
    print round(representation, 2)
