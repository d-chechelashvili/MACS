import operator
import sys

KMERLEN = 6


def read_file(file_name):
    data = []
    for line in open(file_name, 'r'):
        data.append(line[0:-1] if line[-1] == "#" else line)
    return data


def A(kmer_counts, kmer_starts, kmer_conservation_counts):
    sorted_kmer_counts = sorted(kmer_counts.items(), key=operator.itemgetter(1), reverse=True)
    top_50_frequent = []
    for kmer, _ in sorted_kmer_counts:
        top_50_frequent.append(kmer)
        if len(top_50_frequent) == 50:
            break

    kmer_conservations = {}
    for kmer in kmer_conservation_counts:
        conservation_score = 0
        if kmer_conservation_counts[kmer] > 0:
            conservation_score = float(float(kmer_counts[kmer]) / float(kmer_conservation_counts[kmer]))
        kmer_conservations[kmer] = conservation_score
    sorted_kmer_conservations = sorted(kmer_conservations.items(), key=operator.itemgetter(1), reverse=True)
    top_50_conserved = []
    for kmer, _ in sorted_kmer_conservations:
        top_50_conserved.append(kmer)
        if len(top_50_conserved) == 50:
            break
    return top_50_frequent, top_50_conserved


def read_yeast_file(yeast_file_name):
    data = []
    for line in open(yeast_file_name, 'r'):
        if line[-1] == "\n":
            line = line[0:-1]
        motif_name_seq = line.split(" ")
        if len(motif_name_seq) < 2:
            continue
        data.append(motif_name_seq[1])
    return data[1:]


def B(all_kmers, yeast_file_name):
    yeast_sequences = read_yeast_file(yeast_file_name)
    found_kmers = {}
    for kmer in all_kmers:
        for yeast_sequence in yeast_sequences:
            if kmer in yeast_sequence:
                if kmer in found_kmers:
                    if yeast_sequence not in found_kmers[kmer]:
                        found_kmers[kmer].append(yeast_sequence)
                else:
                    found_kmers[kmer] = [yeast_sequence]
    return found_kmers


def main():
    intergenic_file_name = sys.argv[1]
    annotation_file_name = sys.argv[2]
    yeast_file_name = sys.argv[3]

    intergenic_regions = read_file(intergenic_file_name)[0]
    conservation_annotations = read_file(annotation_file_name)[0]

    kmer_counts = {}
    kmer_starts = {}
    kmer_conservation_counts = {}
    for i in xrange(len(intergenic_regions) - KMERLEN + 1):
        kmer = intergenic_regions[i:i + KMERLEN]
        conserved_nucleotide_num = conservation_annotations[i:i + KMERLEN].count("*")
        if kmer not in kmer_counts:
            kmer_counts[kmer] = 1
            kmer_starts[kmer] = [i]
            kmer_conservation_counts[kmer] = conserved_nucleotide_num
        else:
            kmer_counts[kmer] += 1
            kmer_starts[kmer].append(i)
            kmer_conservation_counts[kmer] += conserved_nucleotide_num

    # A
    top_50_frequent, top_50_conserved = A(kmer_counts, kmer_starts, kmer_conservation_counts)

    for i in xrange(50):
        print top_50_frequent[i] + "         |         " + top_50_conserved[i]

    # B
    found_kmers = B(kmer_counts.keys(), yeast_file_name)
    for kmer in found_kmers:
        print kmer + " has been found in:" + str(found_kmers[kmer])


if __name__ == "__main__":
    main()
