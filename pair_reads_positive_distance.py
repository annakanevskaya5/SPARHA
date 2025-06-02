#!pip install pysam

import pandas as pd
import pysam
from bisect import bisect_left

def extract_positions_by_strand(bam_file, strand):
#extracts the read coordinates and their identifiers
    positions = []
    with pysam.alignmentfile(bam_file, "rb") as samfile:
        for read in samfile.fetch():
            if not read.is_unmapped and not read.is_secondary and not read.is_supplementary:  #only primary reads
                if strand == "+" and not read.is_reverse:
                    positions.append((read.query_name, read.reference_start, read.reference_end))
                elif strand == "-" and read.is_reverse:
                    positions.append((read.query_name, read.reference_end, read.reference_start))
    return positions

def find_nearest_pairs_with_binary_search(plus_positions, minus_positions, max_distance):
    pairs = []
    discarded_pairs = []  #list for discarded reads
    used_plus_reads = set()
    used_minus_reads = set()
#reads from the plus and minus strands are used only once
#list of only the start positions of minus strand reads for binary search
    minus_starts = [m_start for _, m_start, _ in minus_positions]

    for p_id, p_start, p_end in plus_positions:
        if p_id in used_plus_reads:
            continue  #skip the plus strand read if it has already been used

        #binary search to determine the index of the range
        left_index = bisect_left(minus_starts, p_start - max_distance)
        right_index = bisect_left(minus_starts, p_start + 1)

        possible_pairs = []  #list of all possible pairs for the current plus strand read

        for idx in range(left_index, right_index):
            m_id, m_start, m_end = minus_positions[idx]
            if m_id in used_minus_reads:
                continue  #skip the minus strand read if it has already been used

            #the start position of the minus strand read can't be greater than the plus strand read
            if m_start > p_start:
                continue

            dist = abs(p_start - m_start)
            if dist <= max_distance:
                possible_pairs.append(((p_id, p_start, p_end), (m_id, m_start, m_end), dist))

        if possible_pairs:
            #the pair with the minimum distance
            best_pair = min(possible_pairs, key=lambda x: x[2])
            pairs.append(best_pair)
            used_plus_reads.add(best_pair[0][0])  #mark the plus strand read as used
            used_minus_reads.add(best_pair[1][0]) #mark the minus strand read as used

        else:
            #if no suitable pair is found, add the read to the discarded list
            discarded_pairs.append((p_id, p_start, p_end))

    return pairs, discarded_pairs

def write_pairs_to_excel(output_file, pairs):
    #create a table for the read pairs
    pair_data = []
    for (p_id, p_start, p_end), (m_id, m_start, m_end), dist in pairs:
        pair_data.append([p_id, f"{p_start}-{p_end}", "+", m_id, f"{m_start}-{m_end}", "-", dist])

    pairs_df = pd.DataFrame(pair_data, columns=["read on + strand", "coordinates + strand", "strand +",
                                                "read on - strand", "coordinates - strand", "strand -", "distance"])

    #write to an excel file
    pairs_df.to_excel(output_file, index=False)

def write_discarded_pairs_to_excel(discarded_file, discarded_pairs):
    discarded_data = []
    for p_id, p_start, p_end in discarded_pairs:
        discarded_data.append([p_id, f"{p_start}-{p_end}", "+"])

    discarded_df = pd.DataFrame(discarded_data, columns=["read on + strand", "coordinates + strand", "strand +"])
    discarded_df.to_excel(discarded_file, index=False)

def calculate_statistics(pairs, max_distance):
    distances = [dist for _, _, dist in pairs]
    total_pairs = len(pairs)
    avg_distance = sum(distances) / total_pairs if total_pairs > 0 else 0
    distance_distribution = {}

    for dist in distances:
        if dist not in distance_distribution:
            distance_distribution[dist] = 0
        distance_distribution[dist] += 1

    #statistics
    stats_data = []
    for i in range(max_distance + 1):
        count = distance_distribution.get(i, 0)
        stats_data.append([i, count, total_pairs, avg_distance])

    stats_df = pd.DataFrame(stats_data, columns=["distance", "count", "total count", "average distance"])
    return stats_df

def write_statistics_to_excel(stats_file, stats_df):
    stats_df.to_excel(stats_file, index=False)

def main():
    bam_file = "file.bam"  #BAM file
    max_distance = NUMBER, int  #maximum distance between reads
    pairs_output_file = "valid_paired.xlsx"  #valid read pairs
    discarded_output_file = "discarded_paired.xlsx"  #discarded read pairs
    stats_output_file = "statistics.xlsx"  #statistics


    # Calculate statistics
    stats_df = calculate_statistics(pairs, max_distance)

if __name__ == "__main__":
    main()
