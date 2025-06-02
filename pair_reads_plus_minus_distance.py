#!pip install pysam

import pandas as pd
import pysam
from bisect import bisect_left

def extract_positions_by_strand(bam_file, strand):
#only primary reads
    positions = []
    with pysam.alignmentfile(bam_file, "rb") as samfile:
        for read in samfile.fetch():
            if not read.is_unmapped and not read.is_secondary and not read.is_supplementary:
                if strand == "+" and not read.is_reverse:
                    positions.append((read.query_name, read.reference_start, read.reference_end))
                elif strand == "-" and read.is_reverse:
                    positions.append((read.query_name, read.reference_end, read.reference_start))
    return positions

def find_all_pairs_within_distance(plus_positions, minus_positions, max_distance):
    pairs = []
    minus_starts = [m_start for _, m_start, _ in minus_positions]

    for p_id, p_start, p_end in plus_positions:
        #binary search
        left_index = bisect_left(minus_starts, p_start - max_distance)
        right_index = bisect_left(minus_starts, p_start + max_distance)

        for idx in range(left_index, right_index):
            m_id, m_start, m_end = minus_positions[idx]
            distance = p_start - m_start


            if distance < 0:
                distance = -abs(distance)
            else:
                distance = abs(distance)

            if abs(distance) <= max_distance:
                pairs.append(((p_id, p_start, p_end), (m_id, m_start, m_end), distance))
    return pairs

def write_pairs_to_excel(output_file, pairs):
    pair_data = []
    for (p_id, p_start, p_end), (m_id, m_start, m_end), dist in pairs:
        pair_data.append([p_id, f"{p_start}-{p_end}", "+", m_id, f"{m_start}-{m_end}", "-", dist])

    pairs_df = pd.DataFrame(pair_data, columns=["+ read ID", "+ coordinates", "+ strand",
                                                "- read ID", "- coordinates", "- strand", "distance"])
    pairs_df.to_excel(output_file, index=False)

def write_discarded_pairs_to_excel(discarded_file, discarded_pairs):
    discarded_data = []
    for p_id, p_start, p_end in discarded_pairs:
        discarded_data.append([p_id, f"{p_start}-{p_end}", "+"])

    discarded_df = pd.DataFrame(discarded_data, columns=["+ read ID", "+ coordinates", "+ strand"])
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

    stats_data = []
    for i in range(-max_distance, max_distance + 1):
        count = distance_distribution.get(i, 0)
        stats_data.append([i, count, total_pairs, avg_distance])

    stats_df = pd.DataFrame(stats_data, columns=["distance", "count", "all pairs", "average distance"])
    return stats_df

def write_statistics_to_excel(stats_file, stats_df):
    stats_df.to_excel(stats_file, index=False)

def main():
    #input and output files
    bam_file = "file.bam"
    max_distance = NUMBER, int #maximum distance between reads
    pairs_output_file = " "
    discarded_output_file = " "
    stats_output_file = " "

    plus_positions = extract_positions_by_strand(bam_file, "+")
    minus_positions = extract_positions_by_strand(bam_file, "-")
    pairs = find_all_pairs_within_distance(plus_positions, minus_positions, max_distance)

    write_pairs_to_excel(pairs_output_file, pairs)
    discarded_pairs = []
    write_discarded_pairs_to_excel(discarded_output_file, discarded_pairs)

    #statistics
    stats_df = calculate_statistics(pairs, max_distance)
    write_statistics_to_excel(stats_output_file, stats_df)

if __name__ == "__main__":
    main()
