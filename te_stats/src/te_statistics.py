#!/usr/bin/env python3
"""
TE Statistics Module for Transposable Element Statistical Analysis
This module provides comprehensive statistical analysis capabilities for transposable elements.
"""

import subprocess
import tempfile
import os
from typing import Dict, List, Tuple, Optional, Any
from dataclasses import dataclass, asdict
from collections import defaultdict, Counter
import pandas as pd
import numpy as np
import logging
from pathlib import Path

from .gff_parser import TERecord, GFFParser
from .classifier import TEClassifier, ClassLevel

logger = logging.getLogger(__name__)


@dataclass
class TESummary:
    """Data class representing TE summary statistics for a category."""
    level1: str
    level2: str
    level3: str
    element_count: int
    raw_length: int
    merged_length: int
    length_median: float
    length_mean: float
    length_std: float
    chromosome_count: int
    density_per_mb: float = 0.0
    genome_percentage: float = 0.0
    compression_ratio: float = 0.0  # merged_length / raw_length


class TEStatistics:
    """
    Statistical analysis engine for transposable elements.

    This class provides comprehensive statistical analysis including:
    - Element counting and length calculations
    - Interval merging using bedtools
    - Density and percentage calculations
    - Multi-level hierarchical statistics
    """

    def __init__(self, te_records: List[TERecord], classifier: TEClassifier,
                 genome_size: Optional[int] = None, bedtools_path: str = "bedtools"):
        """
        Initialize the TE statistics analyzer.

        Args:
            te_records (List[TERecord]): List of TE records
            classifier (TEClassifier): TE classifier instance
            genome_size (Optional[int]): Total genome size in base pairs
            bedtools_path (str): Path to bedtools executable
        """
        self.te_records = te_records
        self.classifier = classifier
        self.genome_size = genome_size
        self.bedtools_path = bedtools_path

        # Verify bedtools availability
        self._check_bedtools()

        # Initialize storage for results
        self.classified_records: Dict[Tuple[str, str, str], List[TERecord]] = defaultdict(list)
        self.merged_intervals: Dict[Tuple[str, str, str], List[Tuple[str, int, int]]] = {}
        self.summary_stats: Dict[Tuple[str, str, str], TESummary] = {}

        # Classify all records
        self._classify_records()

        logger.info(f"Initialized TE statistics with {len(te_records)} records")

    def _check_bedtools(self) -> None:
        """Check if bedtools is available and accessible."""
        try:
            result = subprocess.run([self.bedtools_path, "--version"],
                                  capture_output=True, text=True, timeout=10)
            if result.returncode == 0:
                logger.info(f"Using bedtools: {result.stdout.strip()}")
            else:
                raise RuntimeError("bedtools check failed")
        except (subprocess.TimeoutExpired, FileNotFoundError, RuntimeError) as e:
            logger.error(f"bedtools not found or not accessible: {e}")
            logger.warning("Interval merging functionality will be disabled")
            self.bedtools_path = None

    def _classify_records(self) -> None:
        """Classify all TE records into hierarchical categories."""
        logger.info("Classifying TE records...")

        for record in self.te_records:
            # Get hierarchical classification
            level1, level2, level3 = self.classifier.classify_element(record.te_class)
            classification_key = (level1, level2, level3)

            self.classified_records[classification_key].append(record)

        logger.info(f"Classified {len(self.te_records)} records into "
                   f"{len(self.classified_records)} categories")

    def _merge_intervals(self, records: List[TERecord]) -> List[Tuple[str, int, int]]:
        """
        Merge overlapping intervals using bedtools.

        Args:
            records (List[TERecord]): List of TE records to merge

        Returns:
            List[Tuple[str, int, int]]: Merged intervals (chrom, start, end)
        """
        if not records or not self.bedtools_path:
            return []

        # Sort records by chromosome and start position for bedtools
        sorted_records = sorted(records, key=lambda x: (x.seqid, x.start))

        # Create temporary files for bedtools
        with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as input_file:
            # Write records to BED format (0-based)
            for record in sorted_records:
                bed_start = record.start - 1  # Convert to 0-based
                bed_end = record.end
                input_file.write(f"{record.seqid}\t{bed_start}\t{bed_end}\n")

            input_file_path = input_file.name

        try:
            # Run bedtools merge
            cmd = [self.bedtools_path, "merge", "-i", input_file_path]
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)

            if result.returncode != 0:
                logger.error(f"bedtools merge failed: {result.stderr}")
                return []

            # Parse merged intervals
            merged_intervals = []
            for line in result.stdout.strip().split('\n'):
                if line:
                    parts = line.split('\t')
                    if len(parts) >= 3:
                        chrom = parts[0]
                        start = int(parts[1])
                        end = int(parts[2])
                        merged_intervals.append((chrom, start, end))

            return merged_intervals

        except subprocess.TimeoutExpired:
            logger.error("bedtools merge timed out")
            return []
        except Exception as e:
            logger.error(f"Error during interval merging: {e}")
            return []
        finally:
            # Clean up temporary file
            try:
                os.unlink(input_file_path)
            except OSError:
                pass

    def _calculate_basic_statistics(self, records: List[TERecord]) -> Dict[str, float]:
        """
        Calculate basic statistics for a set of records.

        Args:
            records (List[TERecord]): List of TE records

        Returns:
            Dict[str, float]: Basic statistics
        """
        if not records:
            return {
                'count': 0,
                'total_length': 0,
                'mean_length': 0,
                'median_length': 0,
                'std_length': 0
            }

        lengths = [record.length for record in records]

        return {
            'count': len(records),
            'total_length': sum(lengths),
            'mean_length': np.mean(lengths),
            'median_length': np.median(lengths),
            'std_length': np.std(lengths)
        }

    def calculate_category_statistics(self, classification_key: Tuple[str, str, str]) -> TESummary:
        """
        Calculate comprehensive statistics for a specific category.

        Args:
            classification_key (Tuple[str, str, str]): (level1, level2, level3) category

        Returns:
            TESummary: Summary statistics for the category
        """
        records = self.classified_records[classification_key]
        level1, level2, level3 = classification_key

        if not records:
            return TESummary(
                level1=level1, level2=level2, level3=level3,
                element_count=0, raw_length=0, merged_length=0,
                length_median=0, length_mean=0, length_std=0,
                chromosome_count=0
            )

        # Basic statistics
        basic_stats = self._calculate_basic_statistics(records)

        # Merge intervals (if bedtools is available)
        if self.bedtools_path:
            if classification_key not in self.merged_intervals:
                self.merged_intervals[classification_key] = self._merge_intervals(records)

            merged_intervals = self.merged_intervals[classification_key]
            merged_length = sum(end - start for chrom, start, end in merged_intervals)
        else:
            merged_intervals = []
            merged_length = basic_stats['total_length']  # Fallback to raw length

        # Count unique chromosomes
        chromosomes = set(record.seqid for record in records)

        # Calculate additional metrics
        density_per_mb = 0.0
        genome_percentage = 0.0
        compression_ratio = 0.0

        if self.genome_size:
            density_per_mb = (len(records) / self.genome_size) * 1_000_000
            genome_percentage = (merged_length / self.genome_size) * 100

        if basic_stats['total_length'] > 0:
            compression_ratio = merged_length / basic_stats['total_length']

        return TESummary(
            level1=level1,
            level2=level2,
            level3=level3,
            element_count=basic_stats['count'],
            raw_length=basic_stats['total_length'],
            merged_length=merged_length,
            length_median=basic_stats['median_length'],
            length_mean=basic_stats['mean_length'],
            length_std=basic_stats['std_length'],
            chromosome_count=len(chromosomes),
            density_per_mb=density_per_mb,
            genome_percentage=genome_percentage,
            compression_ratio=compression_ratio
        )

    def calculate_all_statistics(self) -> Dict[Tuple[str, str, str], TESummary]:
        """
        Calculate statistics for all categories.

        Returns:
            Dict[Tuple[str, str, str], TESummary]: Statistics for all categories
        """
        logger.info("Calculating statistics for all categories...")

        for classification_key in self.classified_records.keys():
            self.summary_stats[classification_key] = self.calculate_category_statistics(classification_key)

        logger.info(f"Calculated statistics for {len(self.summary_stats)} categories")
        return self.summary_stats

    def get_hierarchical_summary(self, level: ClassLevel = ClassLevel.LEVEL1) -> pd.DataFrame:
        """
        Get hierarchical summary statistics.

        Args:
            level (ClassLevel): Aggregation level

        Returns:
            pd.DataFrame: Hierarchical summary statistics
        """
        if not self.summary_stats:
            self.calculate_all_statistics()

        # Convert summaries to list of dictionaries
        summary_list = []
        for key, summary in self.summary_stats.items():
            summary_dict = asdict(summary)
            summary_list.append(summary_dict)

        df = pd.DataFrame(summary_list)

        if df.empty:
            return df

        # Aggregate based on the specified level
        if level == ClassLevel.LEVEL1:
            group_cols = ['level1']
        elif level == ClassLevel.LEVEL2:
            group_cols = ['level1', 'level2']
        else:
            group_cols = ['level1', 'level2', 'level3']

        agg_functions = {
            'element_count': 'sum',
            'raw_length': 'sum',
            'merged_length': 'sum',
            'length_median': 'mean',
            'length_mean': 'mean',
            'length_std': 'mean',
            'chromosome_count': 'max',
            'density_per_mb': 'sum',
            'genome_percentage': 'sum',
            'compression_ratio': 'mean'
        }

        return df.groupby(group_cols).agg(agg_functions).reset_index()

    def get_top_categories(self, metric: str = 'merged_length', level: ClassLevel = ClassLevel.LEVEL3,
                          top_n: int = 10) -> pd.DataFrame:
        """
        Get top categories by a specific metric.

        Args:
            metric (str): Metric to sort by (element_count, raw_length, merged_length, etc.)
            level (ClassLevel): Aggregation level
            top_n (int): Number of top categories to return

        Returns:
            pd.DataFrame: Top categories by the specified metric
        """
        df = self.get_hierarchical_summary(level)

        if metric not in df.columns:
            raise ValueError(f"Unknown metric: {metric}")

        return df.nlargest(top_n, metric)

    def get_chromosome_distribution(self) -> Dict[str, Dict[str, Any]]:
        """
        Get TE distribution across chromosomes.

        Returns:
            Dict[str, Dict[str, Any]]: Chromosome-wise TE distribution
        """
        chromosome_stats = defaultdict(lambda: {
            'total_count': 0,
            'total_length': 0,
            'categories': defaultdict(lambda: {'count': 0, 'length': 0})
        })

        for record in self.te_records:
            chrom = record.seqid
            level1, level2, level3 = self.classifier.classify_element(record.te_class)
            category_key = f"{level1}/{level2}"

            chromosome_stats[chrom]['total_count'] += 1
            chromosome_stats[chrom]['total_length'] += record.length
            chromosome_stats[chrom]['categories'][category_key]['count'] += 1
            chromosome_stats[chrom]['categories'][category_key]['length'] += record.length

        return dict(chromosome_stats)

    def get_length_distribution(self, classification_key: Optional[Tuple[str, str, str]] = None) -> Dict:
        """
        Get length distribution statistics.

        Args:
            classification_key (Optional[Tuple[str, str, str]]): Specific category to analyze,
                                                               or None for all TEs

        Returns:
            Dict: Length distribution statistics
        """
        if classification_key:
            records = self.classified_records[classification_key]
        else:
            records = self.te_records

        if not records:
            return {}

        lengths = [record.length for record in records]

        # Calculate length distribution bins
        bins = [0, 100, 500, 1000, 2000, 5000, 10000, float('inf')]
        bin_labels = ['<100', '100-500', '500-1K', '1K-2K', '2K-5K', '5K-10K', '>10K']

        hist, _ = np.histogram(lengths, bins=bins)
        distribution = dict(zip(bin_labels, hist.tolist()))

        return {
            'total_elements': len(lengths),
            'min_length': min(lengths),
            'max_length': max(lengths),
            'mean_length': np.mean(lengths),
            'median_length': np.median(lengths),
            'std_length': np.std(lengths),
            'distribution_bins': distribution
        }

    def export_results(self, output_dir: str, formats: List[str] = ['csv']) -> None:
        """
        Export statistical results to files.

        Args:
            output_dir (str): Output directory path
            formats (List[str]): Export formats (csv, tsv, excel)
        """
        if not self.summary_stats:
            self.calculate_all_statistics()

        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        # Export detailed statistics
        summary_list = []
        for key, summary in self.summary_stats.items():
            summary_dict = asdict(summary)
            summary_dict['category_key'] = '/'.join(key)
            summary_list.append(summary_dict)

        df_detailed = pd.DataFrame(summary_list)

        # Export in requested formats
        for format_name in formats:
            if format_name == 'csv':
                df_detailed.to_csv(output_path / 'te_statistics_detailed.csv', index=False)
            elif format_name == 'tsv':
                df_detailed.to_csv(output_path / 'te_statistics_detailed.tsv', index=False, sep='\t')
            elif format_name == 'excel':
                with pd.ExcelWriter(output_path / 'te_statistics_report.xlsx') as writer:
                    df_detailed.to_excel(writer, sheet_name='Detailed Statistics', index=False)

                    # Add hierarchical summaries
                    for level in [ClassLevel.LEVEL1, ClassLevel.LEVEL2, ClassLevel.LEVEL3]:
                        df_hier = self.get_hierarchical_summary(level)
                        sheet_name = f'Hierarchy_{level.value.replace(" ", "_")}'
                        df_hier.to_excel(writer, sheet_name=sheet_name, index=False)

        # Export chromosome distribution
        chrom_dist = self.get_chromosome_distribution()
        chrom_df = []
        for chrom, stats in chrom_dist.items():
            chrom_df.append({
                'chromosome': chrom,
                'total_count': stats['total_count'],
                'total_length': stats['total_length']
            })

        if chrom_df:
            pd.DataFrame(chrom_df).to_csv(output_path / 'chromosome_distribution.csv', index=False)

        logger.info(f"Exported results to {output_dir} in formats: {formats}")

    def print_summary(self) -> None:
        """Print a formatted summary of the statistical analysis."""
        if not self.summary_stats:
            self.calculate_all_statistics()

        print("\n" + "="*80)
        print("TRANSPOSABLE ELEMENT STATISTICAL ANALYSIS SUMMARY")
        print("="*80)

        if self.genome_size:
            print(f"Genome size: {self.genome_size:,} bp")
        print(f"Total TE records analyzed: {len(self.te_records):,}")
        print(f"Unique categories found: {len(self.summary_stats)}")
        print(f"Analysis date: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}")
        print("-"*80)

        # Print Level 1 summary
        level1_summary = self.get_hierarchical_summary(ClassLevel.LEVEL1)
        if not level1_summary.empty:
            print("\nLEVEL 1 (CLASS LEVEL) SUMMARY:")
            level1_summary = level1_summary.sort_values('merged_length', ascending=False)
            for _, row in level1_summary.iterrows():
                print(f"  {row['level1']:12} | "
                      f"Count: {row['element_count']:6,} | "
                      f"Length: {row['merged_length']:10,} bp | "
                      f"Genome: {row['genome_percentage']:5.2f}%")

        # Print top Level 3 categories
        print("\nTOP 10 LEVEL 3 (FAMILY LEVEL) CATEGORIES BY LENGTH:")
        top_categories = self.get_top_categories('merged_length', ClassLevel.LEVEL3, 10)
        for _, row in top_categories.iterrows():
            print(f"  {row['level1']}/{row['level2']}/{row['level3'][:20]:20} | "
                  f"Count: {row['element_count']:5,} | "
                  f"Length: {row['merged_length']:10,} bp | "
                  f"Mean: {row['length_mean']:6.0f} bp")

        print("\n" + "="*80)


def main():
    """Example usage of TEStatistics."""
    from .gff_parser import GFFParser
    from .classifier import TEClassifier

    # Example with test data
    parser = GFFParser("../test.gff")
    records = parser.parse_gff_file()

    classifier = TEClassifier()
    stats = TEStatistics(records, classifier, genome_size=16_000_000)

    # Calculate all statistics
    stats.calculate_all_statistics()

    # Print summary
    stats.print_summary()

    # Export results
    stats.export_results("results", ['csv', 'excel'])


if __name__ == "__main__":
    main()