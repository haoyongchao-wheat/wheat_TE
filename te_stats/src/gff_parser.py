#!/usr/bin/env python3
"""
GFF Parser Module for Transposable Element Analysis
This module provides functionality to parse GFF3 files and extract transposable element information.
"""

import re
import gzip
from typing import Dict, List, Tuple, Optional, Iterator
from dataclasses import dataclass
from pathlib import Path
import logging

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


@dataclass
class TERecord:
    """Data class representing a transposable element record from GFF3."""
    seqid: str
    source: str
    feature_type: str
    start: int
    end: int
    score: float
    strand: str
    phase: str
    attributes: Dict[str, str]

    def __post_init__(self):
        """Validate and convert numeric fields."""
        self.start = int(self.start)
        self.end = int(self.end)
        self.score = float(self.score) if self.score != '.' else 0.0

    @property
    def length(self) -> int:
        """Calculate the length of the TE record."""
        return self.end - self.start + 1

    @property
    def te_class(self) -> str:
        """Extract TE class from attributes."""
        return self.attributes.get('Class', 'Unknown')

    @property
    def te_target(self) -> str:
        """Extract TE target from attributes."""
        return self.attributes.get('Target', 'Unknown')

    @property
    def te_id(self) -> str:
        """Extract TE ID from attributes."""
        return self.attributes.get('ID', 'Unknown')


class GFFParser:
    """
    Parser for GFF3 files containing transposable element annotations.

    This class handles reading and parsing GFF3 format files, extracting
    transposable element records and their associated attributes.
    """

    def __init__(self, file_path: str):
        """
        Initialize the GFF parser.

        Args:
            file_path (str): Path to the GFF3 file
        """
        self.file_path = Path(file_path)
        self.te_records: List[TERecord] = []
        self.class_counts: Dict[str, int] = {}

        if not self.file_path.exists():
            raise FileNotFoundError(f"GFF3 file not found: {file_path}")

        logger.info(f"Initialized GFF parser for file: {file_path}")

    def _is_gzipped(self) -> bool:
        """Check if the file is gzipped."""
        return self.file_path.suffix.lower() == '.gz'

    def _open_file(self):
        """Open file with appropriate decompression if needed."""
        if self._is_gzipped():
            return gzip.open(self.file_path, 'rt', encoding='utf-8')
        else:
            return open(self.file_path, 'r', encoding='utf-8')

    def _parse_attributes(self, attr_string: str) -> Dict[str, str]:
        """
        Parse GFF3 attribute string into dictionary.

        Args:
            attr_string (str): Attribute string from GFF3 column 9

        Returns:
            Dict[str, str]: Parsed attributes as key-value pairs
        """
        attributes = {}
        if not attr_string or attr_string == '.':
            return attributes

        # Split by semicolon and clean up
        parts = [part.strip() for part in attr_string.split(';') if part.strip()]

        for part in parts:
            if '=' in part:
                key, value = part.split('=', 1)
                # Clean up the value
                value = value.strip('"\'')
                attributes[key] = value
            else:
                # Handle case where there might be no equals sign
                attributes[part] = ''

        return attributes

    def _parse_line(self, line: str) -> Optional[TERecord]:
        """
        Parse a single GFF3 line into a TERecord object.

        Args:
            line (str): Single line from GFF3 file

        Returns:
            Optional[TERecord]: Parsed TE record or None if not a valid record
        """
        line = line.strip()

        # Skip comments and empty lines
        if not line or line.startswith('#'):
            return None

        fields = line.split('\t')

        # Check if we have the right number of fields
        if len(fields) != 9:
            logger.warning(f"Invalid GFF3 line (wrong number of fields): {line}")
            return None

        seqid, source, feature_type, start, end, score, strand, phase, attributes = fields

        # Only process transposon features (allow variations in naming)
        te_feature_types = [
            'Transposon', 'transposon', 'TE', 'transposable_element',
            'transposable_element_gene', 'transposable_element_fragment',
            'retrotransposon', 'DNA_transposon', 'mobile_element',
            'insertion_sequence', 'IS_element'
        ]
        if feature_type not in te_feature_types:
            return None

        try:
            # Parse attributes
            attr_dict = self._parse_attributes(attributes)

            # Create TERecord
            record = TERecord(
                seqid=seqid,
                source=source,
                feature_type=feature_type,
                start=start,
                end=end,
                score=score,
                strand=strand,
                phase=phase,
                attributes=attr_dict
            )

            return record

        except (ValueError, IndexError) as e:
            logger.warning(f"Error parsing line: {line}. Error: {e}")
            return None

    def parse_gff_file(self) -> List[TERecord]:
        """
        Parse the entire GFF3 file and extract TE records.

        Returns:
            List[TERecord]: List of transposable element records
        """
        logger.info(f"Starting to parse GFF3 file: {self.file_path}")

        self.te_records = []
        self.class_counts = {}

        with self._open_file() as file:
            for line_num, line in enumerate(file, 1):
                record = self._parse_line(line)

                if record:
                    self.te_records.append(record)

                    # Update class counts
                    te_class = record.te_class
                    self.class_counts[te_class] = self.class_counts.get(te_class, 0) + 1

                    # Progress reporting for large files
                    if line_num % 10000 == 0:
                        logger.info(f"Processed {line_num} lines, found {len(self.te_records)} TE records")

        logger.info(f"Finished parsing. Found {len(self.te_records)} TE records")
        logger.info(f"Found {len(self.class_counts)} different TE classes")

        return self.te_records

    def get_te_records(self) -> List[TERecord]:
        """
        Get all parsed TE records.

        Returns:
            List[TERecord]: List of TE records
        """
        return self.te_records

    def get_records_by_class(self, te_class: str) -> List[TERecord]:
        """
        Get TE records filtered by class.

        Args:
            te_class (str): TE class to filter by

        Returns:
            List[TERecord]: List of TE records of the specified class
        """
        return [record for record in self.te_records if record.te_class == te_class]

    def get_records_by_chromosome(self, chromosome: str) -> List[TERecord]:
        """
        Get TE records filtered by chromosome.

        Args:
            chromosome (str): Chromosome name to filter by

        Returns:
            List[TERecord]: List of TE records on the specified chromosome
        """
        return [record for record in self.te_records if record.seqid == chromosome]

    def get_chromosome_list(self) -> List[str]:
        """
        Get list of all chromosomes present in the data.

        Returns:
            List[str]: List of chromosome names
        """
        return sorted(list(set(record.seqid for record in self.te_records)))

    def get_class_distribution(self) -> Dict[str, int]:
        """
        Get distribution of TE classes.

        Returns:
            Dict[str, int]: Dictionary with class names and their counts
        """
        return self.class_counts.copy()

    def get_summary_statistics(self) -> Dict:
        """
        Get basic summary statistics of the parsed data.

        Returns:
            Dict: Summary statistics
        """
        if not self.te_records:
            return {}

        total_length = sum(record.length for record in self.te_records)
        chromosomes = self.get_chromosome_list()

        return {
            'total_records': len(self.te_records),
            'total_length': total_length,
            'average_length': total_length / len(self.te_records),
            'num_classes': len(self.class_counts),
            'num_chromosomes': len(chromosomes),
            'chromosomes': chromosomes,
            'class_distribution': self.class_counts
        }

    def to_bed_format(self, records: Optional[List[TERecord]] = None) -> List[Tuple[str, int, int, str]]:
        """
        Convert TE records to BED format.

        Args:
            records (Optional[List[TERecord]]): Records to convert. If None, uses all records.

        Returns:
            List[Tuple[str, int, int, str]]: BED format records (chrom, start, end, name)
        """
        if records is None:
            records = self.te_records

        bed_records = []
        for record in records:
            # Convert to 0-based coordinates for BED format
            bed_start = record.start - 1
            bed_end = record.end
            name = f"{record.te_id}|{record.te_class}"

            bed_records.append((record.seqid, bed_start, bed_end, name))

        return bed_records

    def write_bed_file(self, output_path: str, records: Optional[List[TERecord]] = None) -> None:
        """
        Write TE records to BED file.

        Args:
            output_path (str): Output BED file path
            records (Optional[List[TERecord]]): Records to write. If None, writes all records.
        """
        bed_records = self.to_bed_format(records)

        with open(output_path, 'w') as f:
            for chrom, start, end, name in bed_records:
                f.write(f"{chrom}\t{start}\t{end}\t{name}\n")

        logger.info(f"Wrote {len(bed_records)} records to BED file: {output_path}")


def main():
    """Example usage of GFFParser."""
    import sys

    if len(sys.argv) != 2:
        print("Usage: python gff_parser.py <gff3_file>")
        sys.exit(1)

    gff_file = sys.argv[1]

    try:
        parser = GFFParser(gff_file)
        records = parser.parse_gff_file()

        # Print summary statistics
        stats = parser.get_summary_statistics()
        print("\n=== GFF3 File Summary ===")
        print(f"Total TE records: {stats['total_records']}")
        print(f"Total length: {stats['total_length']:,} bp")
        print(f"Average length: {stats['average_length']:.1f} bp")
        print(f"Number of classes: {stats['num_classes']}")
        print(f"Number of chromosomes: {stats['num_chromosomes']}")

        print("\n=== Class Distribution ===")
        for te_class, count in sorted(stats['class_distribution'].items(),
                                    key=lambda x: x[1], reverse=True):
            print(f"{te_class}: {count}")

    except Exception as e:
        logger.error(f"Error processing file: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()