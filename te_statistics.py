#!/usr/bin/env python3
"""
TE Statistics - Comprehensive Transposable Element Analysis Tool

This script provides a complete pipeline for analyzing transposable elements
from GFF3 annotation files, including hierarchical classification, statistical
analysis, and visualization.

Author: Bioinformatics Team
Version: 1.0.0
"""

import argparse
import sys
import os
import yaml
import logging
from pathlib import Path
from typing import Dict, List, Optional, Any
from datetime import datetime

# Add src directory to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'te_stats', 'src'))

from te_stats.src.gff_parser import GFFParser
from te_stats.src.classifier import TEClassifier, ClassLevel
from te_stats.src.te_statistics import TEStatistics
from te_stats.src.visualizer import TEVisualizer

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('te_statistics.log'),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)


def parse_arguments() -> argparse.Namespace:
    """
    Parse command line arguments.

    Returns:
        argparse.Namespace: Parsed arguments
    """
    parser = argparse.ArgumentParser(
        description='Comprehensive Transposable Element Analysis Tool',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic analysis
  python te_statistics.py -i input.gff -o results/

  # Analysis with genome size and multiple output formats
  python te_statistics.py -i input.gff -o results/ -g 16000000 --format csv excel

  # Generate plots only
  python te_statistics.py -i input.gff -o plots/ --plot-only --plot-type pie bar

  # Use custom configuration
  python te_statistics.py -i input.gff -o results/ --config custom_config.yaml

  # Verbose output with additional statistics
  python te_statistics.py -i input.gff -o results/ --verbose --merge-intervals
        """
    )

    # Input/Output arguments
    parser.add_argument('-i', '--input', required=True,
                       help='Input GFF3 file containing TE annotations')
    parser.add_argument('-o', '--output', required=True,
                       help='Output directory for results')
    parser.add_argument('-c', '--config',
                       help='Configuration file (YAML format)')

    # Analysis options
    parser.add_argument('-g', '--genome-size', type=int,
                       help='Genome size in base pairs (for percentage calculations)')
    parser.add_argument('--min-length', type=int, default=50,
                       help='Minimum TE length to include in analysis (default: 50)')
    parser.add_argument('--exclude-unspecified', action='store_true',
                       help='Exclude unspecified/unknown TEs from analysis')
    parser.add_argument('--merge-intervals', action='store_true', default=True,
                       help='Merge overlapping intervals using bedtools (default: True)')
    parser.add_argument('--bedtools-path', default='bedtools',
                       help='Path to bedtools executable (default: bedtools)')

    # Output format options
    parser.add_argument('--format', nargs='+', default=['csv'],
                       choices=['csv', 'tsv', 'excel'],
                       help='Output format(s) for statistics (default: csv)')
    parser.add_argument('--plot-only', action='store_true',
                       help='Only generate plots, skip statistical analysis')

    # Plot options
    parser.add_argument('--plot-type', nargs='+', default=['all'],
                       choices=['all', 'pie', 'bar', 'length', 'chromosome', 'dashboard'],
                       help='Types of plots to generate (default: all)')
    parser.add_argument('--plot-format', nargs='+', default=['png'],
                       choices=['png', 'svg', 'pdf'],
                       help='Plot output format(s) (default: png)')
    parser.add_argument('--style', default='seaborn-v0_8',
                       help='Matplotlib style for plots (default: seaborn-v0_8)')

    # Other options
    parser.add_argument('--verbose', action='store_true',
                       help='Enable verbose output')
    parser.add_argument('--quiet', action='store_true',
                       help='Suppress non-error output')
    parser.add_argument('--version', action='version', version='TE Statistics 1.0.0')

    return parser.parse_args()


def load_config(config_path: str) -> Dict[str, Any]:
    """
    Load configuration from YAML file.

    Args:
        config_path (str): Path to configuration file

    Returns:
        Dict[str, Any]: Configuration dictionary
    """
    try:
        with open(config_path, 'r') as f:
            config = yaml.safe_load(f)
        logger.info(f"Loaded configuration from {config_path}")
        return config
    except Exception as e:
        logger.error(f"Error loading configuration file {config_path}: {e}")
        sys.exit(1)


def create_default_config(output_path: str) -> str:
    """
    Create a default configuration file.

    Args:
        output_path (str): Path to create the configuration file

    Returns:
        str: Path to created configuration file
    """
    default_config = {
        'classification': {
            'level1': {
                'Class I': ['LTR', 'LINE', 'SINE', 'DIRS'],
                'Class II': ['TIR', 'Helitron', 'Maverick', 'DNA'],
                'Other': ['Unspecified', 'Unknown']
            }
        },
        'analysis': {
            'merge_overlaps': True,
            'min_length': 50,
            'exclude_unspecified': False
        },
        'output': {
            'formats': ['csv', 'tsv', 'excel'],
            'create_plots': True,
            'plot_types': ['pie', 'bar', 'length', 'chromosome', 'dashboard']
        },
        'visualization': {
            'style': 'seaborn-v0_8',
            'dpi': 300,
            'figure_format': ['png', 'svg']
        }
    }

    config_path = Path(output_path) / 'config.yaml'
    with open(config_path, 'w') as f:
        yaml.dump(default_config, f, default_flow_style=False, indent=2)

    logger.info(f"Created default configuration file: {config_path}")
    return str(config_path)


def setup_logging(verbose: bool, quiet: bool) -> None:
    """
    Setup logging based on verbosity options.

    Args:
        verbose (bool): Enable verbose logging
        quiet (bool): Suppress non-error output
    """
    if quiet:
        logging.getLogger().setLevel(logging.ERROR)
    elif verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    else:
        logging.getLogger().setLevel(logging.INFO)


def validate_inputs(args: argparse.Namespace) -> None:
    """
    Validate input arguments and files.

    Args:
        args (argparse.Namespace): Parsed arguments
    """
    # Check input file
    input_path = Path(args.input)
    if not input_path.exists():
        logger.error(f"Input file not found: {args.input}")
        sys.exit(1)

    if not input_path.suffix.lower() in ['.gff', '.gff3', '.gz']:
        logger.warning(f"Input file does not have expected GFF extension: {args.input}")

    # Check/create output directory
    output_path = Path(args.output)
    try:
        output_path.mkdir(parents=True, exist_ok=True)
    except Exception as e:
        logger.error(f"Cannot create output directory {args.output}: {e}")
        sys.exit(1)

    # Check bedtools if interval merging is requested
    if args.merge_intervals:
        try:
            import subprocess
            result = subprocess.run([args.bedtools_path, '--version'],
                                  capture_output=True, text=True, timeout=5)
            if result.returncode != 0:
                raise RuntimeError("bedtools not found")
        except Exception as e:
            logger.warning(f"bedtools not available: {e}. Interval merging will be disabled.")
            args.merge_intervals = False


def run_analysis(args: argparse.Namespace, config: Optional[Dict] = None) -> None:
    """
    Run the complete TE analysis pipeline.

    Args:
        args (argparse.Namespace): Parsed arguments
        config (Optional[Dict]): Configuration dictionary
    """
    logger.info("Starting TE analysis pipeline")
    start_time = datetime.now()

    try:
        # Step 1: Parse GFF file
        logger.info(f"Parsing GFF file: {args.input}")
        parser = GFFParser(args.input)
        te_records = parser.parse_gff_file()

        # Filter by minimum length
        if args.min_length > 0:
            original_count = len(te_records)
            te_records = [record for record in te_records if record.length >= args.min_length]
            filtered_count = original_count - len(te_records)
            if filtered_count > 0:
                logger.info(f"Filtered out {filtered_count} TEs shorter than {args.min_length} bp")

        # Exclude unspecified if requested
        if args.exclude_unspecified:
            original_count = len(te_records)
            te_records = [record for record in te_records
                         if not any(keyword in record.te_class.lower()
                                  for keyword in ['unspecified', 'unknown'])]
            excluded_count = original_count - len(te_records)
            if excluded_count > 0:
                logger.info(f"Excluded {excluded_count} unspecified/unknown TEs")

        logger.info(f"Processing {len(te_records)} TE records")

        # Step 2: Initialize classifier
        logger.info("Initializing TE classifier")
        custom_classification = config.get('classification', {}) if config else None
        classifier = TEClassifier(custom_classification)

        # Print classification summary if verbose
        if args.verbose:
            class_summary = classifier.get_classification_summary()
            logger.info(f"Classification system: {class_summary['num_level1_categories']} level 1, "
                       f"{class_summary['num_level2_categories']} level 2, "
                       f"{class_summary['total_classification_rules']} total rules")

        # Step 3: Initialize statistics engine
        logger.info("Initializing statistics engine")
        stats = TEStatistics(
            te_records=te_records,
            classifier=classifier,
            genome_size=args.genome_size,
            bedtools_path=args.bedtools_path
        )

        # Step 4: Calculate statistics (unless plot-only)
        if not args.plot_only:
            logger.info("Calculating statistical analysis")
            stats.calculate_all_statistics()

            # Print summary
            stats.print_summary()

            # Export statistical results
            logger.info(f"Exporting results to {args.output}")
            stats.export_results(args.output, args.format)

            # Export classification rules
            if args.verbose:
                rules_path = Path(args.output) / 'classification_rules.csv'
                classifier.export_classification_rules(str(rules_path))

        # Step 5: Generate plots
        if 'all' in args.plot_type or args.plot_only:
            plot_types = ['pie', 'bar', 'length', 'chromosome', 'dashboard']
        else:
            plot_types = args.plot_type

        if plot_types:
            logger.info(f"Generating {len(plot_types)} types of plots")
            visualizer = TEVisualizer(stats, style=args.style)

            plots_dir = Path(args.output) / 'plots'
            plots_dir.mkdir(exist_ok=True)

            # Generate individual plots
            for plot_type in plot_types:
                try:
                    if plot_type == 'pie':
                        fig = visualizer.plot_composition_pie(ClassLevel.LEVEL1)
                        for fmt in args.plot_format:
                            fig.savefig(plots_dir / f'composition_class1.{fmt}',
                                      dpi=300, bbox_inches='tight')
                        fig = visualizer.plot_composition_pie(ClassLevel.LEVEL2)
                        for fmt in args.plot_format:
                            fig.savefig(plots_dir / f'composition_class2.{fmt}',
                                      dpi=300, bbox_inches='tight')

                    elif plot_type == 'bar':
                        fig = visualizer.plot_comparison_bar(ClassLevel.LEVEL2)
                        for fmt in args.plot_format:
                            fig.savefig(plots_dir / f'comparison_class2.{fmt}',
                                      dpi=300, bbox_inches='tight')

                    elif plot_type == 'length':
                        fig = visualizer.plot_length_distribution()
                        for fmt in args.plot_format:
                            fig.savefig(plots_dir / f'length_distribution.{fmt}',
                                      dpi=300, bbox_inches='tight')

                    elif plot_type == 'chromosome':
                        fig = visualizer.plot_chromosome_distribution()
                        for fmt in args.plot_format:
                            fig.savefig(plots_dir / f'chromosome_distribution.{fmt}',
                                      dpi=300, bbox_inches='tight')

                    elif plot_type == 'dashboard':
                        fig = visualizer.plot_summary_dashboard()
                        for fmt in args.plot_format:
                            fig.savefig(plots_dir / f'summary_dashboard.{fmt}',
                                      dpi=300, bbox_inches='tight')

                    import matplotlib.pyplot as plt
                    plt.close('all')  # Close figures to free memory

                except Exception as e:
                    logger.error(f"Error generating {plot_type} plot: {e}")

        # Calculate execution time
        end_time = datetime.now()
        execution_time = end_time - start_time
        logger.info(f"Analysis completed successfully in {execution_time}")

        # Save execution summary
        summary_path = Path(args.output) / 'execution_summary.txt'
        with open(summary_path, 'w') as f:
            f.write("TE Statistics Analysis Summary\n")
            f.write("=" * 40 + "\n")
            f.write(f"Input file: {args.input}\n")
            f.write(f"Output directory: {args.output}\n")
            f.write(f"Analysis date: {start_time.strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Execution time: {execution_time}\n")
            f.write(f"TE records processed: {len(te_records)}\n")
            if args.genome_size:
                f.write(f"Genome size: {args.genome_size:,} bp\n")
            f.write(f"Merge intervals: {args.merge_intervals}\n")
            f.write(f"Min length filter: {args.min_length} bp\n")
            f.write(f"Output formats: {', '.join(args.format)}\n")
            f.write(f"Plot formats: {', '.join(args.plot_format)}\n")

    except Exception as e:
        logger.error(f"Analysis failed: {e}")
        if args.verbose:
            import traceback
            logger.error(traceback.format_exc())
        sys.exit(1)


def main():
    """Main function."""
    # Parse arguments
    args = parse_arguments()

    # Setup logging
    setup_logging(args.verbose, args.quiet)

    # Validate inputs
    validate_inputs(args)

    # Load configuration if provided
    config = None
    if args.config:
        config = load_config(args.config)
    elif args.verbose:
        # Create default config for reference
        create_default_config(args.output)

    # Run analysis
    run_analysis(args, config)


if __name__ == "__main__":
    main()