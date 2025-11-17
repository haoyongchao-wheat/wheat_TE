#!/usr/bin/env python3
"""
TE Visualizer Module for Transposable Element Data Visualization
This module provides comprehensive visualization capabilities for TE analysis results.
"""

import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import Patch
import seaborn as sns
import pandas as pd
import numpy as np
from typing import Dict, List, Tuple, Optional, Any
from pathlib import Path
import logging
from dataclasses import asdict

from .te_statistics import TEStatistics, TESummary
from .classifier import ClassLevel

logger = logging.getLogger(__name__)

# Set up matplotlib for publication quality figures
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 12
plt.rcParams['axes.labelsize'] = 12
plt.rcParams['axes.titlesize'] = 14
plt.rcParams['xtick.labelsize'] = 10
plt.rcParams['ytick.labelsize'] = 10
plt.rcParams['legend.fontsize'] = 10
plt.rcParams['figure.titlesize'] = 16
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['savefig.bbox'] = 'tight'
plt.rcParams['axes.grid'] = True
plt.rcParams['grid.alpha'] = 0.3

# Publication-quality color palettes
CLASS1_COLORS = {
    'Class I': '#1f77b4',      # Blue
    'Class II': '#ff7f0e',     # Orange
    'Other': '#2ca02c'         # Green
}

SUPERFAMILY_COLORS = [
    '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
    '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf',
    '#aec7e8', '#ffbb78', '#98df8a', '#ff9896', '#c5b0d5',
    '#c49c94', '#f7b6d2', '#c7c7c7', '#dbdb8d', '#9edae5'
]

LENGTH_DIST_COLORS = {
    '<100': '#f7fbff',
    '100-500': '#deebf7',
    '500-1K': '#c6dbef',
    '1K-2K': '#9ecae1',
    '2K-5K': '#6baed6',
    '5K-10K': '#4292c6',
    '>10K': '#2171b5'
}


class TEVisualizer:
    """
    Visualization engine for transposable element analysis results.

    This class provides various types of plots and visualizations:
    - Pie charts for composition analysis
    - Bar charts for comparisons
    - Histograms for length distributions
    - Heatmaps for chromosome distributions
    - Multi-panel summary figures
    """

    def __init__(self, stats: TEStatistics, style: str = 'seaborn-v0_8'):
        """
        Initialize the TE visualizer.

        Args:
            stats (TEStatistics): TE statistics instance
            style (str): Matplotlib style to use
        """
        self.stats = stats
        self.style = style

        # Set matplotlib style
        try:
            if style in plt.style.available:
                plt.style.use(style)
            else:
                plt.style.use('default')
        except:
            plt.style.use('default')

        logger.info(f"Initialized TE visualizer with style: {style}")

    def _format_number(self, num: int, suffix: str = '') -> str:
        """
        Format numbers with thousand separators and suffix.

        Args:
            num (int): Number to format
            suffix (str): Suffix to add

        Returns:
            str: Formatted number string
        """
        if num >= 1_000_000:
            return f"{num/1_000_000:.1f}M{suffix}"
        elif num >= 1_000:
            return f"{num/1_000:.1f}K{suffix}"
        else:
            return f"{num}{suffix}"

    def _format_bp(self, bp: int) -> str:
        """
        Format base pairs with appropriate units.

        Args:
            bp (int): Base pairs

        Returns:
            str: Formatted base pair string
        """
        if bp >= 1_000_000_000:
            return f"{bp/1_000_000_000:.1f}Gb"
        elif bp >= 1_000_000:
            return f"{bp/1_000_000:.1f}Mb"
        elif bp >= 1_000:
            return f"{bp/1_000:.1f}Kb"
        else:
            return f"{bp}bp"

    def _get_colors(self, categories: List[str], color_type: str = 'default') -> List[str]:
        """
        Get colors for categories.

        Args:
            categories (List[str]): Category names
            color_type (str): Type of color scheme

        Returns:
            List[str]: List of colors
        """
        if color_type == 'class1':
            return [CLASS1_COLORS.get(cat, '#808080') for cat in categories]
        elif color_type == 'length_dist':
            return [LENGTH_DIST_COLORS.get(cat, '#808080') for cat in categories]
        else:
            # Use seaborn color palette
            return sns.color_palette("husl", len(categories)).as_hex()

    def plot_composition_pie(self, level: ClassLevel = ClassLevel.LEVEL1,
                           metric: str = 'merged_length',
                           save_path: Optional[str] = None,
                           figsize: Tuple[int, int] = (10, 8)) -> plt.Figure:
        """
        Create a pie chart showing TE composition at the specified level.

        Args:
            level (ClassLevel): Classification level for the pie chart
            metric (str): Metric to visualize (element_count, raw_length, merged_length)
            save_path (Optional[str]): Path to save the figure
            figsize (Tuple[int, int]): Figure size

        Returns:
            plt.Figure: Matplotlib figure object
        """
        # Get data for the specified level
        df = self.stats.get_hierarchical_summary(level)

        if df.empty or metric not in df.columns:
            raise ValueError(f"No data available for level {level} and metric {metric}")

        # Sort by the metric
        df = df.sort_values(metric, ascending=False)

        # Create figure
        fig, ax = plt.subplots(figsize=figsize)

        # Prepare data
        if level == ClassLevel.LEVEL1:
            labels = df['level1']
            colors = self._get_colors(labels, 'class1')
        elif level == ClassLevel.LEVEL2:
            labels = [f"{row['level1']}/{row['level2']}" for _, row in df.iterrows()]
            colors = self._get_colors(labels, 'default')
        else:
            labels = [f"{row['level2']}/{row['level3']}" for _, row in df.iterrows()]
            colors = self._get_colors(labels, 'default')

        sizes = df[metric]

        # Create pie chart
        wedges, texts, autotexts = ax.pie(sizes, labels=labels, colors=colors,
                                          autopct='%1.1f%%', startangle=90,
                                          textprops={'fontsize': 10})

        # Enhance text formatting
        for autotext in autotexts:
            autotext.set_color('white')
            autotext.set_fontweight('bold')

        # Title
        metric_name = metric.replace('_', ' ').title()
        level_name = level.value.replace('_', ' ').title()
        ax.set_title(f'TE Composition by {level_name}\n(Based on {metric_name})',
                    fontweight='bold', pad=20)

        # Add legend with detailed information
        legend_elements = []
        for i, (label, size) in enumerate(zip(labels, sizes)):
            if metric in ['element_count']:
                label_text = f"{label}: {self._format_number(size, '')}"
            elif 'length' in metric:
                label_text = f"{label}: {self._format_bp(size)}"
            else:
                label_text = f"{label}: {size:.1f}"

            legend_elements.append(Patch(facecolor=colors[i], label=label_text))

        ax.legend(handles=legend_elements, loc='center left', bbox_to_anchor=(1, 0.5))

        plt.tight_layout()

        # Save figure if path provided
        if save_path:
            fig.savefig(save_path, dpi=300, bbox_inches='tight')
            logger.info(f"Saved pie chart to {save_path}")

        return fig

    def plot_comparison_bar(self, level: ClassLevel = ClassLevel.LEVEL2,
                           metrics: List[str] = ['element_count', 'merged_length'],
                           save_path: Optional[str] = None,
                           figsize: Tuple[int, int] = (12, 8)) -> plt.Figure:
        """
        Create a bar chart comparing TE categories by multiple metrics.

        Args:
            level (ClassLevel): Classification level
            metrics (List[str]): Metrics to compare
            save_path (Optional[str]): Path to save the figure
            figsize (Tuple[int, int]): Figure size

        Returns:
            plt.Figure: Matplotlib figure object
        """
        # Get data
        df = self.stats.get_hierarchical_summary(level)

        if df.empty:
            raise ValueError(f"No data available for level {level}")

        # Sort by merged_length
        df = df.sort_values('merged_length', ascending=True)  # Sort ascending for horizontal bars

        # Create figure with subplots
        fig, axes = plt.subplots(1, len(metrics), figsize=figsize, sharey=True)
        if len(metrics) == 1:
            axes = [axes]

        # Prepare labels
        if level == ClassLevel.LEVEL1:
            labels = df['level1']
        elif level == ClassLevel.LEVEL2:
            labels = [f"{row['level1']}/{row['level2']}" for _, row in df.iterrows()]
        else:
            labels = [f"{row['level2']}/{row['level3']}" for _, row in df.iterrows()]

        # Create bar charts
        for i, metric in enumerate(metrics):
            if metric not in df.columns:
                logger.warning(f"Metric {metric} not available, skipping")
                continue

            values = df[metric]
            colors = self._get_colors(labels, 'default')

            # Create horizontal bar chart
            bars = axes[i].barh(range(len(labels)), values, color=colors)

            # Customize labels and title
            axes[i].set_xlabel(metric.replace('_', ' ').title())
            axes[i].set_title(metric.replace('_', ' ').title(), fontweight='bold')
            axes[i].grid(axis='x', alpha=0.3)

            # Format x-axis labels
            if metric in ['element_count']:
                axes[i].set_xscale('log')
                axes[i].xaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: self._format_number(int(x))))
            elif 'length' in metric:
                axes[i].xaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: self._format_bp(int(x))))

            # Add value labels on bars
            for j, (bar, value) in enumerate(zip(bars, values)):
                if metric in ['element_count']:
                    label = self._format_number(int(value), '')
                elif 'length' in metric:
                    label = self._format_bp(int(value))
                else:
                    label = f"{value:.1f}"

                axes[i].text(bar.get_width() * 1.01, bar.get_y() + bar.get_height()/2,
                           label, va='center', fontsize=9)

            # Set y-axis labels only for the leftmost plot
            if i == 0:
                axes[i].set_yticks(range(len(labels)))
                axes[i].set_yticklabels(labels)
            else:
                axes[i].set_yticklabels([])

        plt.tight_layout()

        # Save figure if path provided
        if save_path:
            fig.savefig(save_path, dpi=300, bbox_inches='tight')
            logger.info(f"Saved comparison bar chart to {save_path}")

        return fig

    def plot_length_distribution(self, classification_key: Optional[Tuple[str, str, str]] = None,
                                bins: int = 50,
                                save_path: Optional[str] = None,
                                figsize: Tuple[int, int] = (12, 8)) -> plt.Figure:
        """
        Plot the length distribution of TEs.

        Args:
            classification_key (Optional[Tuple[str, str, str]]): Specific TE category to analyze
            bins (int): Number of histogram bins
            save_path (Optional[str]): Path to save the figure
            figsize (Tuple[int, int]): Figure size

        Returns:
            plt.Figure: Matplotlib figure object
        """
        # Get length distribution data
        length_dist = self.stats.get_length_distribution(classification_key)

        if not length_dist:
            raise ValueError("No length distribution data available")

        # Create figure
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)

        # Prepare title
        if classification_key:
            title = f"Length Distribution: {'/'.join(classification_key)}"
        else:
            title = "Overall TE Length Distribution"

        # Get actual lengths for histogram
        if classification_key:
            records = self.stats.classified_records[classification_key]
        else:
            records = self.stats.te_records

        lengths = [record.length for record in records]

        # Histogram
        ax1.hist(lengths, bins=bins, alpha=0.7, color='steelblue', edgecolor='black')
        ax1.set_xlabel('Length (bp)')
        ax1.set_ylabel('Frequency')
        ax1.set_title('Length Distribution Histogram')
        ax1.grid(True, alpha=0.3)

        # Format x-axis for large values
        ax1.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: self._format_bp(int(x))))

        # Box plot
        ax2.boxplot(lengths, vert=True, patch_artist=True,
                   boxprops=dict(facecolor='lightblue', alpha=0.7),
                   medianprops=dict(color='red', linewidth=2))
        ax2.set_ylabel('Length (bp)')
        ax2.set_title('Length Distribution Box Plot')
        ax2.grid(True, alpha=0.3)
        ax2.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: self._format_bp(int(x))))

        # Add statistics text
        stats_text = (f"Elements: {length_dist['total_elements']:,}\n"
                     f"Mean: {length_dist['mean_length']:.0f} bp\n"
                     f"Median: {length_dist['median_length']:.0f} bp\n"
                     f"Std: {length_dist['std_length']:.0f} bp\n"
                     f"Min: {length_dist['min_length']:,} bp\n"
                     f"Max: {length_dist['max_length']:,} bp")

        ax2.text(1.1, 0.5, stats_text, transform=ax2.transAxes,
                bbox=dict(boxstyle="round,pad=0.3", facecolor="lightyellow"),
                verticalalignment='center', fontsize=10)

        plt.suptitle(title, fontweight='bold', y=1.02)
        plt.tight_layout()

        # Save figure if path provided
        if save_path:
            fig.savefig(save_path, dpi=300, bbox_inches='tight')
            logger.info(f"Saved length distribution plot to {save_path}")

        return fig

    def plot_chromosome_distribution(self, metric: str = 'merged_length',
                                   save_path: Optional[str] = None,
                                   figsize: Tuple[int, int] = (14, 8)) -> plt.Figure:
        """
        Plot TE distribution across chromosomes.

        Args:
            metric (str): Metric to visualize
            save_path (Optional[str]): Path to save the figure
            figsize (Tuple[int, int]): Figure size

        Returns:
            plt.Figure: Matplotlib figure object
        """
        # Get chromosome distribution data
        chrom_dist = self.stats.get_chromosome_distribution()

        if not chrom_dist:
            raise ValueError("No chromosome distribution data available")

        # Prepare data
        chromosomes = sorted(chrom_dist.keys())
        # Handle metric mapping for chromosome distribution
        metric_mapping = {
            'element_count': 'total_count',
            'merged_length': 'total_length',
            'raw_length': 'total_length'
        }
        actual_metric = metric_mapping.get(metric, metric)

        values = [chrom_dist[chrom][actual_metric]
                 if actual_metric in chrom_dist[chrom]
                 else 0 for chrom in chromosomes]

        # Create figure
        fig, ax = plt.subplots(figsize=figsize)

        # Create bar chart
        bars = ax.bar(range(len(chromosomes)), values, color='steelblue', alpha=0.7, edgecolor='black')

        # Customize plot
        ax.set_xlabel('Chromosome')
        ax.set_ylabel(metric.replace('_', ' ').title())
        ax.set_title(f'TE Distribution Across Chromosomes ({metric.replace("_", " ").title()})',
                    fontweight='bold')
        ax.set_xticks(range(len(chromosomes)))
        ax.set_xticklabels(chromosomes, rotation=45, ha='right')
        ax.grid(True, alpha=0.3)

        # Format y-axis
        if metric in ['element_count', 'total_count']:
            ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: self._format_number(int(x))))
        elif 'length' in metric:
            ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: self._format_bp(int(x))))

        # Add value labels on bars
        for i, (bar, value) in enumerate(zip(bars, values)):
            if metric in ['element_count', 'total_count']:
                label = self._format_number(int(value), '')
            elif 'length' in metric:
                label = self._format_bp(int(value))
            else:
                label = f"{value:.1f}"

            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() * 1.01,
                   label, ha='center', va='bottom', fontsize=9)

        plt.tight_layout()

        # Save figure if path provided
        if save_path:
            fig.savefig(save_path, dpi=300, bbox_inches='tight')
            logger.info(f"Saved chromosome distribution plot to {save_path}")

        return fig

    def plot_summary_dashboard(self, save_path: Optional[str] = None,
                             figsize: Tuple[int, int] = (16, 12)) -> plt.Figure:
        """
        Create a comprehensive summary dashboard with multiple plots.

        Args:
            save_path (Optional[str]): Path to save the figure
            figsize (Tuple[int, int]): Figure size

        Returns:
            plt.Figure: Matplotlib figure object
        """
        # Create figure with subplots
        fig = plt.figure(figsize=figsize)

        # Define subplot layout
        gs = fig.add_gridspec(3, 3, hspace=0.3, wspace=0.3)

        # 1. Class Level Pie Chart (top left)
        ax1 = fig.add_subplot(gs[0, 0])
        df_level1 = self.stats.get_hierarchical_summary(ClassLevel.LEVEL1)
        if not df_level1.empty:
            df_level1 = df_level1.sort_values('merged_length', ascending=False)
            colors = self._get_colors(df_level1['level1'], 'class1')
            wedges, texts, autotexts = ax1.pie(df_level1['merged_length'],
                                              labels=df_level1['level1'],
                                              colors=colors, autopct='%1.1f%%')
            ax1.set_title('Class Level Composition', fontweight='bold')
            for autotext in autotexts:
                autotext.set_color('white')
                autotext.set_fontweight('bold')

        # 2. Superfamily Level Bar Chart (top middle and right)
        ax2 = fig.add_subplot(gs[0, 1:])
        df_level2 = self.stats.get_hierarchical_summary(ClassLevel.LEVEL2)
        if not df_level2.empty:
            df_level2 = df_level2.sort_values('merged_length', ascending=True)
            labels = [f"{row['level1']}/{row['level2']}" for _, row in df_level2.iterrows()]
            bars = ax2.barh(range(len(labels)), df_level2['merged_length'],
                           color= sns.color_palette("husl", len(labels)))
            ax2.set_yticks(range(len(labels)))
            ax2.set_yticklabels(labels, fontsize=10)
            ax2.set_xlabel('Merged Length (bp)')
            ax2.set_title('Superfamily Level Coverage', fontweight='bold')
            ax2.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: self._format_bp(int(x))))

        # 3. Length Distribution Histogram (middle left)
        ax3 = fig.add_subplot(gs[1, 0])
        length_dist = self.stats.get_length_distribution()
        if length_dist:
            lengths = [record.length for record in self.stats.te_records]
            ax3.hist(lengths, bins=30, alpha=0.7, color='lightcoral', edgecolor='black')
            ax3.set_xlabel('Length (bp)')
            ax3.set_ylabel('Frequency')
            ax3.set_title('Overall Length Distribution', fontweight='bold')
            ax3.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: self._format_bp(int(x))))

        # 4. Chromosome Distribution (middle middle and right)
        ax4 = fig.add_subplot(gs[1, 1:])
        chrom_dist = self.stats.get_chromosome_distribution()
        if chrom_dist:
            chromosomes = sorted(chrom_dist.keys())
            values = [chrom_dist[chrom]['total_length'] for chrom in chromosomes]
            bars = ax4.bar(range(len(chromosomes)), values, color='lightgreen', alpha=0.7)
            ax4.set_xlabel('Chromosome')
            ax4.set_ylabel('Total TE Length (bp)')
            ax4.set_title('Chromosome Distribution', fontweight='bold')
            ax4.set_xticks(range(len(chromosomes)))
            ax4.set_xticklabels(chromosomes, rotation=45, ha='right')
            ax4.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: self._format_bp(int(x))))

        # 5. Top Families Table (bottom)
        ax5 = fig.add_subplot(gs[2, :])
        ax5.axis('off')  # Hide axes for table

        # Get top families
        top_families = self.stats.get_top_categories('merged_length', ClassLevel.LEVEL3, 10)
        if not top_families.empty:
            # Create table data
            table_data = []
            for _, row in top_families.iterrows():
                family_name = f"{row['level2']}/{row['level3']}"
                table_data.append([
                    family_name[:30],  # Limit family name length
                    f"{row['element_count']:,}",
                    self._format_bp(row['merged_length']),
                    f"{row['length_mean']:.0f} bp",
                    f"{row['genome_percentage']:.2f}%"
                ])

            # Create table
            table = ax5.table(cellText=table_data,
                            colLabels=['Family', 'Count', 'Length', 'Mean Len', 'Genome %'],
                            cellLoc='left',
                            loc='center',
                            colWidths=[0.35, 0.15, 0.2, 0.15, 0.15])

            table.auto_set_font_size(False)
            table.set_fontsize(10)
            table.scale(1, 1.5)

            # Style the table
            for i in range(len(table_data) + 1):  # +1 for header
                for j in range(5):
                    cell = table[i, j]
                    if i == 0:  # Header row
                        cell.set_facecolor('#4CAF50')
                        cell.set_text_props(weight='bold', color='white')
                    else:
                        cell.set_facecolor('#f0f0f0' if i % 2 == 0 else 'white')

            ax5.set_title('Top 10 TE Families by Coverage', fontweight='bold', pad=20)

        # Main title
        fig.suptitle('Transposable Element Analysis Dashboard', fontsize=18, fontweight='bold', y=0.95)

        # Add summary statistics as text
        total_elements = len(self.stats.te_records)
        total_length = sum(summary.merged_length for summary in self.stats.summary_stats.values())
        unique_families = len(self.stats.summary_stats)

        summary_text = (f"Total Elements: {self._format_number(total_elements)} | "
                       f"Total Length: {self._format_bp(total_length)} | "
                       f"Unique Families: {unique_families}")

        fig.text(0.5, 0.02, summary_text, ha='center', fontsize=12,
                bbox=dict(boxstyle="round,pad=0.5", facecolor="lightblue", alpha=0.8))

        plt.tight_layout(rect=[0, 0.05, 1, 0.93])  # Adjust for title and summary text

        # Save figure if path provided
        if save_path:
            fig.savefig(save_path, dpi=300, bbox_inches='tight')
            logger.info(f"Saved summary dashboard to {save_path}")

        return fig

    def create_all_plots(self, output_dir: str, formats: List[str] = ['png', 'svg']) -> None:
        """
        Generate all standard plots and save them to the output directory.

        Args:
            output_dir (str): Output directory path
            formats (List[str]): Output formats (png, svg, pdf)
        """
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        logger.info(f"Generating all plots in {output_dir}")

        # Generate all plots
        plot_functions = [
            (self.plot_composition_pie, [ClassLevel.LEVEL1], 'composition_class1'),
            (self.plot_composition_pie, [ClassLevel.LEVEL2], 'composition_class2'),
            (self.plot_comparison_bar, [ClassLevel.LEVEL2], 'comparison_class2'),
            (self.plot_length_distribution, [], 'length_distribution'),
            (self.plot_chromosome_distribution, [], 'chromosome_distribution'),
            (self.plot_summary_dashboard, [], 'summary_dashboard')
        ]

        for plot_func, args, plot_name in plot_functions:
            try:
                fig = plot_func(*args)

                # Save in all requested formats
                for fmt in formats:
                    file_path = output_path / f"{plot_name}.{fmt}"
                    fig.savefig(file_path, dpi=300, bbox_inches='tight')
                    logger.info(f"Saved {plot_name} as {fmt}")

                plt.close(fig)  # Close figure to free memory

            except Exception as e:
                logger.error(f"Error generating {plot_name}: {e}")

        logger.info("Plot generation completed")


def main():
    """Example usage of TEVisualizer."""
    from .gff_parser import GFFParser
    from .classifier import TEClassifier
    from .te_statistics import TEStatistics

    # Example with test data
    parser = GFFParser("../test.gff")
    records = parser.parse_gff_file()

    classifier = TEClassifier()
    stats = TEStatistics(records, classifier, genome_size=16_000_000)
    stats.calculate_all_statistics()

    # Create visualizer
    visualizer = TEVisualizer(stats)

    # Generate all plots
    visualizer.create_all_plots("plots", ['png', 'svg'])


if __name__ == "__main__":
    main()