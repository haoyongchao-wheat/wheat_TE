"""
TE Statistics Package
A comprehensive tool for transposable element statistics and analysis.
"""

__version__ = "1.0.0"
__author__ = "Bioinformatics Team"

from .gff_parser import GFFParser
from .classifier import TEClassifier
from .te_statistics import TEStatistics
from .visualizer import TEVisualizer

__all__ = ['GFFParser', 'TEClassifier', 'TEStatistics', 'TEVisualizer']