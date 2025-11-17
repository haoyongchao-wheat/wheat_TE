#!/usr/bin/env python3
"""
TE Classifier Module for Transposable Element Hierarchical Classification
This module provides a hierarchical classification system for transposable elements.
"""

import re
from typing import Dict, List, Tuple, Optional, Set
from enum import Enum
import logging
from dataclasses import dataclass

logger = logging.getLogger(__name__)


class ClassLevel(Enum):
    """Enumeration for classification levels."""
    LEVEL1 = "Class Level"      # Class I, Class II, Other
    LEVEL2 = "Superfamily"      # LTR, LINE, SINE, TIR, etc.
    LEVEL3 = "Family"           # Specific families like LTR/Copia, LINE/L1, etc.


@dataclass
class TECategory:
    """Data class representing a TE category with its hierarchical information."""
    name: str
    level: ClassLevel
    parent: Optional[str] = None
    children: Optional[List[str]] = None
    description: Optional[str] = None

    def __post_init__(self):
        if self.children is None:
            self.children = []


class TEClassifier:
    """
    Hierarchical classifier for transposable elements.

    This class implements a three-level classification system:
    Level 1: Class I, Class II, Other
    Level 2: Superfamily (LTR, LINE, SINE, TIR, etc.)
    Level 3: Specific families
    """

    def __init__(self, custom_classification: Optional[Dict] = None):
        """
        Initialize the TE classifier.

        Args:
            custom_classification (Optional[Dict]): Custom classification rules
        """
        self.classification_tree = self._build_classification_tree(custom_classification)
        self.class_lookup = self._build_class_lookup()
        self.reverse_lookup = self._build_reverse_lookup()

        logger.info(f"Initialized TE classifier with {len(self.class_lookup)} classification rules")

    def _build_classification_tree(self, custom_classification: Optional[Dict] = None) -> Dict:
        """
        Build the hierarchical classification tree.

        Args:
            custom_classification (Optional[Dict]): Custom classification rules to override defaults

        Returns:
            Dict: Classification tree structure
        """
        # Default classification based on standard TE classification
        default_classification = {
            "Class I": {
                "description": "Retrotransposons - RNA intermediate transposons",
                "children": {
                    "LTR": {
                        "description": "Long Terminal Repeat retrotransposons",
                        "children": {
                            "LTR/Copia": "Copia-like LTR retrotransposons",
                            "LTR/Gypsy": "Gypsy-like LTR retrotransposons",
                            "LTR/ERV1": "Endogenous retrovirus-like LTR",
                            "LTR/ERV2": "Endogenous retrovirus-like LTR type 2",
                            "LTR/Other": "Other LTR retrotransposons"
                        }
                    },
                    "LINE": {
                        "description": "Long Interspersed Nuclear Elements",
                        "children": {
                            "LINE/L1": "LINE1-like elements",
                            "LINE/L2": "LINE2-like elements",
                            "LINE/Jockey": "Jockey-like elements",
                            "LINE/RTE": "RTE-like elements",
                            "LINE/I": "LINE I-like elements",
                            "LINE/Other": "Other LINE elements"
                        }
                    },
                    "SINE": {
                        "description": "Short Interspersed Nuclear Elements",
                        "children": {
                            "SINE/Alu": "Alu-like SINE elements",
                            "SINE/tRNA": "tRNA-derived SINE elements",
                            "SINE/Other": "Other SINE elements"
                        }
                    },
                    "DIRS": {
                        "description": "DIRS-like elements",
                        "children": {
                            "DIRS/Other": "Other DIRS-like elements"
                        }
                    }
                }
            },
            "Class II": {
                "description": "DNA transposons - DNA-to-DNA transposons",
                "children": {
                    "TIR": {
                        "description": "Terminal Inverted Repeat transposons",
                        "children": {
                            "DNA/TcMar": "Tc1/Mariner DNA transposons",
                            "DNA/hAT": "hAT family DNA transposons",
                            "DNA/MuDR": "Mutator family DNA transposons",
                            "DNA/CACTA": "CACTA family DNA transposons",
                            "DNA/EnSpm": "EnSpm/CACTA family DNA transposons",
                            "DNA/Pif": "Pif family DNA transposons",
                            "DNA/PIF-harbinger": "PIF-Harbinger DNA transposons",
                            "DNA/Tc1-Mariner": "Tc1-Mariner DNA transposons",
                            "DNA/Merlin": "Merlin-like DNA transposons",
                            "DNA/TIR-Other": "Other TIR DNA transposons"
                        }
                    },
                    "Helitron": {
                        "description": "Rolling-circle DNA transposons",
                        "children": {
                            "DNA/Helitron": "Helitron DNA transposons",
                            "DNA/Helitron-Other": "Other Helitron-like elements"
                        }
                    },
                    "Maverick": {
                        "description": "Maverick/Polinton DNA transposons",
                        "children": {
                            "DNA/Maverick": "Maverick/Polinton DNA transposons",
                            "DNA/Maverick-Other": "Other Maverick-like elements"
                        }
                    },
                    "DNA": {
                        "description": "Other DNA transposons",
                        "children": {
                            "DNA/DNA": "Generic DNA transposons",
                            "DNA/Other": "Other DNA transposons",
                            "DNA/Unspecified": "Unspecified DNA transposons",
                            "DNA/Ginger": "Ginger family DNA transposons",
                            "DNA/MULE": "Mutator-like elements (MULEs)"
                        }
                    }
                }
            },
            "Other": {
                "description": "Unclassified or unknown elements",
                "children": {
                    "Unspecified": {
                        "description": "Unspecified transposons",
                        "children": {
                            "Unspecified": "Generic unspecified elements"
                        }
                    },
                    "Unknown": {
                        "description": "Unknown transposons",
                        "children": {
                            "Unknown": "Generic unknown elements"
                        }
                    },
                    "Satellite": {
                        "description": "Satellite DNA",
                        "children": {
                            "Satellite": "Satellite repeats"
                        }
                    },
                    "Simple_repeat": {
                        "description": "Simple repeats",
                        "children": {
                            "Simple_repeat": "Simple sequence repeats"
                        }
                    },
                    "rRNA": {
                        "description": "Ribosomal RNA repeats",
                        "children": {
                            "rRNA": "rRNA gene repeats"
                        }
                    },
                    "snRNA": {
                        "description": "Small nuclear RNA repeats",
                        "children": {
                            "snRNA": "snRNA gene repeats"
                        }
                    },
                    "tRNA": {
                        "description": "Transfer RNA repeats",
                        "children": {
                            "tRNA": "tRNA gene repeats"
                        }
                    },
                    "Other": {
                        "description": "Other repetitive elements",
                        "children": {
                            "Other": "Other repetitive elements"
                        }
                    }
                }
            }
        }

        # Override with custom classification if provided
        if custom_classification:
            classification = custom_classification
        else:
            classification = default_classification

        return classification

    def _build_class_lookup(self) -> Dict[str, Tuple[str, str, str]]:
        """
        Build a lookup table for quick classification.

        Returns:
            Dict[str, Tuple[str, str, str]]: Maps class names to (level1, level2, level3) hierarchy
        """
        lookup = {}

        def process_tree(tree: Dict, level1: str = "", level2: str = "") -> None:
            """
            Recursively process the classification tree to build lookup.

            Args:
                tree (Dict): Current tree level
                level1 (str): Current level1 category
                level2 (str): Current level2 category
            """
            for key, value in tree.items():
                if key == "description":
                    continue

                if key == "children":
                    # Process children
                    for child_key, child_value in value.items():
                        if isinstance(child_value, dict):
                            if level2:
                                # This is level3
                                lookup[child_key] = (level1, level2, child_key)
                                process_tree(child_value, level1, level2)
                            elif level1:
                                # This is level2
                                process_tree(child_value, level1, child_key)
                        else:
                            # Leaf node with description
                            lookup[child_key] = (level1, level2, child_key)
                elif isinstance(value, dict) and "children" in value:
                    # This is a category with children
                    if not level1:
                        # This is level1
                        process_tree(value, key, "")
                    elif not level2:
                        # This is level2
                        process_tree(value, level1, key)

        process_tree(self.classification_tree)

        # Add common patterns and variations
        variations = {
            # Class I variations
            "LTR": ("Class I", "LTR", "LTR"),
            "LINE": ("Class I", "LINE", "LINE"),
            "SINE": ("Class I", "SINE", "SINE"),
            "retrotransposon": ("Class I", "Other", "Other"),
            "Retrotransposon": ("Class I", "Other", "Other"),

            # Class II variations
            "DNA": ("Class II", "DNA", "DNA"),
            "DNA/DNA": ("Class II", "DNA", "DNA"),
            "DNA_transposon": ("Class II", "Other", "Other"),
            "DNA_transposase": ("Class II", "Other", "Other"),
            "Helitron": ("Class II", "Helitron", "Helitron"),
            "TIR": ("Class II", "TIR", "TIR"),

            # Handle common patterns with slashes
            # Will be processed by the pattern matching in classify_element
        }

        lookup.update(variations)

        return lookup

    def _build_reverse_lookup(self) -> Dict[str, Set[str]]:
        """
        Build a reverse lookup for pattern matching.

        Returns:
            Dict[str, Set[str]]: Maps keywords to possible categories
        """
        reverse_lookup = {}

        for class_name, (level1, level2, level3) in self.class_lookup.items():
            # Add the full class name
            if class_name not in reverse_lookup:
                reverse_lookup[class_name] = set()
            reverse_lookup[class_name].add(class_name)

            # Add parts of the name
            parts = class_name.split('/')
            for part in parts:
                if part not in reverse_lookup:
                    reverse_lookup[part] = set()
                reverse_lookup[part].add(class_name)

            # Add common variations and substrings
            if 'copia' in class_name.lower():
                if 'copia' not in reverse_lookup:
                    reverse_lookup['copia'] = set()
                reverse_lookup['copia'].add(class_name)
            if 'gypsy' in class_name.lower():
                if 'gypsy' not in reverse_lookup:
                    reverse_lookup['gypsy'] = set()
                reverse_lookup['gypsy'].add(class_name)
            if 'line' in class_name.lower():
                if 'line' not in reverse_lookup:
                    reverse_lookup['line'] = set()
                reverse_lookup['line'].add(class_name)
            if 'sine' in class_name.lower():
                if 'sine' not in reverse_lookup:
                    reverse_lookup['sine'] = set()
                reverse_lookup['sine'].add(class_name)

        return reverse_lookup

    def classify_element(self, class_name: str) -> Tuple[str, str, str]:
        """
        Classify a TE element into the three-level hierarchy.

        Args:
            class_name (str): Original class name from GFF

        Returns:
            Tuple[str, str, str]: (level1, level2, level3) classification
        """
        if not class_name or class_name == '.':
            return ("Other", "Unspecified", "Unspecified")

        # Clean up the class name
        clean_class = class_name.strip()

        # Direct lookup first
        if clean_class in self.class_lookup:
            return self.class_lookup[clean_class]

        # Try exact match with different case
        for lookup_key in self.class_lookup:
            if clean_class.lower() == lookup_key.lower():
                return self.class_lookup[lookup_key]

        # Pattern matching for entries with slashes
        if '/' in clean_class:
            parts = clean_class.split('/')

            # Try to match patterns like "LTR/Copia", "LINE/L1", etc.
            for lookup_key in self.class_lookup:
                lookup_parts = lookup_key.split('/')

                # Check if all parts match (in order)
                if len(parts) <= len(lookup_parts):
                    match = True
                    for i, part in enumerate(parts):
                        if i < len(lookup_parts):
                            if part.lower() != lookup_parts[i].lower():
                                match = False
                                break

                    if match:
                        return self.class_lookup[lookup_key]

        # Fuzzy matching with partial matches
        clean_lower = clean_class.lower()

        # Check for key substrings
        for key, categories in self.reverse_lookup.items():
            if key.lower() in clean_lower or clean_lower in key.lower():
                # Return the first match (categories is a set, so we need to get an element)
                for category in categories:
                    if category in self.class_lookup:
                        return self.class_lookup[category]

        # Default classification based on keywords (more specific matching to avoid false positives)
        if any(keyword in clean_lower for keyword in [' ltr ', 'retrotransposon', ' copia ', ' gypsy ']):
            return ("Class I", "LTR", "LTR")
        elif any(keyword in clean_lower for keyword in [' line ', ' line1', ' line2', ' jockey ', 'rte']):
            return ("Class I", "LINE", "LINE")
        elif any(keyword in clean_lower for keyword in [' sine ', ' alu ', ' 7sl ', ' 5s ']):
            return ("Class I", "SINE", "SINE")
        elif any(keyword in clean_lower for keyword in [' tir ', 'tc1', 'mariner', ' hat ', 'helifield']):
            return ("Class II", "TIR", "TIR")
        elif any(keyword in clean_lower for keyword in [' helitron ', 'rolling_circle']):
            return ("Class II", "Helitron", "Helitron")
        elif any(keyword in clean_lower for keyword in ['maverick', 'polinton']):
            return ("Class II", "Maverick", "Maverick")
        elif any(keyword in clean_lower for keyword in ['dna_transposon', 'dna/transposon']):
            return ("Class II", "DNA", "DNA")
        elif any(keyword in clean_lower for keyword in ['unspecified', 'unknown']):
            return ("Other", "Unspecified", "Unspecified")
        elif any(keyword in clean_lower for keyword in ['satellite', 'repeat', 'rrna', 'snrna', 'trna']):
            return ("Other", "Other", "Other")
        else:
            logger.warning(f"Unable to classify element: {class_name}, defaulting to Other/Unspecified")
            return ("Other", "Unspecified", "Unspecified")

    def get_classification_summary(self) -> Dict:
        """
        Get a summary of the classification system.

        Returns:
            Dict: Summary of classification levels and categories
        """
        level1_categories = set()
        level2_categories = set()
        level3_categories = set()

        for level1, level2, level3 in self.class_lookup.values():
            level1_categories.add(level1)
            level2_categories.add(level2)
            level3_categories.add(level3)

        return {
            'num_level1_categories': len(level1_categories),
            'num_level2_categories': len(level2_categories),
            'num_level3_categories': len(level3_categories),
            'level1_categories': sorted(list(level1_categories)),
            'level2_categories': sorted(list(level2_categories)),
            'total_classification_rules': len(self.class_lookup)
        }

    def get_category_children(self, category: str, level: ClassLevel) -> List[str]:
        """
        Get children categories for a given category and level.

        Args:
            category (str): Parent category name
            level (ClassLevel): The level of the parent category

        Returns:
            List[str]: List of child category names
        """
        children = []

        def find_children(tree: Dict, target_category: str, current_level: int) -> None:
            if current_level >= level.value:
                return

            for key, value in tree.items():
                if key == target_category and "children" in value:
                    if isinstance(value["children"], dict):
                        children.extend(value["children"].keys())
                elif isinstance(value, dict):
                    find_children(value, target_category, current_level + 1)

        find_children(self.classification_tree, category, 0)
        return sorted(children)

    def export_classification_rules(self, output_file: str) -> None:
        """
        Export classification rules to a CSV file.

        Args:
            output_file (str): Path to output CSV file
        """
        with open(output_file, 'w') as f:
            f.write("Original_Class,Level1,Level2,Level3\n")

            for original_class, (level1, level2, level3) in sorted(self.class_lookup.items()):
                f.write(f"{original_class},{level1},{level2},{level3}\n")

        logger.info(f"Exported {len(self.class_lookup)} classification rules to {output_file}")


def main():
    """Example usage of TEClassifier."""
    classifier = TEClassifier()

    # Test classification
    test_classes = [
        "LTR/Copia",
        "LINE/L1",
        "DNA/TcMar",
        "DNA/DNA",
        "Unspecified",
        "Unknown",
        "Helitron",
        "SINE/Alu"
    ]

    print("=== TE Classification Test ===")
    for test_class in test_classes:
        level1, level2, level3 = classifier.classify_element(test_class)
        print(f"{test_class:20} -> {level1:10} | {level2:15} | {level3}")

    # Print classification summary
    summary = classifier.get_classification_summary()
    print(f"\n=== Classification System Summary ===")
    print(f"Level 1 categories: {summary['num_level1_categories']}")
    print(f"Level 2 categories: {summary['num_level2_categories']}")
    print(f"Level 3 categories: {summary['num_level3_categories']}")
    print(f"Total rules: {summary['total_classification_rules']}")


if __name__ == "__main__":
    main()