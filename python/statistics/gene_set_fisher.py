# -*- coding: utf-8 -*-
# Script: gene_set_fisher.py

"""
Gene Set Fisher's Exact Test

This script performs Fisher's exact test to determine if a set of genes is enriched
in a pathway or gene set compared to a background set.

Usage:
    python gene_set_fisher.py [--input <input_file>] [--output <output_file>]

Arguments:
    --input: Path to input file (Excel or CSV) containing gene lists
    --output: Path to save the output results
    --format: Input file format (excel or csv, default: excel)

Environment Variables (can be set in .env file):
    INPUT_FILE: Default path to input file
    OUTPUT_FILE: Default path to output file
    RESULTS_DIR: Base directory for results (used if specific paths not provided)
"""

import os
import pandas as pd
import numpy as np
from scipy.stats import fisher_exact
import argparse
import logging
from dotenv import load_dotenv

# Load environment variables from .env file if it exists
load_dotenv()

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def show_input_template():
    """Show a template for the input file format."""
    template = """
# Input file template (Excel or CSV)
# Required columns: background, up-regulated, down-regulated, pathway
# Each column should contain gene IDs

background    | up-regulated | down-regulated | pathway
----------------|--------------|----------------|------------
Gene1          | Gene1        | Gene5          | Gene1
Gene2          | Gene3        | Gene7          | Gene2
Gene3          | Gene4        | Gene9          | Gene4
...            | ...          | ...            | ...
"""
    logger.info("Expected input file format:")
    for line in template.strip().split('\n'):
        logger.info(line)

def load_gene_data(file_path, file_format='excel'):
    """Load gene data from an Excel or CSV file."""
    try:
        if file_format.lower() == 'excel':
            df = pd.read_excel(file_path)
        elif file_format.lower() == 'csv':
            df = pd.read_csv(file_path)
        else:
            logger.error(f"Unsupported file format: {file_format}. Use 'excel' or 'csv'.")
            show_input_template()
            raise ValueError(f"Unsupported file format: {file_format}")
            
        logger.info(f"Loaded gene data from {file_path}: {df.shape[0]} rows, {df.shape[1]} columns")
        
        # Show actual columns found in the file
        logger.info(f"Columns found in the file: {', '.join(df.columns.tolist())}")
        
        # Check if required columns exist
        required_cols = ['background', 'up-regulated', 'down-regulated', 'pathway']
        missing_cols = [col for col in required_cols if col not in df.columns]
        
        # Check for case-insensitive matches or common variations
        if missing_cols:
            # Create a mapping of lowercase column names to actual column names
            col_map = {col.lower().replace(' ', '-').replace('_', '-'): col for col in df.columns}
            
            # Try to find case-insensitive matches for missing columns
            fixed_cols = {}
            still_missing = []
            for col in missing_cols:
                col_key = col.lower()
                if col_key in col_map:
                    fixed_cols[col] = col_map[col_key]
                # Try common variations
                elif col_key == 'up-regulated' and 'up' in col_map:
                    fixed_cols[col] = col_map['up']
                elif col_key == 'down-regulated' and 'down' in col_map:
                    fixed_cols[col] = col_map['down']
                elif col_key == 'background' and 'all' in col_map:
                    fixed_cols[col] = col_map['all']
                elif col_key == 'pathway' and 'gene_set' in col_map:
                    fixed_cols[col] = col_map['gene_set']
                else:
                    still_missing.append(col)
            
            if fixed_cols:
                logger.info(f"Found potential column matches: {', '.join([f'{k} -> {v}' for k, v in fixed_cols.items()])}")
                logger.info("Using these columns as substitutes.")
                
                # Use the matched columns
                for req_col, actual_col in fixed_cols.items():
                    if req_col == 'background':
                        background_genes = set(df[actual_col].dropna())
                    elif req_col == 'up-regulated':
                        upregulated_genes = set(df[actual_col].dropna())
                    elif req_col == 'down-regulated':
                        downregulated_genes = set(df[actual_col].dropna())
                    elif req_col == 'pathway':
                        pathway_genes = set(df[actual_col].dropna())
                
                # If we still have missing columns, raise an error
                if still_missing:
                    logger.error(f"Required columns still not found: {', '.join(still_missing)}")
                    show_input_template()
                    raise ValueError(f"Required columns not found in input file: {', '.join(still_missing)}")
            else:
                logger.error(f"Required columns not found in input file: {', '.join(missing_cols)}")
                show_input_template()
                raise ValueError(f"Required columns not found in input file: {', '.join(missing_cols)}")
        else:
            # Convert columns to sets, removing NaN values
            background_genes = set(df['background'].dropna())
            upregulated_genes = set(df['up-regulated'].dropna())
            downregulated_genes = set(df['down-regulated'].dropna())
            pathway_genes = set(df['pathway'].dropna())
        
        logger.info(f"Loaded {len(background_genes)} background genes")
        logger.info(f"Loaded {len(upregulated_genes)} upregulated genes")
        logger.info(f"Loaded {len(downregulated_genes)} downregulated genes")
        logger.info(f"Loaded {len(pathway_genes)} pathway genes")
        
        return background_genes, upregulated_genes, downregulated_genes, pathway_genes
    
    except Exception as e:
        logger.error(f"Error loading gene data from {file_path}: {e}")
        show_input_template()
        raise

def pathway_enrichment_test(diff_genes, pathway_genes, background_genes):
    """Perform Fisher's exact test for pathway enrichment.
    
    Args:
        diff_genes: Set of differentially expressed genes (up or down)
        pathway_genes: Set of genes in the pathway
        background_genes: Set of all genes (universe)
        
    Returns:
        dict: Dictionary containing test results and statistics
    """
    # Only consider pathway genes that are in background
    pathway_in_background = pathway_genes.intersection(background_genes)
    
    # Create 2x2 contingency table
    a = len(diff_genes.intersection(pathway_in_background))  # DE genes in pathway
    b = len(diff_genes) - a  # DE genes not in pathway
    c = len(pathway_in_background) - a  # Non-DE genes in pathway
    d = len(background_genes) - a - b - c  # Non-DE genes not in pathway
    
    # Create contingency matrix
    contingency_table = np.array([[a, b], [c, d]])
    
    # Perform Fisher's exact test
    odds_ratio, p_value = fisher_exact(contingency_table, alternative='greater')
    
    # Get the list of overlapping genes
    overlap_genes = diff_genes.intersection(pathway_in_background)
    
    return {
        'genes_in_pathway': a,
        'total_diff_genes': len(diff_genes),
        'pathway_size': len(pathway_in_background),
        'background_size': len(background_genes),
        'odds_ratio': odds_ratio,
        'p_value': p_value,
        'contingency_table': contingency_table,
        'overlap_genes': overlap_genes
    }

def print_help():
    """Print a more detailed help message with examples."""
    help_text = """
==============================================================================
GENE SET FISHER'S EXACT TEST
==============================================================================

This script performs Fisher's exact test to determine if a set of genes is enriched
in a pathway or gene set compared to a background set.

USAGE:
    python gene_set_fisher.py [--input <file>] [--output <file>] [--format <format>]

ARGUMENTS:
    --input       Path to input file (Excel or CSV) containing gene lists
    --output      Path to save the output results (CSV)
    --format      Input file format: 'excel' or 'csv' (default: excel)
    --help, -h    Show this help message and exit

ENVIRONMENT VARIABLES (can be set in .env file):
    INPUT_FILE    Default path to input file
    OUTPUT_FILE   Default path to output file
    RESULTS_DIR   Base directory for results (used if specific paths not provided)

EXAMPLES:
    # Run with command-line arguments
    python gene_set_fisher.py --input gene_lists.xlsx --output results.csv

    # Run with CSV input
    python gene_set_fisher.py --input gene_lists.csv --format csv --output results.csv

    # Run with environment variables (set in .env file)
    python gene_set_fisher.py

INPUT FILE FORMAT:
"""
    print(help_text)
    
    # Show file format template
    template = """
# Input file template (Excel or CSV)
# Required columns: background, up-regulated, down-regulated, pathway
# Each column should contain gene IDs

background    | up-regulated | down-regulated | pathway
----------------|--------------|----------------|------------
Gene1          | Gene1        | Gene5          | Gene1
Gene2          | Gene3        | Gene7          | Gene2
Gene3          | Gene4        | Gene9          | Gene4
...            | ...          | ...            | ...
"""
    print(template)
    print("\n==============================================================================\n")

def save_results(up_result, down_result, output_file):
    """Save the results to a CSV file."""
    # Create a DataFrame for the results
    results = pd.DataFrame([
        {
            'gene_set': 'up-regulated',
            'genes_in_pathway': up_result['genes_in_pathway'],
            'total_genes': up_result['total_diff_genes'],
            'pathway_size': up_result['pathway_size'],
            'background_size': up_result['background_size'],
            'odds_ratio': up_result['odds_ratio'],
            'p_value': up_result['p_value'],
            'overlap_genes': ','.join(up_result['overlap_genes']) if up_result['overlap_genes'] else ''
        },
        {
            'gene_set': 'down-regulated',
            'genes_in_pathway': down_result['genes_in_pathway'],
            'total_genes': down_result['total_diff_genes'],
            'pathway_size': down_result['pathway_size'],
            'background_size': down_result['background_size'],
            'odds_ratio': down_result['odds_ratio'],
            'p_value': down_result['p_value'],
            'overlap_genes': ','.join(down_result['overlap_genes']) if down_result['overlap_genes'] else ''
        }
    ])
    
    # Create output directory if it doesn't exist
    output_dir = os.path.dirname(output_file)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
        logger.info(f"Created output directory: {output_dir}")
    
    # Save to CSV
    results.to_csv(output_file, index=False)
    logger.info(f"Results saved to {output_file}")
    
    return results

def main():
    """Main function to run the analysis."""
    
    # Get default paths from environment variables
    results_dir = os.getenv('RESULTS_DIR')
    default_input_file = os.getenv('INPUT_FILE')
    default_output_file = os.getenv('OUTPUT_FILE')
    
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Gene Set Fisher\'s Exact Test', add_help=False)
    parser.add_argument('--input', help='Path to input file (Excel or CSV)')
    parser.add_argument('--output', help='Path to save output CSV file')
    parser.add_argument('--format', default='excel', help='Input file format: excel or csv')
    parser.add_argument('--help', '-h', action='store_true', help='Show detailed help message')
    args = parser.parse_args()
    
    # Show help and exit if requested
    if args.help:
        print_help()
        return
    
    # Command-line arguments take precedence over environment variables
    input_file = args.input or default_input_file
    output_file = args.output or default_output_file
    file_format = args.format
    
    # If results_dir is set and relative paths are provided, make them absolute
    if results_dir:
        if input_file and not os.path.isabs(input_file):
            input_file = os.path.join(results_dir, input_file)
        if output_file and not os.path.isabs(output_file):
            output_file = os.path.join(results_dir, output_file)
    
    # Check if required files are specified
    if not input_file:
        print("\n" + "=" * 80)
        print("ERROR: INPUT FILE PATH NOT SPECIFIED")
        print("=" * 80)
        print("\nPlease provide an input file path using one of these methods:")
        print("  1. Command-line argument: --input <file>")
        print("  2. Environment variable in a .env file: INPUT_FILE=<file>")
        print("\nFor more information, run:")
        print("  python gene_set_fisher.py --help")
        print("=" * 80)
        logger.error("Input file path not specified")
        return 1
    
    if not output_file:
        # Generate default output file name based on input file
        output_file = os.path.splitext(input_file)[0] + "_results.csv"
        logger.info(f"Output file not specified, using: {output_file}")
    
    # Check if input file exists
    if not os.path.exists(input_file):
        print("\n" + "=" * 80)
        print("ERROR: INPUT FILE NOT FOUND")
        print("=" * 80)
        print(f"\nThe input file could not be found: {input_file}")
        print("\nPlease check that the file path is correct and the file exists.")
        print("\nExpected input file format:")
        show_input_template()
        print("=" * 80)
        logger.error(f"Input file not found: {input_file}")
        return 1
    
    try:
        # Load gene data
        logger.info("Loading gene data...")
        background_genes, upregulated_genes, downregulated_genes, pathway_genes = \
            load_gene_data(input_file, file_format)
        
        # Perform Fisher's exact test for upregulated genes
        logger.info("Performing Fisher's exact test for upregulated genes...")
        up_result = pathway_enrichment_test(upregulated_genes, pathway_genes, background_genes)
        
        # Perform Fisher's exact test for downregulated genes
        logger.info("Performing Fisher's exact test for downregulated genes...")
        down_result = pathway_enrichment_test(downregulated_genes, pathway_genes, background_genes)
        
        # Save results to file
        results = save_results(up_result, down_result, output_file)
        
        # Print summary to console
        print("\n" + "=" * 80)
        print("GENE SET ENRICHMENT RESULTS")
        print("=" * 80)
        print("\nUpregulated genes enrichment:")
        print(f"P-value: {up_result['p_value']:.6f}")
        print(f"Odds ratio: {up_result['odds_ratio']:.3f}")
        print(f"Genes in pathway: {up_result['genes_in_pathway']}/{up_result['total_diff_genes']}")
        if up_result['overlap_genes']:
            print(f"Overlapping genes: {', '.join(up_result['overlap_genes'])}")
        
        print("\nDownregulated genes enrichment:")
        print(f"P-value: {down_result['p_value']:.6f}")
        print(f"Odds ratio: {down_result['odds_ratio']:.3f}")
        print(f"Genes in pathway: {down_result['genes_in_pathway']}/{down_result['total_diff_genes']}")
        if down_result['overlap_genes']:
            print(f"Overlapping genes: {', '.join(down_result['overlap_genes'])}")
        
        print(f"\nResults saved to: {output_file}")
        print("=" * 80)
        
        return 0
    
    except Exception as e:
        logger.error(f"Error during analysis: {e}")
        return 1

if __name__ == "__main__":
    try:
        exit_code = main()
        import sys
        sys.exit(exit_code)
    except Exception as e:
        logger.error(f"Unhandled exception: {e}")
        import sys
        sys.exit(1)
