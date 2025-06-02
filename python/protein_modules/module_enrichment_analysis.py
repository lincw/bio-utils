# -*- coding: utf-8 -*-
# Script: module_enrichment_analysis.py
# date_created: 2025-05-27T09:18:18+02:00
# date_modified: 2025-06-02T15:00:32+02:00

"""
Module Enrichment Analysis Script

This script calculates the statistical significance (using Fisher's exact test) 
of the overlap between a gene list and protein modules.

Usage:
    python module_enrichment_analysis.py [--gene_list <gene_list_file>] [--module_file <module_file>] [--output <output_file>]

Arguments:
    --gene_list: Path to a text file containing one gene per line (header will be skipped)
    --module_file: Path to a CSV file containing module data with columns: ID, process_name, Members
    --output: Path to save the output CSV file

Environment Variables (can be set in .env file):
    GENE_LIST_FILE: Default path to gene list file
    MODULE_FILE: Default path to module file
    OUTPUT_FILE: Default path to output file
    RESULTS_DIR: Base directory for results (used if specific paths not provided)
"""

import os
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
import logging
import argparse
from dotenv import load_dotenv

# Load environment variables from .env file if it exists
load_dotenv()

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def show_gene_list_template():
    """Show a template for the gene list file format."""
    template = """
# Gene list file template
# First line is treated as header and skipped
# Each subsequent line should contain one gene ID
Gene_ID
Gene1
Gene2
Gene3
...
"""
    logger.info("Expected gene list format:")
    for line in template.strip().split('\n'):
        logger.info(line)

def load_gene_list(file_path):
    """Load a gene list from a text file, one gene per line."""
    try:
        with open(file_path, 'r') as f:
            # Skip header
            next(f)
            genes = {line.strip() for line in f if line.strip()}
        logger.info(f"Loaded {len(genes)} genes from {file_path}")
        return genes
    except Exception as e:
        logger.error(f"Error loading gene list from {file_path}: {e}")
        show_gene_list_template()
        raise

def show_module_file_template():
    """Show a template for the module file format."""
    template = """
# Module file template (CSV format)
# Required columns: ID, process_name, Members
# Members column should contain comma-separated gene IDs

ID,process_name,Members
M1,Process 1,Gene1,Gene2,Gene3
M2,Process 2,Gene4,Gene5,Gene6
M3,Process 3,Gene1,Gene7,Gene8
...
"""
    logger.info("Expected module file format (CSV):")
    for line in template.strip().split('\n'):
        logger.info(line)

def load_modules(file_path):
    """Load module data from CSV file."""
    try:
        df = pd.read_csv(file_path)
        logger.info(f"Loaded module data from {file_path}: {df.shape[0]} rows, {df.shape[1]} columns")

        # Check if required columns exist
        required_cols = ['ID', 'process_name', 'Members']
        missing_cols = [col for col in required_cols if col not in df.columns]
        
        if missing_cols:
            logger.error(f"Required columns not found in module file: {', '.join(missing_cols)}")
            show_module_file_template()
            raise ValueError(f"Required columns not found in module file: {', '.join(missing_cols)}")
        
        # Use the Members column which contains the gene lists
        gene_col = 'Members'
        logger.info(f"Using '{gene_col}' as the gene column")
        
        return df, gene_col
    except Exception as e:
        logger.error(f"Error loading module data from {file_path}: {e}")
        if not str(e).startswith("Required columns not found"):
            show_module_file_template()
        raise

def get_module_genes(df, module_id, gene_col):
    """Get the set of genes for a specific module."""
    module_rows = df[df['ID'] == module_id]
    if module_rows.empty:
        return set()
    
    # Extract genes from the Members column (comma-separated list)
    genes_str = module_rows.iloc[0][gene_col]
    if pd.isna(genes_str) or not genes_str:
        return set()
    
    # Split by comma and strip whitespace
    genes = {gene.strip() for gene in genes_str.split(',')}
    return genes

def perform_fisher_test(module_genes, test_genes, background_genes):
    """
    Perform Fisher's exact test for enrichment.
    
    Args:
        module_genes: Set of genes in the module
        test_genes: Set of test genes (up/down/merged)
        background_genes: Set of all genes (universe)
    
    Returns:
        tuple: (overlap_count, p_value, odds_ratio)
    """
    # Calculate overlap
    overlap = module_genes.intersection(test_genes)
    overlap_count = len(overlap)
    
    # Construct contingency table
    a = overlap_count  # In module and in test set
    b = len(test_genes) - a  # In test set but not in module
    c = len(module_genes) - a  # In module but not in test set
    d = len(background_genes) - a - b - c  # Neither in module nor in test set
    
    # Perform Fisher's exact test
    contingency_table = [[a, b], [c, d]]
    odds_ratio, p_value = fisher_exact(contingency_table, alternative='greater')
    
    return overlap_count, p_value, odds_ratio, list(overlap)

def print_help():
    """Print a more detailed help message with examples."""
    help_text = """
==============================================================================
MODULE ENRICHMENT ANALYSIS
==============================================================================

This script performs enrichment analysis of gene lists against modules using
Fisher's exact test.

USAGE:
  python module_enrichment_analysis.py [options]

OPTIONS:
  --gene_list    Path to gene list file (one gene per line, header skipped)
  --module_file  Path to module file (CSV with ID, process_name, Members columns)
  --output       Path to save output CSV file
  --help         Show this help message

ENVIRONMENT VARIABLES (can be set in .env file):
  GENE_LIST_FILE  Default path to gene list file
  MODULE_FILE     Default path to module file
  OUTPUT_FILE     Default path to output file
  RESULTS_DIR     Base directory for results (used if specific paths not provided)

EXAMPLES:
  # Using command-line arguments
  python module_enrichment_analysis.py \
    --gene_list data/my_genes.txt \
    --module_file data/modules.csv \
    --output results/enrichment_results.csv

  # Using .env file
  # Create a .env file with the following content:
  #   GENE_LIST_FILE=data/my_genes.txt
  #   MODULE_FILE=data/modules.csv
  #   OUTPUT_FILE=results/enrichment_results.csv
  python module_enrichment_analysis.py

  # Mix of .env and command-line (command-line takes precedence)
  python module_enrichment_analysis.py --gene_list data/different_genes.txt

FILE FORMATS:
  Gene List: Text file with one gene per line (first line skipped as header)
  Module File: CSV with columns ID, process_name, and Members
              (Members should contain comma-separated gene IDs)
==============================================================================
"""
    print(help_text)

def main():
    """Main function to run the analysis."""
    
    # Get default paths from environment variables
    results_dir = os.getenv('RESULTS_DIR')
    default_gene_file = os.getenv('GENE_LIST_FILE')
    default_module_file = os.getenv('MODULE_FILE')
    default_output_file = os.getenv('OUTPUT_FILE')
    
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Module Enrichment Analysis', add_help=False)
    parser.add_argument('--gene_list', help='Path to gene list file (one gene per line)')
    parser.add_argument('--module_file', help='Path to module file (CSV with ID, process_name, Members columns)')
    parser.add_argument('--output', help='Path to save output CSV file')
    parser.add_argument('--help', '-h', action='store_true', help='Show detailed help message')
    args = parser.parse_args()
    
    # Show help and exit if requested
    if args.help:
        print_help()
        return
    
    # Command-line arguments take precedence over environment variables
    gene_file = args.gene_list or default_gene_file
    module_file = args.module_file or default_module_file
    output_file = args.output or default_output_file
    
    # If results_dir is set and relative paths are provided, make them absolute
    if results_dir:
        if gene_file and not os.path.isabs(gene_file):
            gene_file = os.path.join(results_dir, gene_file)
        if module_file and not os.path.isabs(module_file):
            module_file = os.path.join(results_dir, module_file)
        if output_file and not os.path.isabs(output_file):
            output_file = os.path.join(results_dir, output_file)
    
    # Check if required files are specified
    missing_paths = []
    if not gene_file:
        missing_paths.append("gene_list")
    if not module_file:
        missing_paths.append("module_file")
    if not output_file:
        missing_paths.append("output_file")
        
    if missing_paths:
        print("\n" + "=" * 80)
        print("ERROR: MISSING REQUIRED FILE PATHS")
        print("=" * 80)
        print("\nThe following required paths are missing:")
        for path in missing_paths:
            print(f"  - {path}")
            
        print("\nYou can provide these paths in two ways:")
        print("  1. Command-line arguments:")
        print("     python module_enrichment_analysis.py --gene_list <file> --module_file <file> --output <file>")
        print("\n  2. Environment variables in a .env file:")
        print("     GENE_LIST_FILE=<file>")
        print("     MODULE_FILE=<file>")
        print("     OUTPUT_FILE=<file>")
        print("     RESULTS_DIR=<base_directory>  # Optional, for relative paths")
        
        print("\nFor more information on file formats, run:")
        print("  python module_enrichment_analysis.py --help")
        print("=" * 80)
        
        # Still log the error for debugging purposes
        logger.error(f"Missing required paths: {', '.join(missing_paths)}")
        raise ValueError(f"Missing required paths: {', '.join(missing_paths)}")
    
    # Validate file paths
    missing_files = []
    if not os.path.exists(gene_file):
        missing_files.append(("gene_list", gene_file))
    if not os.path.exists(module_file):
        missing_files.append(("module_file", module_file))
        
    if missing_files:
        print("\n" + "=" * 80)
        print("ERROR: INPUT FILES NOT FOUND")
        print("=" * 80)
        print("\nThe following input files could not be found:")
        for file_type, path in missing_files:
            print(f"  - {file_type}: {path}")
        
        print("\nPlease check that the file paths are correct and the files exist.")
        
        # Show templates for missing files
        if not os.path.exists(gene_file):
            print("\n" + "-" * 40)
            print("GENE LIST TEMPLATE:")
            print("-" * 40)
            template = """
# Gene list file template
# First line is treated as header and skipped
# Each subsequent line should contain one gene ID
Gene_ID
Gene1
Gene2
Gene3
..."""
            print(template)
        
        if not os.path.exists(module_file):
            print("\n" + "-" * 40)
            print("MODULE FILE TEMPLATE (CSV):")
            print("-" * 40)
            template = """
# Module file template (CSV format)
# Required columns: ID, process_name, Members
# Members column should contain comma-separated gene IDs

ID,process_name,Members
M1,Process 1,Gene1,Gene2,Gene3
M2,Process 2,Gene4,Gene5,Gene6
M3,Process 3,Gene1,Gene7,Gene8
..."""
            print(template)
            
        print("=" * 80)
        
        # Still log the error for debugging purposes
        logger.error(f"Input files not found: {[path for _, path in missing_files]}")
        raise ValueError("Input files not found")
    
    # Create output directory if it doesn't exist
    output_dir = os.path.dirname(output_file)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
        logger.info(f"Created output directory: {output_dir}")
    
    # Load data
    logger.info("Loading data...")
    query_genes = load_gene_list(gene_file)
    logger.info(f"Loaded {len(query_genes)} genes from {gene_file}")
    
    module_df, gene_col = load_modules(module_file)
    
    # Get all unique modules
    unique_modules = module_df['ID'].unique()
    logger.info(f"Found {len(unique_modules)} unique modules")
    
    # Extract all genes from the Members column to create the background set
    background_genes = set()
    for _, row in module_df.iterrows():
        members_str = row[gene_col]
        if pd.notna(members_str) and members_str:
            genes = {gene.strip() for gene in members_str.split(',')}
            background_genes.update(genes)
    logger.info(f"Background set contains {len(background_genes)} genes")
    
    # Check if query genes are in background
    query_in_bg = query_genes.intersection(background_genes)
    logger.info(f"{len(query_in_bg)}/{len(query_genes)} query genes found in background")
    
    # Only use genes that are in the background
    query_genes = query_in_bg
    
    # Initialize results list
    results = []
    
    # Process each module
    logger.info("Processing modules...")
    for module_id in unique_modules:
        # Get module name
        module_name = module_df[module_df['ID'] == module_id]['process_name'].iloc[0]
        
        # Get genes in this module
        module_genes = get_module_genes(module_df, module_id, gene_col)
        
        # Skip modules with no genes
        if not module_genes:
            logger.warning(f"Module {module_id} ({module_name}) has no genes, skipping")
            continue
        
        # Perform Fisher's exact test for query genes
        overlap_count, p_value, odds_ratio, overlap_genes = perform_fisher_test(
            module_genes, query_genes, background_genes
        )
        
        # Add to results
        results.append({
            'module_id': module_id,
            'module_name': module_name,
            'module_size': len(module_genes),
            'overlap_count': overlap_count,
            'p_value': p_value,
            'odds_ratio': odds_ratio,
            'overlap_genes': ';'.join(overlap_genes) if overlap_genes else ''
        })
    
    # Convert to DataFrame
    results_df = pd.DataFrame(results)
    
    # Apply FDR correction
    if not results_df.empty:
        # FDR correction
        _, fdr_values, _, _ = multipletests(results_df['p_value'].values, method='fdr_bh')
        results_df['fdr'] = fdr_values
    else:
        logger.warning("No results to apply FDR correction")
    
    # Sort by p-value
    results_df = results_df.sort_values('p_value')
    
    # Save results
    results_df.to_csv(output_file, index=False)
    logger.info(f"Results saved to {output_file}")
    
    # Save significant modules (p-value < 0.05)
    sig_file = os.path.splitext(output_file)[0] + '_significant_pvalue0.05.csv'
    sig_df = results_df[results_df['p_value'] < 0.05].copy()
    if not sig_df.empty:
        sig_df.to_csv(sig_file, index=False)
        logger.info(f"Saved {len(sig_df)} significant modules to {sig_file}")
    else:
        logger.warning("No modules with p-value < 0.05")
    
    # Print summary
    logger.info("\nSummary of significant results:")
    logger.info(f"Modules with p-value < 0.05: {(results_df['p_value'] < 0.05).sum()}")
    logger.info(f"Modules with FDR < 0.05: {(results_df['fdr'] < 0.05).sum()}")

if __name__ == "__main__":
    try:
        main()
    except ValueError as e:
        # The detailed error message has already been printed
        # Just exit with error code
        import sys
        sys.exit(1)

