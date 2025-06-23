# -*- coding: utf-8 -*-
# script: convert_uniprot_to_gene.py
# date_created: 2025-06-23T11:02:41+02:00
# date_modified: 2025-06-23T12:39:11+02:00

"""
UniProt ID to Gene Symbol Converter

This script converts UniProt IDs to gene symbols in OmniPath interaction files
using the UniProt API. It includes a progress bar to track conversion progress.

Usage:
    python convert_uniprot_to_gene.py input.csv [output.csv] --organism ORGANISM [options]

Arguments:
    input.csv       Path to the input CSV file with UniProt IDs
    output.csv      Optional path to the output CSV file (default: input_file_genesymbols.csv)
    --organism      Organism name (default: mouse)
    --keep-both     Keep both UniProt IDs and gene symbols in the output
    --verbose       Enable verbose logging of conversion results
    --list-organisms List all available organisms in OmniPath and exit
"""

import os
import sys
import time
import argparse
import logging
import re
import pandas as pd
import requests
from datetime import datetime
from tqdm import tqdm
from io import StringIO


# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)


class UniProtConverter:
    """Class for converting UniProt IDs to gene symbols using UniProt API."""
    
    # Base URLs for UniProt API
    UNIPROT_API_URL = "https://rest.uniprot.org/uniprotkb/search"
    
    # Regular expression to match UniProt IDs
    UNIPROT_ID_PATTERN = re.compile(r'^[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$')
    
    def __init__(self, organism="mouse", taxonomy_id=None, verbose=False):
        """
        Initialize the converter.
        
        Args:
            organism: Organism name (e.g., human, mouse, rat)
            verbose: Whether to enable verbose logging
        """
        self.organism = organism.lower()
        # Use provided taxonomy ID if given, else fetch dynamically
        self.taxonomy_id = taxonomy_id or self._get_taxonomy_id(self.organism)
        self.verbose = verbose
        
        if not self.taxonomy_id:
            raise ValueError(f"Could not find taxonomy ID for organism: {organism}")
        
        # Cache for UniProt ID to gene symbol mappings
        self.id_to_gene_cache = {}
        
        # Statistics for conversion results
        self.stats = {
            "total_ids": 0,
            "converted_ids": 0,
            "already_gene_symbols": 0,
            "failed_ids": 0,
            "not_found_ids": []
        }
    
    def _get_taxonomy_id(self, organism):
        """
        Get the taxonomy ID for the given organism using UniProt's taxonomy API.
        
        Args:
            organism: Organism name (e.g., human, mouse, rat)
            
        Returns:
            Taxonomy ID as a string, or None if not found
            
        Raises:
            ValueError: If the organism cannot be found or if there's an API error
        """
        try:
            # Query the UniProt taxonomy API with more specific parameters
            url = "https://rest.uniprot.org/taxonomy/search"
            
            # Use a more specific query format
            query = f"name:{organism}"
            
            params = {
                "query": query,
                "fields": "id,scientific_name,common_name",
                "format": "tsv"
            }
            
            logger.debug(f"Querying UniProt taxonomy API for organism: {organism}")
            response = requests.get(url, params=params, timeout=10)
            response.raise_for_status()
            
            # Parse the response
            lines = response.text.strip().split('\n')
            if len(lines) <= 1:  # Only header or empty
                error_msg = f"No taxonomy found for organism: {organism}. Please check the organism name and try again."
                logger.error(error_msg)
                raise ValueError(error_msg)
                
            # Skip the header line and get the first result
            # Format: Taxon ID\tScientific name\tCommon name
            parts = lines[1].split('\t')
            if len(parts) >= 1:
                taxonomy_id = parts[0]
                scientific_name = parts[1] if len(parts) > 1 else organism
                logger.info(f"Found taxonomy ID for {organism} ({scientific_name}): {taxonomy_id}")
                return taxonomy_id
            
            error_msg = f"Invalid response format from UniProt taxonomy API for organism: {organism}"
            logger.error(error_msg)
            raise ValueError(error_msg)
            
        except requests.exceptions.RequestException as e:
            error_msg = f"Error connecting to UniProt taxonomy API: {e}. Please check your internet connection."
            logger.error(error_msg)
            raise ValueError(error_msg)
    
    # This method is now implemented above
    
    def convert_ids_batch(self, uniprot_ids):
        """
        Convert a batch of UniProt IDs to gene symbols using UniProt's search API.
        This is a more reliable method than the ID mapping service for our use case.
        
        Args:
            uniprot_ids: List of UniProt IDs to convert
            
        Returns:
            Dictionary mapping UniProt IDs to gene symbols
        """
        # Filter out IDs that are already in the cache
        ids_to_query = [id for id in uniprot_ids if id not in self.id_to_gene_cache]
        
        # Validate and separate UniProt IDs from potential gene symbols
        valid_uniprot_ids = []
        potential_gene_symbols = []
        
        for id in ids_to_query:
            if id and isinstance(id, str) and self.UNIPROT_ID_PATTERN.match(id):
                valid_uniprot_ids.append(id)
            else:
                # If it doesn't match UniProt ID pattern, assume it might be a gene symbol already
                if id and isinstance(id, str):
                    potential_gene_symbols.append(id)
                    self.id_to_gene_cache[id] = id  # Keep as is
                    self.stats["already_gene_symbols"] += 1
                    if self.verbose:
                        logger.debug(f"ID '{id}' does not match UniProt ID pattern, treating as gene symbol")
        
        # Process valid UniProt IDs
        if valid_uniprot_ids:
            # Update statistics for valid UniProt IDs
            self.stats["total_ids"] += len(valid_uniprot_ids)
            
            # Query the API for valid UniProt IDs
            result = self._query_uniprot_search_api(valid_uniprot_ids)
            
            # Update the cache with the results
            self.id_to_gene_cache.update(result)
            
            # Update statistics
            self.stats["converted_ids"] += len(result)
            self.stats["failed_ids"] += len(valid_uniprot_ids) - len(result)
            
            # Track IDs that couldn't be converted
            not_found = [id for id in valid_uniprot_ids if id not in result]
            if not_found and self.verbose:
                self.stats["not_found_ids"].extend(not_found)
                logger.debug(f"Could not convert {len(not_found)} UniProt IDs: {', '.join(not_found[:5])}" + 
                           (f" and {len(not_found)-5} more" if len(not_found) > 5 else ""))
        
        # Combine cache results with new results
        combined_results = {}
        for id in uniprot_ids:
            if id in self.id_to_gene_cache:
                combined_results[id] = self.id_to_gene_cache[id]
            elif id and isinstance(id, str):
                combined_results[id] = id  # Use ID as fallback
        
        return combined_results
    
    def _query_uniprot_search_api(self, uniprot_ids):
        """
        Query the UniProt search API for gene symbols.
        
        Args:
            uniprot_ids: List of UniProt IDs to query
            
        Returns:
            Dictionary mapping UniProt IDs to gene symbols
        """
        result = {}
        
        # Process in smaller batches to avoid API limitations
        max_ids_per_request = 10  # Smaller batch size to avoid errors
        
        for i in range(0, len(uniprot_ids), max_ids_per_request):
            batch = uniprot_ids[i:i+max_ids_per_request]
            
            # Construct the query - use OR for multiple IDs
            if len(batch) == 1:
                query = f"accession:{batch[0]}"
            else:
                query = "(" + " OR ".join([f"accession:{id}" for id in batch]) + ")"
            
            # Add taxonomy filter if available
            if self.taxonomy_id:
                query += f" AND organism_id:{self.taxonomy_id}"
            
            # Set up the parameters
            params = {
                "query": query,
                "format": "tsv",
                "fields": "accession,gene_names",
                "size": max_ids_per_request
            }
            
            retry_count = 0
            max_retries = 3
            
            while retry_count < max_retries:
                try:
                    # Make the request
                    response = requests.get(self.UNIPROT_API_URL, params=params)
                    response.raise_for_status()
                    
                    # Parse the TSV response
                    lines = response.text.strip().split('\n')
                    if len(lines) <= 1:  # Only header or empty
                        logger.debug(f"No results found for query: {query}")
                        break
                    
                    # Skip the header line
                    for line in lines[1:]:
                        parts = line.split('\t')
                        if len(parts) >= 2:
                            uniprot_id = parts[0]
                            # Take the first gene name if multiple are returned
                            gene_names = parts[1].strip()
                            if gene_names:
                                # Use the first gene name (primary gene symbol)
                                gene_symbol = gene_names.split()[0]
                                result[uniprot_id] = gene_symbol
                            else:
                                # If no gene name, use the UniProt ID as fallback
                                result[uniprot_id] = uniprot_id
                    
                    # Success, break the retry loop
                    break
                    
                except requests.exceptions.RequestException as e:
                    retry_count += 1
                    if retry_count >= max_retries:
                        logger.error(f"Error querying UniProt API after {max_retries} attempts: {e}")
                    else:
                        logger.warning(f"Retry {retry_count}/{max_retries} after error: {e}")
                        time.sleep(2 * retry_count)  # Exponential backoff
            
            # Add a small delay between batches to avoid rate limiting
            time.sleep(1)
        
        return result
    
    # Removed the ID mapping job methods as we're now using the search API directly
    
    def convert_dataframe(self, df, source_col="source", target_col="target", batch_size=50):
        """
        Convert UniProt IDs in a DataFrame to gene symbols.
        
        Args:
            df: DataFrame containing UniProt IDs
            source_col: Column name for source UniProt IDs
            target_col: Column name for target UniProt IDs
            batch_size: Number of IDs to convert in each batch
            
        Returns:
            Tuple of (gene_df, both_df) where:
              - gene_df: DataFrame with UniProt IDs replaced by gene symbols
              - both_df: DataFrame with both UniProt IDs and gene symbols
        """
        # Get unique UniProt IDs from both source and target columns, handling non-string values
        all_ids = set()
        for col in [source_col, target_col]:
            # Convert column to string type to handle mixed types
            unique_values = df[col].astype(str).unique()
            # Filter out NaN, None, etc.
            valid_ids = [id for id in unique_values if id and id.lower() not in ['nan', 'none', 'null', '']]  
            all_ids.update(valid_ids)
        
        unique_ids = list(all_ids)
        logger.info(f"Found {len(unique_ids)} unique IDs to convert")
        
        # Convert all unique IDs in batches with progress bar
        id_to_gene = {}
        with tqdm(total=len(unique_ids), desc="Converting UniProt IDs", unit="ids") as pbar:
            for i in range(0, len(unique_ids), batch_size):
                batch = unique_ids[i:i+batch_size]
                batch_result = self.convert_ids_batch(batch)
                id_to_gene.update(batch_result)
                pbar.update(len(batch))
        
        # Create a copy of the DataFrame for gene symbols
        gene_df = df.copy()
        
        # Replace UniProt IDs with gene symbols, handling non-string values
        def safe_map(x):
            if isinstance(x, str):
                return id_to_gene.get(x, x)
            return x
        
        gene_df[source_col] = gene_df[source_col].map(safe_map)
        gene_df[target_col] = gene_df[target_col].map(safe_map)
        
        # Create a copy with both UniProt IDs and gene symbols
        both_df = df.copy()
        both_df[f"{source_col}_gene"] = both_df[source_col].map(safe_map)
        both_df[f"{target_col}_gene"] = both_df[target_col].map(safe_map)
        
        # Log conversion statistics
        logger.info(f"Conversion statistics:")
        logger.info(f"  Total unique IDs: {self.stats['total_ids']}")
        logger.info(f"  Already gene symbols: {self.stats['already_gene_symbols']}")
        
        if self.stats['total_ids'] > 0:
            conversion_rate = (self.stats['converted_ids'] / self.stats['total_ids']) * 100
            failure_rate = (self.stats['failed_ids'] / self.stats['total_ids']) * 100
            logger.info(f"  Successfully converted: {self.stats['converted_ids']} ({conversion_rate:.1f}% success rate)")
            logger.info(f"  Failed to convert: {self.stats['failed_ids']} ({failure_rate:.1f}% failure rate)")
        
        if self.verbose and self.stats["not_found_ids"]:
            sample_ids = self.stats["not_found_ids"][:5]
            remaining = len(self.stats["not_found_ids"]) - len(sample_ids)
            logger.info(f"  Sample of IDs that couldn't be converted: {', '.join(sample_ids)}" + 
                      (f" and {remaining} more" if remaining > 0 else ""))
            
        return gene_df, both_df

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Convert UniProt IDs to gene symbols in interaction files")
    
    # Create a parent parser for common arguments
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument("-g", "--organism", default="mouse", help="Organism name (default: mouse)")
    parent_parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose output")
    
    
    # Add subparsers for different modes
    subparsers = parser.add_subparsers(dest="command", help="Command to run")
    
    # Convert command (default)
    convert_parser = subparsers.add_parser("convert", parents=[parent_parser], help="Convert UniProt IDs to gene symbols")
    convert_parser.add_argument("input_file", help="Path to the input CSV file with UniProt IDs")
    convert_parser.add_argument("-o", "--output-file", help="Path to the output CSV file with gene symbols")
    convert_parser.add_argument("-k", "--keep-both", action="store_true", help="Keep both UniProt IDs and gene symbols in output")
    convert_parser.add_argument("-b", "--batch-size", type=int, default=50, help="Batch size for API requests (default: 50)")
    convert_parser.add_argument("-t", "--taxonomy-id", help="NCBI taxonomy ID to use directly (skips lookup)")
    
    # Parse arguments
    args = parser.parse_args()
    
    # If no command is specified, default to convert
    if not args.command:
        parser.error("the following arguments are required: input_file")
    
    return args


def main():
    """Main function."""
    args = parse_arguments()
    
    # Determine the output file path if not provided
    if not args.output_file:
        # Simply add _genesymbols before the extension
        base_name, ext = os.path.splitext(args.input_file)
        args.output_file = f"{base_name}_genesymbols{ext}"
    
    logger.info(f"Converting UniProt IDs to gene symbols for organism: {args.organism}")
    logger.info(f"Input file: {args.input_file}")
    logger.info(f"Output file: {args.output_file}")
    logger.info(f"Verbose mode: {'enabled' if args.verbose else 'disabled'}")
    logger.info(f"Batch size: {args.batch_size}")
    logger.info(f"Keep both IDs and symbols: {'yes' if args.keep_both else 'no'}")
    
    # Read the input CSV file
    try:
        logger.info(f"Reading input file: {args.input_file}")
        df = pd.read_csv(args.input_file, low_memory=False)
        logger.info(f"Read {len(df)} rows from input file")
    except Exception as e:
        logger.error(f"Error reading input file: {e}")
        sys.exit(1)
        
    # Check if the required columns exist
    required_columns = ["source", "target"]
    missing_columns = [col for col in required_columns if col not in df.columns]
    if missing_columns:
        logger.error(f"Missing required columns in input file: {missing_columns}")
        logger.error(f"The input file must contain 'source' and 'target' columns with UniProt IDs")
        sys.exit(1)
    
    # Check for empty dataframe
    if df.empty:
        logger.error("Input file contains no data")
        sys.exit(1)
        
    # Display sample of input data
    if args.verbose:
        logger.info("Sample of input data (first 5 rows):")
        logger.info(df.head(5))
        
    # Initialize the converter
    try:
        logger.info(f"Initializing UniProt converter for organism: {args.organism}")
        converter = UniProtConverter(organism=args.organism, taxonomy_id=args.taxonomy_id, verbose=args.verbose)
        logger.info(f"Successfully initialized converter with taxonomy ID: {converter.taxonomy_id}")
    except ValueError as e:
        logger.error(f"Error initializing converter: {e}")
        logger.error("Please check your internet connection and the organism name")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Unexpected error initializing converter: {e}")
        sys.exit(1)
    
    # Convert UniProt IDs to gene symbols
    try:
        logger.info("Starting conversion of UniProt IDs to gene symbols...")
        gene_df, both_df = converter.convert_dataframe(df, batch_size=args.batch_size)
        logger.info("Conversion completed successfully")
    except requests.exceptions.RequestException as e:
        logger.error(f"Network error during conversion: {e}")
        logger.error("Please check your internet connection and try again")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Error during conversion: {e}")
        sys.exit(1)
        
    # Save the converted DataFrame to a CSV file
    try:
        logger.info(f"Saving converted data to: {args.output_file}")
        if args.keep_both:
            both_df.to_csv(args.output_file, index=False)
            logger.info(f"Saved interactions with both UniProt IDs and gene symbols to {args.output_file}")
        else:
            gene_df.to_csv(args.output_file, index=False)
            logger.info(f"Saved interactions with gene symbols to {args.output_file}")
        
        # Display sample of output data
        if args.verbose:
            logger.info("Sample of output data (first 5 rows):")
            if args.keep_both:
                logger.info(both_df.head(5))
            else:
                logger.info(gene_df.head(5))
            
    except Exception as e:
        logger.error(f"Error saving output file: {e}")
        sys.exit(1)


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        logger.info("\nConversion interrupted by user")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Unexpected error: {e}")
        import traceback
        logger.error(traceback.format_exc())
        sys.exit(1)
