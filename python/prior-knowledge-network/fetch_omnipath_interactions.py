# -*- coding: utf-8 -*-
# script: fetch_omnipath_interactions.py
# date_created: 2025-06-23T09:33:22+02:00
# date_modified: 2025-06-23T09:35:01+02:00

"""
OmniPath Interaction Downloader

This script provides a modular interface for:
1. Exploring available interaction sources in OmniPath
2. Exploring available organisms in OmniPath
3. Downloading interaction data with custom filters

Usage:
    python fetch_omnipath_interactions.py list-sources
    python fetch_omnipath_interactions.py list-organisms
    python fetch_omnipath_interactions.py download --organism human --datasets omnipath,dorothea --output my_interactions.csv
    python fetch_omnipath_interactions.py download-all --organism human --output all_interactions.csv
    python fetch_omnipath_interactions.py download-pt --organism human --output pt_interactions.csv
"""

import sys
import inspect
import argparse
import importlib
import pandas as pd
from typing import List, Optional, Union, Dict, Type, Any
from omnipath.constants import InteractionDataset, Organism
import omnipath.interactions as interactions_module
from omnipath.interactions import AllInteractions, PostTranslational

# Dynamically discover all interaction dataset classes
def discover_dataset_classes() -> Dict[str, Type]:
    """Dynamically discover all interaction dataset classes from the omnipath.interactions module."""
    # Get all classes from the interactions module
    classes = inspect.getmembers(interactions_module, inspect.isclass)
    
    # Filter classes that are defined in the interactions module (not imported)
    dataset_classes = {}
    for name, cls in classes:
        # Skip abstract classes and base classes
        if inspect.isabstract(cls) or name in ('ABC', 'InteractionRequest', 'CommonParamFilter'):
            continue
        
        # Check if the class has a get method
        if hasattr(cls, 'get') and callable(getattr(cls, 'get')):
            # Try to get the dataset name from the class
            try:
                # Create an instance to check if it works
                instance = cls()
                # Use the class name in lowercase as the key
                key = name.lower()
                dataset_classes[key] = cls
            except Exception:
                # Skip classes that can't be instantiated
                pass
    
    return dataset_classes

# Dynamically discover dataset classes
DATASET_CLASSES = discover_dataset_classes()

class OmnipathExplorer:
    """Class for exploring OmniPath resources."""
    
    @staticmethod
    def list_interaction_sources() -> None:
        """List all available interaction sources in OmniPath."""
        print("Available interaction sources in OmniPath:")
        sources = []
        for source in InteractionDataset:
            sources.append((source.name, source.value))
            
        # Sort alphabetically by name for better readability
        sources.sort(key=lambda x: x[0])
        
        for name, value in sources:
            print(f"- {name}: {value}")
        print(f"\nTotal number of interaction sources: {len(InteractionDataset.__members__)}")
        
        print("\nAvailable dataset-specific downloaders:")
        for name in sorted(DATASET_CLASSES.keys()):
            print(f"- {name}")
    
    @staticmethod
    def list_available_organisms() -> None:
        """List all available organisms in OmniPath."""
        print("\nAvailable organisms in OmniPath:")
        for organism in Organism:
            print(f"- {organism.name}: {organism.value} (NCBI Taxonomy ID: {organism.code})")
        print(f"\nTotal number of organisms: {len(Organism.__members__)}")
    
    @staticmethod
    def get_organism_by_name(name: str) -> Organism:
        """
        Get Organism enum by name (case-insensitive).
        
        Parameters:
        -----------
        name : str
            Name of the organism (e.g., 'human', 'mouse', 'rat')
            
        Returns:
        --------
        Organism
            The corresponding Organism enum
            
        Raises:
        -------
        ValueError
            If the organism name is not valid
        """
        name = name.upper()
        try:
            return Organism[name]
        except KeyError:
            valid_organisms = [o.name.lower() for o in Organism]
            raise ValueError(f"Invalid organism: {name.lower()}. Valid options are: {', '.join(valid_organisms)}")
    
    @staticmethod
    def get_datasets_by_names(names: List[str]) -> List[InteractionDataset]:
        """
        Get InteractionDataset enums by names.
        
        Parameters:
        -----------
        names : List[str]
            List of dataset names
            
        Returns:
        --------
        List[InteractionDataset]
            List of corresponding InteractionDataset enums
            
        Raises:
        -------
        ValueError
            If any dataset name is not valid
        """
        datasets = []
        valid_datasets = {ds.value: ds for ds in InteractionDataset}
        
        for name in names:
            if name in valid_datasets:
                datasets.append(valid_datasets[name])
            else:
                raise ValueError(f"Invalid dataset: {name}. Valid options are: {', '.join(valid_datasets.keys())}")
        
        return datasets


class OmnipathDownloader:
    """Class for downloading interaction data from OmniPath."""
    
    @staticmethod
    def download_interactions(organism: str, 
                           dataset_type: str = "all",
                           specific_dataset: Optional[str] = None,
                           include_datasets: Optional[List[str]] = None, 
                           exclude_datasets: Optional[List[str]] = None, 
                           output_file: Optional[str] = None) -> pd.DataFrame:
        """Unified method to download interactions with flexible options.
        
        Args:
            organism: Organism name (human, mouse, rat)
            dataset_type: Type of dataset to download ('all', 'specific', 'post_translational')
            specific_dataset: Name of a specific dataset to download (used when dataset_type='specific')
            include_datasets: List of dataset names to include (used when dataset_type='all')
            exclude_datasets: List of dataset names to exclude (used when dataset_type='all')
            output_file: Directory path to save the downloaded data (optional)
            
        Returns:
            DataFrame containing the downloaded interactions
        """
        # Validate organism name
        organism_enum = OmnipathExplorer.get_organism_by_name(organism)
        if not organism_enum:
            print(f"Error: Invalid organism name '{organism}'. Valid options are: {', '.join([o.value for o in Organism])}")
            return pd.DataFrame()
        
        interactions = pd.DataFrame()
        
        # Handle different dataset types
        if dataset_type == "specific" and specific_dataset:
            # Validate dataset name
            if specific_dataset.lower() not in DATASET_CLASSES:
                print(f"Error: Invalid dataset name '{specific_dataset}'. Valid options are: {', '.join(DATASET_CLASSES.keys())}")
                return pd.DataFrame()
            
            # Get the dataset class
            dataset_class = DATASET_CLASSES[specific_dataset.lower()]
            
            print(f"Downloading {specific_dataset} interactions for {organism_enum.name}...")
            
            # Download interactions - first instantiate the class, then call get() with organism parameter
            interactions = dataset_class().get(organism=organism_enum.value)
            
            print(f"Downloaded {len(interactions)} interactions from {specific_dataset}.")
            
        elif dataset_type == "post_translational":
            print(f"Downloading post-translational interactions for {organism_enum.name}...")
            
            # Download interactions
            interactions = PostTranslational().get(organism=organism_enum.value)
            
            print(f"Downloaded {len(interactions)} post-translational interactions.")
            
        else:  # dataset_type == "all" or any other value
            # Convert include dataset names to enum values if provided
            include_datasets_enum = None
            if include_datasets:
                include_datasets_enum = OmnipathExplorer.get_datasets_by_names(include_datasets)
                if not include_datasets_enum:
                    return pd.DataFrame()
            
            # Convert exclude dataset names to enum values if provided
            exclude_datasets_enum = None
            if exclude_datasets:
                exclude_datasets_enum = OmnipathExplorer.get_datasets_by_names(exclude_datasets)
                if not exclude_datasets_enum and exclude_datasets:
                    return pd.DataFrame()
            
            print(f"Downloading interactions for {organism_enum.name}...")
            if include_datasets_enum:
                print(f"Including datasets: {', '.join([d.name for d in include_datasets_enum])}")
            if exclude_datasets_enum:
                print(f"Excluding datasets: {', '.join([d.name for d in exclude_datasets_enum])}")
            
            # Download interactions
            interactions = AllInteractions().get(
                organism=organism_enum.value,
                include=include_datasets_enum,
                exclude=exclude_datasets_enum
            )
            
            print(f"Downloaded {len(interactions)} interactions.")
        
        if not interactions.empty:
            print(f"Columns in the interactions dataframe: {', '.join(interactions.columns)}")
            
            # Generate filename with timestamp
            import os
            from datetime import datetime
            
            # Generate timestamp in format YYYYMMDDTHHMMSS
            timestamp = datetime.now().strftime("%Y%m%dT%H%M%S")
            # Create filename based on organism and dataset type
            if dataset_type == "specific" and specific_dataset:
                file_prefix = f"{organism.lower()}_{specific_dataset.lower()}"
            elif dataset_type == "post_translational":
                file_prefix = f"{organism.lower()}_post_translational"
            else:
                file_prefix = f"{organism.lower()}_all"
            
            filename = f"{file_prefix}_interactions_{timestamp}.csv"
            
            # Determine the output directory
            if output_file:
                # Check if output_file is a directory or a file path
                if os.path.isdir(output_file) or output_file.endswith('/'):
                    # It's a directory, use it as the output directory
                    output_dir = output_file
                    file_path = os.path.join(output_dir, filename)
            else:
                # No output specified, use current directory
                file_path = os.path.join(os.getcwd(), filename)
            
            # Save the interactions to the file
            interactions.to_csv(file_path, index=False)
            print(f"Data saved to: {file_path}")
        
        return interactions


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="OmniPath Interaction Downloader")
    subparsers = parser.add_subparsers(dest="command", help="Command to execute")
    
    # List sources command
    subparsers.add_parser("list-sources", help="List all available interaction sources")
    
    # List organisms command
    subparsers.add_parser("list-organisms", help="List all available organisms")
    
    # Download command with custom filters
    download_parser = subparsers.add_parser("download", help="Download interactions with custom filters")
    download_parser.add_argument("--organism", type=str, default="human", help="Organism name (human, mouse, rat)")
    download_parser.add_argument("--include", type=str, help="Comma-separated list of datasets to include")
    download_parser.add_argument("--exclude", type=str, help="Comma-separated list of datasets to exclude")
    download_parser.add_argument("--output", type=str, help="Directory path to save the downloaded data")
    
    # Download specific dataset command
    download_dataset_parser = subparsers.add_parser("download-dataset", help="Download interactions from a specific dataset")
    download_dataset_parser.add_argument("--dataset", type=str, required=True, help="Dataset name")
    download_dataset_parser.add_argument("--organism", type=str, default="human", help="Organism name (human, mouse, rat)")
    download_dataset_parser.add_argument("--output", type=str, help="Directory path to save the downloaded data")
    
    # Download all interactions command
    download_all_parser = subparsers.add_parser("download-all", help="Download all interactions")
    download_all_parser.add_argument("--organism", type=str, default="human", help="Organism name (human, mouse, rat)")
    download_all_parser.add_argument("--output", type=str, help="Directory path to save the downloaded data")
    
    # Download post-translational interactions command
    download_pt_parser = subparsers.add_parser("download-pt", help="Download post-translational interactions")
    download_pt_parser.add_argument("--organism", type=str, default="human", help="Organism name (human, mouse, rat)")
    download_pt_parser.add_argument("--output", type=str, help="Directory path to save the downloaded data")
    
    return parser.parse_args()


def main():
    """Main function."""
    args = parse_arguments()
    
    if args.command == "list-sources":
        OmnipathExplorer.list_interaction_sources()
    
    elif args.command == "list-organisms":
        OmnipathExplorer.list_available_organisms()
    
    elif args.command == "download":
        include_datasets = args.include.split(",") if args.include else None
        exclude_datasets = args.exclude.split(",") if args.exclude else None
        OmnipathDownloader.download_interactions(
            organism=args.organism,
            dataset_type="all",
            include_datasets=include_datasets,
            exclude_datasets=exclude_datasets,
            output_file=args.output
        )
    
    elif args.command == "download-dataset":
        OmnipathDownloader.download_interactions(
            organism=args.organism,
            dataset_type="specific",
            specific_dataset=args.dataset,
            output_file=args.output
        )
    
    elif args.command == "download-all":
        OmnipathDownloader.download_interactions(
            organism=args.organism,
            dataset_type="all",
            output_file=args.output
        )
    
    elif args.command == "download-pt":
        OmnipathDownloader.download_interactions(
            organism=args.organism,
            dataset_type="post_translational",
            output_file=args.output
        )
    
    else:
        print("Please specify a command. Use --help for more information.")
        sys.exit(1)


if __name__ == "__main__":
    main()
