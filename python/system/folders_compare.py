# -*- coding: utf-8 -*-
# script: folders_compare.py
# date_created: 2025-06-07T13:40:10+02:00
# date_modified: 2025-06-07T13:40:13+02:00

"""
Compare two directories efficiently using filecmp
"""

import filecmp
import os
from pathlib import Path
import argparse

def compare_directories_detailed(dir1, dir2):
    """
    Compare two directories with detailed reporting
    """
    print(f"Comparing:\n  {dir1}\n  {dir2}\n")
    
    # Deep comparison of directory trees
    comparison = filecmp.dircmp(dir1, dir2)
    
    def report_comparison(dcmp, level=0):
        indent = "  " * level
        
        # Files only in dir1
        if dcmp.left_only:
            print(f"{indent}Only in {dcmp.left}:")
            for file in dcmp.left_only:
                print(f"{indent}  - {file}")
        
        # Files only in dir2
        if dcmp.right_only:
            print(f"{indent}Only in {dcmp.right}:")
            for file in dcmp.right_only:
                print(f"{indent}  + {file}")
        
        # Files that differ
        if dcmp.diff_files:
            print(f"{indent}Different files:")
            for file in dcmp.diff_files:
                print(f"{indent}  ≠ {file}")
        
        # Files that are identical
        if dcmp.same_files:
            print(f"{indent}Identical files: {len(dcmp.same_files)}")
        
        # Recursively check subdirectories
        for name, sub_dcmp in dcmp.subdirs.items():
            print(f"{indent}Subdirectory: {name}")
            report_comparison(sub_dcmp, level + 1)
    
    report_comparison(comparison)
    
    # Summary
    identical = (not comparison.left_only and 
                not comparison.right_only and 
                not comparison.diff_files and
                not comparison.funny_files)
    
    return identical

def quick_directory_comparison(dir1, dir2):
    """
    Quick boolean check if directories are identical
    """
    comparison = filecmp.dircmp(dir1, dir2)
    
    def are_dirs_equal(dcmp):
        if (dcmp.left_only or dcmp.right_only or 
            dcmp.diff_files or dcmp.funny_files):
            return False
        
        for sub_dcmp in dcmp.subdirs.values():
            if not are_dirs_equal(sub_dcmp):
                return False
        
        return True
    
    return are_dirs_equal(comparison)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Compare two directories and report their differences.',
        formatter_class=argparse.RawTextHelpFormatter,
        epilog='Example:\n  python %(prog)s /path/to/dir1 /path/to/dir2'
    )
    parser.add_argument('dir1', metavar='DIR1', type=str, help='Path to the first directory')
    parser.add_argument('dir2', metavar='DIR2', type=str, help='Path to the second directory')
    
    args = parser.parse_args()

    # Validate paths
    path1 = Path(args.dir1)
    path2 = Path(args.dir2)

    if not path1.is_dir():
        print(f"Error: The first path '{args.dir1}' is not a valid directory or does not exist.")
        exit(1)
    if not path2.is_dir():
        print(f"Error: The second path '{args.dir2}' is not a valid directory or does not exist.")
        exit(1)

    # Quick check
    print(f"Performing quick comparison between:\n  {path1}\n  {path2}\n")
    if quick_directory_comparison(str(path1), str(path2)):
        print("✓ Directories are identical (quick check).")
    else:
        print("✗ Directories differ (quick check).")
        print("\nPerforming detailed comparison:")
        compare_directories_detailed(str(path1), str(path2))