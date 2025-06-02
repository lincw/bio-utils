# Python Utilities

This directory contains Python modules and functions for data analysis.

## Structure

Organize Python tools by functionality:

- `network_analysis/` - Network topology and graph analysis tools
- `data_processing/` - Data cleaning and preprocessing utilities  
- `statistics/` - Statistical analysis functions
- `visualization/` - Plotting and visualization tools
- `io_utilities/` - File input/output helpers

## Coding Standards

- Use type hints for function parameters and return values
- Include docstrings with parameter descriptions
- Follow PEP 8 style guidelines
- No hard-coded file paths - accept paths as parameters
- Include error handling for common edge cases

## Example Function Template

```python
from typing import Union, List, Dict, Optional
import pandas as pd

def example_analysis_function(
    data: pd.DataFrame,
    parameter1: str,
    parameter2: Optional[int] = None,
    output_path: Optional[str] = None
) -> Dict[str, Union[float, List]]:
    """
    Brief description of what this function does.
    
    Parameters:
    -----------
    data : pd.DataFrame
        Description of input data format
    parameter1 : str
        Description of parameter1
    parameter2 : int, optional
        Description of parameter2 (default: None)
    output_path : str, optional
        Path to save results (default: None, returns results only)
        
    Returns:
    --------
    Dict[str, Union[float, List]]
        Description of return value structure
    """
    # Implementation here
    pass
```
