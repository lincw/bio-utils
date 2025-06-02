# Integration Guide

This guide explains how to integrate the bio-utils toolkit into your specific research projects using Git submodules and sparse-checkout features.

## Method 1: Using Git Submodule

### Basic Submodule Integration

Add the entire bio-utils toolkit as a submodule to your project:

```bash
# Navigate to your project directory
cd /path/to/your/project

# Add bio-utils as a submodule
git submodule add https://github.com/username/bio-utils.git tools

# Initialize and update the submodule
git submodule update --init --recursive

# Commit the submodule addition
git add .gitmodules tools
git commit -m "Add bio-utils as submodule"
```

### Updating Submodule to Latest Version

```bash
# Update submodule to latest commit
git submodule update --remote tools

# Commit the updated submodule reference
git add tools
git commit -m "Update bio-utils submodule to latest version"
```

### Working with Submodules in Team Projects

When team members clone your project:

```bash
# Clone project with submodules
git clone --recurse-submodules https://github.com/username/your-project.git

# Or if already cloned, initialize submodules
git submodule update --init --recursive
```

## Method 2: Using Git Sparse-Checkout

### Sparse-Checkout for Specific Components

If you only need specific parts of bio-utils, use sparse-checkout:

```bash
# Clone without checking out files
git clone --no-checkout https://github.com/username/bio-utils.git tools
cd tools

# Enable sparse-checkout
git config core.sparseCheckout true

# Create sparse-checkout specification
echo "python/network_analysis/" > .git/info/sparse-checkout
echo "R/statistics/" >> .git/info/sparse-checkout
echo "utils/" >> .git/info/sparse-checkout
echo "docs/" >> .git/info/sparse-checkout
echo "README.md" >> .git/info/sparse-checkout
echo "LICENSE" >> .git/info/sparse-checkout

# Checkout only specified directories/files
git checkout
```

### Common Sparse-Checkout Patterns

**For Python-only projects:**
```bash
echo "python/" > .git/info/sparse-checkout
echo "utils/" >> .git/info/sparse-checkout
echo "docs/" >> .git/info/sparse-checkout
echo "README.md" >> .git/info/sparse-checkout
git checkout
```

**For R-only projects:**
```bash
echo "R/" > .git/info/sparse-checkout
echo "utils/" >> .git/info/sparse-checkout
echo "docs/" >> .git/info/sparse-checkout
echo "README.md" >> .git/info/sparse-checkout
git checkout
```

**For specific analysis types:**
```bash
echo "python/network_analysis/" > .git/info/sparse-checkout
echo "R/bioinformatics/" >> .git/info/sparse-checkout
echo "utils/config_templates/" >> .git/info/sparse-checkout
echo "docs/user_guides/" >> .git/info/sparse-checkout
git checkout
```

## Method 3: Combining Submodule with Sparse-Checkout

For projects requiring only specific components as a submodule:

```bash
# Add submodule
git submodule add https://github.com/username/bio-utils.git tools
cd tools

# Configure sparse-checkout within submodule
git config core.sparseCheckout true
echo "python/network_analysis/" > .git/info/sparse-checkout
echo "docs/api_reference/" >> .git/info/sparse-checkout
echo "README.md" >> .git/info/sparse-checkout

# Apply sparse-checkout
git read-tree -m -u HEAD

# Return to parent project and commit
cd ..
git add .gitmodules tools
git commit -m "Add bio-utils submodule with sparse-checkout configuration"
```

## Project Structure Examples

### Full Integration
```
your-project/
├── data/
├── scripts/
│   ├── analysis.py
│   └── process_data.R
├── tools/                  # Complete bio-utils submodule
│   ├── python/
│   ├── R/
│   ├── utils/
│   ├── docs/
│   └── README.md
├── results/
└── README.md
```

### Sparse Integration (Python-only)
```
your-project/
├── data/
├── scripts/
│   └── analysis.py
├── tools/                  # Sparse bio-utils submodule
│   ├── python/
│   ├── utils/
│   ├── docs/
│   └── README.md
├── results/
└── README.md
```

## Code Integration Examples

### Python Integration

```python
# Add tools to Python path
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), 'tools', 'python'))

# Import specific modules
from network_analysis import protein_networks
from data_processing import transcriptomics_utils

# Use the tools
network = protein_networks.load_ppi_network(data_file)
results = transcriptomics_utils.integrate_expression_data(expression_data, network)
```

### R Integration

```r
# Source required functions
source(file.path("tools", "R", "bioinformatics", "network_functions.R"))
source(file.path("tools", "R", "statistics", "enrichment_analysis.R"))

# Use the functions
network_data <- load_protein_network(network_file)
enrichment_results <- perform_go_enrichment(gene_list, background)
```

## Best Practices

### Version Management
- Pin submodules to specific commits for reproducible research
- Use semantic versioning tags when available
- Document which version of bio-utils was used in your analysis

### Working with Git Tags in Submodules

To use a specific tagged version of bio-utils in your project:

```bash
# Navigate to the submodule directory
cd tools/bio-utils

# Fetch all tags from remote
git fetch --tags

# List available tags
git tag -l

# Checkout a specific tag (e.g., v1.0.0)
git checkout tags/v1.0.0

# Or in a single command
git -C tools/bio-utils checkout $(git -C tools/bio-utils describe --tags $(git -C tools/bio-utils rev-list --tags --max-count=1))

# Return to parent directory and commit the submodule state
cd ../..
git add tools/bio-utils
git commit -m "Pin bio-utils to version v1.0.0"
```

To update to the latest tag:

```bash
# Navigate to the submodule
cd tools/bio-utils

# Fetch all tags and updates
git fetch --all --tags

# Get the latest tag name
LATEST_TAG=$(git describe --tags $(git rev-list --tags --max-count=1))

# Checkout the latest tag
git checkout $LATEST_TAG

# Return to parent directory and commit the update
cd ../..
git add tools/bio-utils
git commit -m "Update bio-utils to latest version $LATEST_TAG"
```

### Path Management
- Use relative paths when referencing submodule contents
- Set up environment variables for commonly used tool paths
- Create wrapper scripts to simplify tool access

### Documentation
- Document which bio-utils components your project uses
- Include submodule update instructions in project README
- Maintain compatibility notes when updating bio-utils version

## Troubleshooting

### Common Issues

**Submodule not initialized:**
```bash
git submodule update --init --recursive
```

**Sparse-checkout not working:**
```bash
# Verify sparse-checkout is enabled
git config core.sparseCheckout

# Check sparse-checkout file
cat .git/info/sparse-checkout

# Re-apply sparse-checkout
git read-tree -m -u HEAD
```

**Submodule pointing to wrong commit:**
```bash
# Update to latest
git submodule update --remote

# Or checkout specific commit
cd tools
git checkout <specific-commit-hash>
cd ..
git add tools
git commit -m "Update tools to specific version"
```

### Team Collaboration Tips

1. Always commit submodule updates explicitly
2. Use `.gitmodules` to specify branch tracking if needed
3. Document sparse-checkout patterns in project documentation
4. Consider using Git hooks to automate submodule updates
