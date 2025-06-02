# Quick Reference

Fast commands for integrating bio-utils into your projects.

## Quick Start Commands

### Full Integration
```bash
# Add complete bio-utils as submodule
git submodule add https://github.com/username/bio-utils.git tools
git submodule update --init --recursive
```

### Python-Only Integration
```bash
git clone --no-checkout https://github.com/username/bio-utils.git tools
cd tools
git config core.sparseCheckout true
echo "python/" > .git/info/sparse-checkout
echo "utils/" >> .git/info/sparse-checkout
echo "README.md" >> .git/info/sparse-checkout
git checkout
```

### R-Only Integration
```bash
git clone --no-checkout https://github.com/username/bio-utils.git tools
cd tools
git config core.sparseCheckout true
echo "R/" > .git/info/sparse-checkout
echo "utils/" >> .git/info/sparse-checkout
echo "README.md" >> .git/info/sparse-checkout
git checkout
```

### Network Analysis Only
```bash
git clone --no-checkout https://github.com/username/bio-utils.git tools
cd tools
git config core.sparseCheckout true
echo "python/network_analysis/" > .git/info/sparse-checkout
echo "R/network_analysis/" >> .git/info/sparse-checkout
echo "docs/" >> .git/info/sparse-checkout
git checkout
```

## Common Import Patterns

### Python
```python
import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), 'tools', 'python'))
from network_analysis import protein_networks
```

### R
```r
source(file.path("tools", "R", "network_analysis", "functions.R"))
```

## Update Commands

```bash
# Update submodule
git submodule update --remote tools

# Update and commit
git submodule update --remote tools
git add tools
git commit -m "Update bio-utils to latest version"
```

For detailed instructions, see [integration_guide.md](integration_guide.md).
