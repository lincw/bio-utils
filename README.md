# Bio-Utils

A collection of portable and reusable data analysis utilities designed for bioinformatics research, with a focus on protein interaction network analysis and transcriptomics integration.

## Overview

This repository contains modular analysis tools that can be easily integrated into different projects without hard-coded file paths or project-specific dependencies. All tools are designed to be path-agnostic and accept data inputs through parameters.

## Project Structure

## Directory Structure

- [python/](python) - Python utilities and modules
- [R/](R) - R functions and scripts
- [utils/](utils) - Cross-platform utilities and configurations
- [docs/](docs) - Documentation and integration guides
- [README.md](README.md) - This file
- [LICENSE](LICENSE) - License information

## Design Principles

- **Path Agnostic**: No hard-coded file paths; all inputs specified via parameters
- **Modular**: Each tool serves a specific analytical purpose
- **Portable**: Easy to copy and integrate into different projects
- **Well-Documented**: Clear function signatures and usage examples
- **Language Flexible**: Support for both R and Python workflows

## Target Applications

- Protein interaction network analysis
- Transcriptomics data integration
- Network topology analysis
- Functional enrichment analysis
- Multi-omics data integration
- Statistical analysis and visualization

## Integration

For instructions on how to integrate bio-utils into your projects using Git submodules or sparse-checkout, see the [integration documentation](docs/integration_guide.md).

Quick start:
```bash
# Add as submodule
git submodule add https://github.com/username/bio-utils.git tools

# Or see docs/quick_reference.md for more options
```

## Contributing

When adding new tools:
1. Ensure no hard-coded paths
2. Include comprehensive documentation
3. Provide usage examples
4. Follow consistent naming conventions
5. Test with different input formats

## Date Format

All date outputs use YYYY-MM-DD format consistently across tools.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

**Note**: This repository focuses on reusable analytical components rather than complete analysis pipelines. For project-specific workflows, create separate repositories that import these tools.
