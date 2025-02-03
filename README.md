# WES SNV Discovery Workflow

This repository contains a workflow used in-house for discovering single nucleotide variants (SNVs) from whole exome sequencing (WES) data.

## Table of Contents
- [Introduction](#introduction)
- [Installation](#installation)
- [Usage](#usage)
- [Workflow Overview](#workflow-overview)
- [Contributing](#contributing)
- [License](#license)

## Introduction
This workflow is designed to identify SNVs from WES data generated in-house. It leverages various bioinformatics tools packaged in GATK Suite and follows GATK Best Practices. 

## Installation
Activate the `gatk` Mamba environment in this AMI before running the Bash workflow.

## Usage
To execute the workflow, just run the following Bash script.
```
bash WES-DNA-SNV-CRMY.sh
```

## Workflow Overview
The workflow consists of the following steps:
1. **Pre-Processing**: Assess the quality of raw sequencing data.
2. **Alignment**: Align reads to the reference genome.
3. **Variant Calling**: Identify SNVs from the aligned reads.
4. **Annotation**: Annotate the identified variants.

## Contributing
Contributions are welcome! Please submit a pull request or open an issue to discuss your ideas.

## License
This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.