# MetaDbGen

**Meta**barcoding **D**ata**b**ase **Gen**erator (MetaDbGen), an automated workflow for the curation of high-quality and reproducible reference sequence databases for metabarcoding applications.

A manuscript is currently in preparation.

## Installation

This application uses a Conda environment for dependency management. Please [install Conda](https://docs.conda.io/projects/conda/en/stable/index.html) if you don't have it already on your device.

### Supported Platforms

Linux and macOS are supported on x64 architecture; macOS is also supported on ARM (aarch64). While we try to support Windows as well, this is no longer possible because the dependence [HMMER](https://hmmer.org/) doesn't support Windows. However, you can still run MetaDbGen inside [WSL2](https://learn.microsoft.com/windows/wsl/install).

### Creating a Conda Environment

After downloading the source code, open a terminal in the directory that contains the source code.

> [!IMPORTANT]  
> On **Windows**, please make sure to use a Conda terminal window. If the command `conda --help` succeeds, you are good to go.

You can create a conda environment following the example below.
Please replace `[ENV_NAME]` in the commands below with a meaningful name.

```Console
conda env create -n [ENV_NAME] --file environment.yml -y
conda activate [ENV_NAME]
```

> [!IMPORTANT]  
> Make sure to run `conda activate [ENV_NAME]` each time before using the workflow.

## Usage

The workflow is implemented using the workflow management system [Snakemake](https://snakemake.readthedocs.io/en/v9.16.3/). Snakemake takes care of efficiently running the steps of the workflow. We included Snakemake in the conda environment due to the high frequency of breaking changes. Although other versions of Snakemake than the currently supported version 9.16.3 may also work, we cannot guarantee it.

Modify our example [config.yml](./config.yml) file to configure the workflow for your barcoding locus and taxonomic group of interest.

You can execute the workflow from within your conda environment using this command:

```Console
snakemake -s cox2-db.smk --configfile config.yml --cores 8
```

The flag `--cores` sets the maximum number of cores snakemake will use. Please see the [Snakemake documentation](https://snakemake.readthedocs.io/en/v9.16.3/executing/cli.html) for additional command line parameters.
