FROM condaforge/miniforge3:latest

COPY environment.yml requirements.txt /tmp/conda-tmp/
RUN conda env create -f /tmp/conda-tmp/environment.yml && conda clean -afy && rm -rf /tmp/conda-tmp

# Make the env the default without needing `conda activate` everywhere
SHELL ["conda", "run", "--no-capture-output", "-n", "cox2-db", "/bin/bash", "-c"]
ENV PATH=/opt/conda/envs/cox2-db/bin:$PATH
