# Pull base image
FROM condaforge/miniforge3:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="0904ae5f2ad79893ef3b81149de1c29922ebb2f91e32ef0b541f04094fb44195"

# Copy
COPY app/ /opt/app/
COPY pset/ /opt/pset/
COPY workflow/ /opt/workflow/
WORKDIR /opt

# Environment
ENV PYTHONPATH="${PYTHONPATH}:$(pwd)"
RUN mamba env create -f /opt/workflow/envs/env.yml && \
    mamba clean --all -y

# Ports
EXPOSE 8000

# Run the command
CMD ["mamba", "run", "-n", "pset", "bash", "-x", "/opt/app/app.sh"]
