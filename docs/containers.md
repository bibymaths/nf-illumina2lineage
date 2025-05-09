# Container Support

This pipeline supports containerized execution using **Docker** and **Singularity**, which ensures reproducibility and simplifies environment setup.

---

## Docker

To run the pipeline with Docker:
```bash
nextflow run main.nf -profile docker
```

Make sure Docker is installed and the current user has permission to run Docker commands (e.g., add user to `docker` group).

### Dockerfile
A `Dockerfile` is provided to build a custom container image:
```bash
docker build -t nf-illumina2lineage .
```

Use your image with:
```bash
nextflow run main.nf -with-docker nf-illumina2lineage
```

---

## Singularity

To run with Singularity (useful on HPCs):
```bash
nextflow run main.nf -profile singularity
```

Make sure `singularity` is installed (or `apptainer` on newer systems).

To convert the Docker image:
```bash
singularity build nf-illumina2lineage.sif docker-daemon://nf-illumina2lineage:latest
```

---

## Profile Configuration (Optional)
Customize container settings in `nextflow.config`:
```groovy
process.container = 'nf-illumina2lineage'
docker.enabled = true
```

> 📦 Using containers is recommended for dependency consistency and reproducibility across systems.
