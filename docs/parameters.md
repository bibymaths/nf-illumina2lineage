# Pipeline Parameters

This section lists the configurable parameters used in the pipeline. These can be adjusted using a custom configuration file or passed directly via the command line.

---

## Core Parameters

| Parameter         | Default               | Description |
|------------------|-----------------------|-------------|
| `reads`          | `*.fastq.gz`          | Glob pattern to match input FASTQ files |
| `reference`      | `reference.fasta`     | Reference genome in FASTA format |
| `threads`        | `8`                   | Number of threads to use for multithreaded tools |

---

## Optional Profiles

You can use execution profiles to control behavior, e.g.:
```bash
nextflow run main.nf -profile docker
```

### Available Profiles

- **`docker`**: Runs processes inside Docker containers using the specified image.
- **`singularity`**: Runs in Singularity containers (if configured).
- **`local`**: Default execution on local machine without containers.

---

## Custom Configuration
To override parameters, use `-params-file`:
```bash
nextflow run main.nf -params-file custom_params.json
```

Example `custom_params.json`:
```json
{
  "reads": "data/*.fq.gz",
  "reference": "data/NC_045512.2.fasta",
  "threads": 16
}
```

> 💡 All parameters can also be set directly via `-params` flag or in `nextflow.config`.
