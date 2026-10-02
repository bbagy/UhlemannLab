# Docker images and Apptainer on HPC

Keep Dockerfiles, workflows, and execution/conversion scripts in Git. Store
Docker images, exported archives, and Apptainer SIF files separately on the
workstation or HPC. Replace `USER@HPC_HOST` below with your HPC login address.

## 1. Locate the existing Docker images

Run these commands on the workstation where the images were built:

```bash
docker image ls
docker context show
docker info --format '{{.DockerRootDir}}'
```

Use the REPOSITORY and TAG from `docker image ls` to select an image.
DockerRootDir is the storage directory of the active Docker daemon. On Linux
it is commonly `/var/lib/docker`; with Docker Desktop it is inside its VM.
If the selected Docker context is remote, the directory is on that remote host.
Export individual images with `docker save` rather than copying Docker's
internal storage directories.

These are the image names currently used by the launchers. Check the actual
workstation image list before exporting; the images themselves are not in Git.

| Pipeline | Default Docker image |
|---|---|
| shortWGS | `shortwgs` (implicit `latest` tag) |
| longWGS | `longwgs` (implicit `latest` tag) |
| RNake | `rnake:1.0` |
| PFsnake | `pf-snake:1.0` |
| KBracken | `kbracken:1.0` |
| Humann | `humann:1.0` |
| MAGs QC | `mags-qc:1.0` |
| MAGs Assembly | `mags-assembly:1.0` |
| MAGs Annotation | `mags-annotation:1.0` |
| TnSeq_ONT | `flank-pipeline:1.0` (README example) |

## 2. Export a Docker image

On the workstation, change to the repository root and run:

```bash
bash common/Go_container_image.sh export shortwgs \
  container-images/shortwgs.docker.tar
```

The script creates the output directory and saves the archive plus an
`.image.txt` sidecar containing the source image ID, OS/CPU architecture, and
tags. It refuses to overwrite an existing output file. For another pipeline,
change the image name and output filename.

Equivalent manual commands:

```bash
mkdir -p container-images
docker save --output container-images/shortwgs.docker.tar shortwgs
```

Use `docker save`, not `docker export`: the latter exports a container's
filesystem rather than the image with its configuration. The image CPU
architecture must match the target HPC.

## 3. Transfer the archive to HPC

Run on the workstation. `rsync -avP` displays progress and retains partial
files so an interrupted transfer can be retried.

```bash
ssh USER@HPC_HOST 'mkdir -p ~/containers'
rsync -avP container-images/shortwgs.docker.tar \
  container-images/shortwgs.docker.tar.image.txt \
  USER@HPC_HOST:~/containers/
```

Prepare the code on HPC as well. For a first-time setup:

```bash
# Run on HPC.
git clone https://github.com/bbagy/UhlemannLab.git
cd UhlemannLab
```

For an existing checkout, check for local changes and update with
`git pull --ff-only`. If your institution provides a transfer node or shared
storage location, use its recommended host and path instead.

## 4. Convert the archive to SIF

Run on a machine with Apptainer installed. On HPC, load the appropriate module
if needed, and perform conversion on a node or within a job allocation where
your cluster permits it.

```bash
module load apptainer  # Site-specific; omit if Apptainer is already available.
cd ~/UhlemannLab
bash common/Go_container_image.sh build \
  "$HOME/containers/shortwgs.docker.tar" \
  "$HOME/containers/shortwgs.sif"
```

Equivalent manual command:

```bash
apptainer build "$HOME/containers/shortwgs.sif" \
  "docker-archive:$HOME/containers/shortwgs.docker.tar"
```

**Use `docker-archive:` followed by a file path, not `docker-archive://`.**
See the [official Apptainer archive conversion guide](https://apptainer.org/docs/user/main/docker_and_oci.html#containers-in-docker-archive-files).
This method does not require a Docker daemon or registry download on HPC.
Allow enough space for conversion temporary files and the final SIF. If your
cluster does not permit conversion, build the SIF on another Linux machine
with Apptainer, then transfer the SIF using the same rsync procedure.
`docker load` can restore the archive into Docker on a workstation, but it is
not needed for Apptainer conversion on HPC.

## 5. Check the tools and run the pipeline

```bash
apptainer exec --cleanenv "$HOME/containers/shortwgs.sif" snakemake --version
apptainer exec --cleanenv "$HOME/containers/shortwgs.sif" fastp --version

# Replace input/DB paths with real HPC paths and start with a dry-run.
bash ~/UhlemannLab/shortWGS/Go_shortWGS.sh \
  --container apptainer \
  --container-image "$HOME/containers/shortwgs.sif" \
  -i /shared/fastq -o shortwgs_out \
  -d /shared/db/wgs -k /shared/db/kraken2 -r /shared/GoWGS \
  -c 8 -n
```

All input paths must exist **on HPC**. Transferring the container image does
not transfer FASTQ files, external databases, or the GoWGS directory; prepare
those separately in HPC storage. After the dry-run passes, remove `-n` and
execute within an appropriate HPC allocation. Match the core count to your
allocated resources.

## Versions and storage

Give validated images versioned filenames, such as `shortwgs_v1.sif`, and
retain the source `.image.txt`. You can record the SIF checksum with
`sha256sum shortwgs_v1.sif`. Recording the Git commit, source Docker image ID,
and final SIF checksum identifies the environment used for an analysis.
Images have not yet been built or transferred as part of this change; these
are preparation and validation instructions.

`.gitignore` excludes `container-images/`, `*.docker.tar`, the corresponding
image metadata, and `*.sif`. Deploy `common/container.sh` along with the
pipeline launchers; cloning the complete repository includes it automatically.
For runtime selection, see the [main README](../README.md#docker--apptainer-selection).
