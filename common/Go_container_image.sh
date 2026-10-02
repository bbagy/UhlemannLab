#!/usr/bin/env bash
set -euo pipefail

usage(){
  cat <<'EOF'
Usage:
  bash common/Go_container_image.sh export IMAGE OUTPUT.docker.tar
  bash common/Go_container_image.sh build INPUT.docker.tar OUTPUT.sif

export: Run on the Docker workstation. Save the image and its ID/platform.
build:  Run where Apptainer is available. Convert the archive to a SIF file.
Transfer files separately using rsync -avP (see common/README.md).
EOF
}

if [ "${1:-}" = --help ] || [ "${1:-}" = -h ]; then usage; exit 0; fi
[ "$#" -eq 3 ] || { usage >&2; exit 2; }
ACTION="$1"
INPUT="$2"
OUTPUT="$3"
case "$ACTION" in
  export|build) ;;
  *) usage >&2; exit 2 ;;
esac
[ ! -e "$OUTPUT" ] || { echo "Output already exists: $OUTPUT" >&2; exit 1; }

case "$ACTION" in
  export)
    command -v docker >/dev/null 2>&1 || { echo "Docker is not available" >&2; exit 1; }
    IMAGE_INFO="$(docker image inspect --format 'id={{.Id}} platform={{.Os}}/{{.Architecture}} tags={{json .RepoTags}}' "$INPUT")"
    mkdir -p "$(dirname "$OUTPUT")"
    docker save --output "$OUTPUT" "$INPUT"
    printf '%s\n' "$IMAGE_INFO" > "$OUTPUT.image.txt"
    echo "Archive: $OUTPUT"
    echo "Image ID/platform: $OUTPUT.image.txt"
    ;;
  build)
    command -v apptainer >/dev/null 2>&1 || { echo "Apptainer is not available; load the HPC module first" >&2; exit 1; }
    [ -f "$INPUT" ] || { echo "Archive not found: $INPUT" >&2; exit 1; }
    ARCHIVE="$(cd "$(dirname "$INPUT")" && pwd -P)/$(basename "$INPUT")"
    mkdir -p "$(dirname "$OUTPUT")"
    apptainer build "$OUTPUT" "docker-archive:$ARCHIVE"
    echo "SIF: $OUTPUT"
    ;;
esac
