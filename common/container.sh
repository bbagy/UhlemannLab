#!/usr/bin/env bash
# Shared runtime selection for the UhlemannLab launchers (Bash 3 compatible).
CONTAINER_RUNTIME=docker
CONTAINER_IMAGE=""
CONTAINER_ARGS=()
while [ "$#" -gt 0 ]; do
  case "$1" in
    --container|--container-image)
      if [ "$#" -lt 2 ] || [ -z "$2" ]; then
        echo "[container] $1 requires a value" >&2
        exit 2
      fi
      if [ "$1" = --container ]; then CONTAINER_RUNTIME="$2"; else CONTAINER_IMAGE="$2"; fi
      shift 2 ;;
    --container=*) CONTAINER_RUNTIME="${1#*=}"; shift ;;
    --container-image=*)
      CONTAINER_IMAGE="${1#*=}"
      [ -n "$CONTAINER_IMAGE" ] || { echo "[container] --container-image requires a value" >&2; exit 2; }
      shift ;;
    --) CONTAINER_ARGS+=("$@"); break ;;
    *) CONTAINER_ARGS+=("$1"); shift ;;
  esac
done
case "$CONTAINER_RUNTIME" in
  docker|apptainer) ;;
  *) echo "[container] --container must be docker or apptainer" >&2; exit 2 ;;
esac

container_require_image(){
  local image="${CONTAINER_IMAGE:-$1}"
  command -v "$CONTAINER_RUNTIME" >/dev/null 2>&1 || {
    echo "[container] $CONTAINER_RUNTIME is not available; install it or load the HPC module." >&2
    return 1
  }
  if [ "$CONTAINER_RUNTIME" = docker ]; then
    docker image inspect "$image" >/dev/null 2>&1 || {
      echo "[container] Docker image unavailable: $image. Build using the pipeline Dockerfile or select --container-image IMAGE." >&2
      return 1
    }
  else
    [ -f "$image" ] && [ -r "$image" ] || {
      echo "[container] SIF file unavailable: $image. Specify --container-image /path/to/image.sif (or -m where supported)." >&2
      return 1
    }
  fi
}

# Accept the launchers' existing mount/env options; preserve Docker arguments.
container_run(){
  if [ "$CONTAINER_RUNTIME" = docker ]; then
    local -a docker_options=()
    while [ "$#" -gt 0 ]; do
      case "$1" in
        --rm) docker_options+=("$1"); shift ;;
        -u|-v|-w|-e|--network) docker_options+=("$1" "$2"); shift 2 ;;
        -*) echo "[container] Unsupported launch option: $1" >&2; return 2 ;;
        *) break ;;
      esac
    done
    local image="${CONTAINER_IMAGE:-$1}"
    shift
    container_require_image "$image" || return
    docker run ${docker_options[@]+"${docker_options[@]}"} "$image" "$@"
  else
    local -a apptainer_options=(exec --cleanenv --no-home)
    while [ "$#" -gt 0 ]; do
      case "$1" in
        --rm) shift ;;
        -u) shift 2 ;; # Apptainer runs with the caller's UID/GID.
        -v) apptainer_options+=(--bind "$2"); shift 2 ;;
        -w) apptainer_options+=(--pwd "$2"); shift 2 ;;
        -e) apptainer_options+=(--env "$2"); shift 2 ;;
        --network)
          [ "$2" = host ] || { echo "[container] Only host networking is supported" >&2; return 2; }
          shift 2 ;; # Host networking is already the Apptainer default.
        -*) echo "[container] Unsupported launch option: $1" >&2; return 2 ;;
        *) break ;;
      esac
    done
    local image="${CONTAINER_IMAGE:-$1}"
    shift
    container_require_image "$image" || return
    apptainer "${apptainer_options[@]}" "$image" "$@"
  fi
}
