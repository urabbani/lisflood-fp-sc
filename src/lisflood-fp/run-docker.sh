#!/bin/bash
# LISFLOOD-FP Docker Run Script

# Default variables
INPUT_DIR="$(pwd)/data/input"
OUTPUT_DIR="$(pwd)/data/output"
USE_GPU=false

# Help function
function show_help {
    echo "LISFLOOD-FP Docker Run Script"
    echo "Usage: ./run-docker.sh [OPTIONS]"
    echo ""
    echo "Options:"
    echo "  -i, --input DIR    Input directory path (default: ./data/input)"
    echo "  -o, --output DIR   Output directory path (default: ./data/output)"
    echo "  -g, --gpu          Enable GPU acceleration"
    echo "  -p, --paramfile    Parameter file name (relative to input directory)"
    echo "  -h, --help         Show this help message"
    echo ""
    echo "Example:"
    echo "  ./run-docker.sh -i ./my_input -o ./my_output -g -p simulation.par"
}

# Parse arguments
POSITIONAL=()
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        -i|--input)
        INPUT_DIR="$2"
        shift
        shift
        ;;
        -o|--output)
        OUTPUT_DIR="$2"
        shift
        shift
        ;;
        -g|--gpu)
        USE_GPU=true
        shift
        ;;
        -p|--paramfile)
        PARAM_FILE="$2"
        shift
        shift
        ;;
        -h|--help)
        show_help
        exit 0
        ;;
        *)
        POSITIONAL+=("$1")
        shift
        ;;
    esac
done
set -- "${POSITIONAL[@]}"

# Check if input directory exists
if [ ! -d "$INPUT_DIR" ]; then
    echo "Error: Input directory does not exist: $INPUT_DIR"
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Set GPU options
if [ "$USE_GPU" = true ]; then
    GPU_OPTS="--gpus all"
    LISFLOOD_OPTS="-cuda"
else
    GPU_OPTS=""
    LISFLOOD_OPTS=""
fi

# Build the docker command
DOCKER_CMD="docker run --rm -i \
    -v $INPUT_DIR:/data/input \
    -v $OUTPUT_DIR:/data/output \
    $GPU_OPTS \
    lisflood-fp"

# Add parameter file if specified
if [ ! -z "$PARAM_FILE" ]; then
    DOCKER_CMD="$DOCKER_CMD /app/lisflood $LISFLOOD_OPTS /data/input/$PARAM_FILE"
fi

# Execute docker command
echo "Running LISFLOOD-FP with the following configuration:"
echo "  Input directory: $INPUT_DIR"
echo "  Output directory: $OUTPUT_DIR"
echo "  GPU enabled: $USE_GPU"
if [ ! -z "$PARAM_FILE" ]; then
    echo "  Parameter file: $PARAM_FILE"
fi
echo ""
echo "Docker command: $DOCKER_CMD"
echo ""

# Execute the command
eval $DOCKER_CMD