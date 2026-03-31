# LISFLOOD-FP Docker Build
# Multi-stage build for optimized image size

# ---- Build stage ----
FROM nvidia/cuda:12.4.0-devel-ubuntu22.04 AS builder

# Install dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    libnuma-dev \
    libnetcdf-dev \
    git \
    && rm -rf /var/lib/apt/lists/*

# Set up work directory
WORKDIR /app

# Copy source code
COPY . .

# Create build directory and compile
RUN mkdir -p build && cd build     && cmake .. -D_CONFIG=../config/local     && make -j $(nproc) VERBOSE=1 && ls -l

# ---- Runtime stage ----
FROM nvidia/cuda:12.4.0-runtime-ubuntu22.04

# Install runtime dependencies
RUN apt-get update && apt-get install -y \
    libnuma-dev \
    libnetcdf-dev \
    && rm -rf /var/lib/apt/lists/*

# Create app directory
WORKDIR /app

# Copy only the built executable and necessary files from builder
COPY --from=builder /app/build/lisflood /app/
COPY examples/ /app/examples/
COPY docs/ /app/docs/

# Create directories for mounting input/output
RUN mkdir -p /data/input /data/output

# Set environment variables
ENV LD_LIBRARY_PATH=/usr/local/cuda/lib64:$LD_LIBRARY_PATH

# Entrypoint script
COPY docker-entrypoint.sh /usr/local/bin/
RUN chmod +x /usr/local/bin/docker-entrypoint.sh
ENTRYPOINT ["docker-entrypoint.sh"]

# Default command
CMD ["/app/lisflood", "-h"]