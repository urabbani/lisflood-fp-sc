#!/bin/bash
set -e

# If first arg starts with hyphen, prepend lisflood to it
if [ "${1:0:1}" = "-" ]; then
    set -- /app/lisflood "$@"
fi

# Execute command
exec "$@"