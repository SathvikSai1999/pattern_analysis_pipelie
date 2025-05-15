#!/bin/bash

# Remove old log files
echo "Cleaning up old log files..."
rm -f size3ES_and_P_*.log
rm -f pipeline_*.log
rm -f pipeline.log

# Create a new logs directory if it doesn't exist
mkdir -p logs

echo "Cleanup completed!" 