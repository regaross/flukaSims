#!/bin/bash

# Get the IP address of the host
ip=$(ifconfig en0 | grep inet | awk '$1=="inet" {print $2}')

# Allow connections to X server from host
xhost + $ip

local_dir="./"

# Run Docker container with necessary configurations
docker run -it --rm -u root --privileged --name flair_singularity_image -e DISPLAY=$ip:0 -v /tmp/.X11-unix:/tmp/.X11-unix -v $local_dir:/root flair_singularity_container
