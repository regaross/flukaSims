# Use a base image with a suitable Linux distribution
FROM ubuntu:23.10

# Set the working directory
WORKDIR /app

# Update and install necessary packages
RUN apt-get update && \
    apt-get install -y \
        software-properties-common \
        wget \
        python3.8 \
        python3-tk

# Install any other dependencies you need for your GUI application

# Install X11 and other necessary tools for GUI
RUN apt-get update && \
    apt-get install -y \
        xauth \
        x11-apps \
        x11-utils

# RUN touch /root/.Xauthority && \
#     xauth generate :0 . trusted && \
#     chmod 600 /root/.Xauthority

RUN wget -q -O - https://flair.web.cern.ch/flair/download/ubuntu/flair.gpg |  tee /etc/apt/trusted.gpg.d/flair.asc && \
 add-apt-repository 'deb [arch=all,amd64] https://cern.ch/flair/download/ubuntu/23.10 /'

RUN apt-get install -y flair
# &&  apt update &&  apt upgrade

# # Copy your application csode into the container
# COPY . /app

# Set up environment variables if needed
# ENV DISPLAY=docker.for.mac.localhost:0.0
# # Expose any necessary ports
# EXPOSE 80

# # Command to run your GUI application
# CMD ["<your_command_to_run_the_app>"]
