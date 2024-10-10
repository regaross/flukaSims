# Docker_tools
## Containerization

Originally, these simulations were run on a SLAC cluster, [S3DF](https://s3df.slac.stanford.edu/#/). As seems to be quite common on such scientific computing clusters, as a simple user without root privileges, you cannot install software directly onto the cluster, but must work through code written inside a container. In particular, a [Singularity](https://sylabs.io/singularity/) container. 

There was another catch, too (there always is). You cannot create a Singularity container except for on a Linux environment. So, as a silly Mac user, I made a [Docker](https://www.Docker.com) container to boot up on my Mac, within which I could compile the Singularity container with all the necessary dependencies so I could run the 3