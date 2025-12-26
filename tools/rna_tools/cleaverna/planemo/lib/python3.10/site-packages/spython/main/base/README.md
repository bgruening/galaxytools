# Robot Generator

This folder contains a (sub-application) for a robot name generator.

## Docker
It's built on Docker Hub, so you can run as:

```
docker run vanessa/robotname
```

or build locally first, with the present working directory as this folder.

```
docker build -t vanessa/robotname .
```

```
for i in `seq 1 10`;
     do
     docker run vanessa/robotname
done
boopy-peanut-butter-7311
blank-snack-0334
hello-buttface-6320
chocolate-bicycle-9982
frigid-frito-9511
doopy-soup-7712
phat-pancake-4952
wobbly-kitty-3213
lovely-mango-1987
milky-poo-7960
```

## Singularity

To build your image:

```
sudo singularity build robotname Singularity
```

or pull from Docker Hub :)

```
singularity pull --name robotname docker://vanessa/robotname
sregistry pull docker://vanessa/robotname
```
