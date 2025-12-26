# OCI Development

Here I'll write how I created an OCI bundle using Singularity to help with
development of the client. First, notice the [config.json](config.json)
in the present working directory - it's a general configuration for OCI
runtime specification (version 1.0) that we can put into a bundle (a folder
that will serve as the root of a container). Here is how I did that.

First, we are developing with the first release of Singularity that supports
OCI:

```bash
$ singularity --version
singularity version 3.1.0-rc2.28.ga72e427
```

Next, we are going to create a bundle. A bundle is a folder that we will
treat as the root of our container filesystem. The easiest way to do
this is to dump a filesystem there from another container.

```bash
$ singularity build --sandbox /tmp/bundle docker://ubuntu:18.04
$ cp config.json /tmp/bundle
```

The purpose of the build is only to dump a complete operating system into the
bundle folder. The configuration file then is to conform to the Oci
Runtime specification.

We can then test interaction with the OCI client of Singularity python
by providing the bundle directory at /tmp/bundle.
