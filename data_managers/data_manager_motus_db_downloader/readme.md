# Serve locally

* Install conda/mamba

```
mamba create -n planemo-env python=3.7
mamba activate planemo-env
pip install -U planemo
```

## Run DM togehter with tool

```
cd <DM-path>/motus-DM
planemo serve data_manager/motus_db_fetcher.xml <motus_profiler PATH>/motus_profiler.xml --biocontainers --galaxy_root ~/git/galaxy
```

## Check if tool and DM work together

In Galaxy go to:

* Admin
* Local Data
* Check if DM is in **Installed Data Managers**
* Click it
* Run the tool with Database type 3.1.0
* Admin
* Data Tables
* motus_db_versioned
* Check if new table was made
* Go to the tool
* Check if new table can be found by: `A pre-installed mOTUs database`