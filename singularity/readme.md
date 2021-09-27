# Information to run container via Singularity on HPC

## You could build the provided Docker container (defined in the `../docker/` directory yourself in Singularity, using ``--remote`` if you have no sudo access like below, and with `--sandbox` or `--writable` to allow interaction with the filesystem,
```bash
singularity build --remote --sandbox ccpp-scm-v5-tutorial docker://dtcenter/ccpp-scm:v5.0.0-tutorial
```
## or you could build it from the provided `Singularity.def` file in this directory, again with `--sandbox` or `--writable` to allow interaction with the filesystem,  
```bash
singularity build --remote --sandbox ccpp-scm-v5-tutorial Singularity.def
```
and that is it --- though you need to make sure you have set up your tokens in order to use the ``--remote`` option. (Ideally, the DTCenter should make their own Singularity cloud account and add a functioning image to theirs, following the build procedure. Another option: You can build directly on the webpage: https://cloud.sylabs.io/builder.) To run the container through an interactive shell, use the command below with `--writable` to allow interaction with the filesystem, otherwise it will run read-only.
```bash
singularity shell --writable ccpp-scm-v5-tutorial
```

Finally, to run the [tutorial](https://dtcenter.org/ccpp-scm-online-tutorial), you could something like this:
- `cd /comsoftware/ccpp-scm/scm/bin/`
- `./multi_run_scm.py -c twpice`

(Note: The environmental variables seem to be wrong in the tutorial, so working with absolute paths may be preferable at this point.)

## Troubleshooting:
The `Singularity.def` in this directory only does the following steps, so try going through them step-by-step to figure out what's causing your problem...

### The `common-community-container:gnu9` container: 
- `singularity build --remote --sandbox test docker://dtcenter/common-community-container:gnu9`

- `singularity shell --writable test`

### The steps inside (unsure if all are needed, but the build seems to result in failure without some of them...): 
- `cd /comsoftware && git clone --recursive -b v5.0.0 https://github.com/NCAR/ccpp-scm`

- `cd /comsoftware/ccpp-scm/ && . contrib/get_thompson_tables.sh` 

- `cd /comsoftware/ccpp-scm/ && . contrib/get_mg_inccn_data.sh`

- `cd /comsoftware/ccpp-scm/scm/ && . etc/CENTOS_docker_setup.sh`
- `bacio_ROOT=/comsoftware/ccpp-scm/nceplibs` 
- `sp_ROOT=/comsoftware/ccpp-scm/nceplibs`
- `w3nco_ROOT=/comsoftware/ccpp-scm/nceplibs`

- `cd /comsoftware/ccpp-scm/scm && mkdir bin && cd bin`

- `export bacio_ROOT=/comsoftware/ccpp-scm/nceplibs`
- `export sp_ROOT=/comsoftware/ccpp-scm/nceplibs`
- `export w3nco_ROOT=/comsoftware/ccpp-scm/nceplibs`

- `ln -s /usr/bin/python3 /usr/local/bin/python`

- `cd /comsoftware/ccpp-scm/scm/bin && cmake ../src && make -j4`

### Running the [tutorial](https://dtcenter.org/ccpp-scm-online-tutorial) to test if it actually works:
- `cd /comsoftware/ccpp-scm/scm/bin/`
- `./multi_run_scm.py -c twpice` 


## Your next step is to figure out how to bind your container to your filesystem and make _scientific_ use of it.  