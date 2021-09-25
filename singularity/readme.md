# Information to run container via Singularity on HPC
## If you just want to pull it from the [Singularity hub](https://cloud.sylabs.io/library/ngam/remote-builds/rb-614f66b3a267b3f6a3832fbe),
```bash
singularity pull library://ngam/remote-builds/rb-614f5842a267b3f6a3832fbb:sha256.3c3b09f362ec69c56776c556d31fbb76f1f1435b198ec4fcebe01452d709450b
```

## Or, actually, preferably you should to build it yourself, using ``--remote`` if you have no sudo access like below,
```bash
singularity build --remote --sandbox ccpp-scm-v5-tutorial Singularity.def
```
and that is it --- though you need to make sure you have set up your tokens in order to use the ``--remote`` option. (Ideally, the DTCenter should make their own Singularity cloud account and add a functioning image to theirs, following the build procedure. Another option: one can build directly on the webpage: https://cloud.sylabs.io/builder.) To run the container, use the command below.
```bash
singularity shell --writable ccpp-scm-v5-tutorial
```

Finally, to run [tutorial](https://dtcenter.org/ccpp-scm-online-tutorial), you could something like this:
- `cd /comsoftware/ccpp-scm/scm/bin/`
- `./multi_run_scm.py -c twpice`  

## TROUBLESHOOTING:
The `Singularity.def` in this directory only does the following steps, so try going through them step-by-step to figure out what's causing your problem...

### The `common-community-container:gnu9` container: 
- `singularity build --remote --sandbox test docker://dtcenter/common-community-container:gnu9`

- `singularity shell --writable test`

### The steps inside: 
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