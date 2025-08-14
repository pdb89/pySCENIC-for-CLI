# **pySCENIC-for-CLI**

### STILL IN TESTING

## A standalone Python script to run a pySCENIC pipeline from the command line with user-defined inputs:

#### usage: `pySCENIC_steps.py [-h] wdir exp_matrix f_tfs f_db_dir f_motif_path outdir`

#### positional arguments:  
  `wdir:          working directory`  
  `exp_matrix:    matrix of gene expression values exported from seurat or anndata object (rows = cells, cols = genes)`  
  `f_tfs:         TF database file`  
  `f_db_dir:      ranking database (feather files) directory`  
  `f_motif_path:  motif database file from tf2motif`  
  `outdir:        output directory`  

## CRITICAL: exp_matrix MUST be a TSV file with rows = cells and columns = genes
### User must manually create folder "output" for now


### MUST create an environment with the following dependencies (output of conda list for my environment):  
packages in environment "pySCENIC":  

Name                       Version          Build            Channel  
_libgcc_mutex                0.1              main  
_openmp_mutex                5.1              1_gnu  
aiohappyeyeballs             2.6.1            pypi_0           pypi  
aiohttp                      3.12.13          pypi_0           pypi  
aiosignal                    1.4.0            pypi_0           pypi  
anndata                      0.11.4           pypi_0           pypi  
anyio                        4.7.0            py310h06a4308_0  
arboreto                     0.1.6            pypi_0           pypi  
argon2-cffi                  21.3.0           pyhd3eb1b0_0  
argon2-cffi-bindings         21.2.0           py310h5eee18b_1  
array-api-compat             1.12.0           pypi_0           pypi  
asttokens                    3.0.0            py310h06a4308_0  
async-lru                    2.0.4            py310h06a4308_0  
async-timeout                5.0.1            pypi_0           pypi  
attrs                        24.3.0           py310h06a4308_0  
babel                        2.16.0           py310h06a4308_0  
backcall                     0.2.0            pypi_0           pypi  
beautifulsoup4               4.12.3           py310h06a4308_0  
blas                         1.0              mkl  
bleach                       6.2.0            py310h06a4308_0  
bokeh                        3.7.3            pypi_0           pypi  
boltons                      25.0.0           pypi_0           pypi  
bottleneck                   1.4.2            py310ha9d4c09_0  
brotlicffi                   1.0.9.2          py310h6a678d5_1  
bzip2                        1.0.8            h5eee18b_6  
c-ares                       1.19.1           h5eee18b_0  
ca-certificates              2025.2.25        h06a4308_0  
certifi                      2025.6.15        py310h06a4308_0  
cffi                         1.17.1           py310h1fdaa30_1  
charset-normalizer           3.3.2            pyhd3eb1b0_0  
click                        8.2.1            pypi_0           pypi  
cloudpickle                  3.1.1            pypi_0           pypi  
comm                         0.2.1            py310h06a4308_0  
contourpy                    1.3.2            pypi_0           pypi  
ctxcore                      0.2.0            pypi_0           pypi  
cycler                       0.12.1           pypi_0           pypi  
cytoolz                      1.0.1            pypi_0           pypi  
dask                         2025.5.1         pypi_0           pypi  
debugpy                      1.8.11           py310h6a678d5_0  
decorator                    5.1.1            pyhd3eb1b0_0  
defusedxml                   0.7.1            pyhd3eb1b0_0  
dill                         0.4.0            pypi_0           pypi  
distributed                  2025.5.1         pypi_0           pypi  
docopt                       0.6.2            pypi_0           pypi  
exceptiongroup               1.2.0            py310h06a4308_0  
executing                    0.8.3            pyhd3eb1b0_0  
expat                        2.7.1            h6a678d5_0  
fonttools                    4.58.5           pypi_0           pypi  
frozendict                   2.4.6            pypi_0           pypi  
frozenlist                   1.7.0            pypi_0           pypi  
fsspec                       2025.5.1         pypi_0           pypi  
h11                          0.16.0           py310h06a4308_0  
h5py                         3.14.0           py310he0d80d8_0  
hdf5                         1.14.5           h2b7332f_2  
httpcore                     1.0.9            py310h06a4308_0  
httpx                        0.28.1           py310h06a4308_0  
idna                         3.7              py310h06a4308_0  
importlib-metadata           8.5.0            py310h06a4308_0  
intel-openmp                 2023.1.0         hdb19cb5_46306  
interlap                     0.2.7            pyh9f0ad1d_0     conda-forge  
ipykernel                    6.29.5           py310h06a4308_1  
ipython                      8.12.3           pypi_0           pypi  
jedi                         0.19.2           py310h06a4308_0  
jinja2                       3.1.6            py310h06a4308_0  
joblib                       1.5.1            pypi_0           pypi  
json5                        0.9.25           py310h06a4308_0  
jsonschema                   4.23.0           py310h06a4308_0  
jsonschema-specifications    2023.7.1         py310h06a4308_0  
jupyter-lsp                  2.2.5            py310h06a4308_0  
jupyter_client               8.6.3            py310h06a4308_0  
jupyter_core                 5.8.1            pyh31011fe_0     conda-forge  
jupyter_events               0.12.0           py310h06a4308_0  
jupyter_server               2.15.0           py310h06a4308_0  
jupyter_server_terminals     0.5.3            py310h06a4308_0  
jupyterlab                   4.4.4            pyhd8ed1ab_0     conda-forge  
jupyterlab_pygments          0.3.0            py310h06a4308_0  
jupyterlab_server            2.27.3           py310h06a4308_0  
kiwisolver                   1.4.8            pypi_0           pypi  
krb5                         1.20.1           h143b758_1  
ld_impl_linux-64             2.40             h12ee557_0  
legacy-api-wrap              1.4.1            pypi_0           pypi  
libcurl                      8.12.1           hc9e6f67_0  
libedit                      3.1.20230828     h5eee18b_0  
libev                        4.33             h7f8727e_1  
libffi                       3.4.4            h6a678d5_1  
libgcc-ng                    11.2.0           h1234567_1  
libgfortran-ng               11.2.0           h00389a5_1  
libgfortran5                 11.2.0           h1234567_1  
libgomp                      11.2.0           h1234567_1  
libnghttp2                   1.57.0           h2d74bed_0  
libsodium                    1.0.18           h7b6447c_0  
libssh2                      1.11.1           h251f7ec_0  
libstdcxx-ng                 11.2.0           h1234567_1  
libuuid                      1.41.5           h5eee18b_0  
libxcb                       1.17.0           h9b100fa_0  
llvmlite                     0.44.0           py310hc1e8f15_1  
locket                       1.0.0            pypi_0           pypi  
loompy                       2.0.16           py_0             bioconda  
lz4                          4.4.4            pypi_0           pypi  
lz4-c                        1.9.4            h6a678d5_1  
markupsafe                   3.0.2            py310h5eee18b_0  
matplotlib                   3.10.3           pypi_0           pypi  
matplotlib-inline            0.1.6            py310h06a4308_0  
mistune                      3.1.2            py310h06a4308_0  
mkl                          2023.1.0         h213fc3f_46344  
mkl-service                  2.4.0            py310h5eee18b_2  
mkl_fft                      1.3.11           py310h5eee18b_0  
mkl_random                   1.2.8            py310h1128e8f_0  
mpi                          1.0              mpich  
mpi4py                       4.0.3            py310hb6b6513_0  
mpich                        4.1.1            hbae89fd_0  
msgpack                      1.1.1            pypi_0           pypi  
multicoretsne                0.1              pypi_0           pypi  
multidict                    6.6.3            pypi_0           pypi  
multiprocessing-on-dill      3.5.0a4          pypi_0           pypi  
narwhals                     1.46.0           pypi_0           pypi  
natsort                      8.4.0            pypi_0           pypi  
nbclient                     0.10.2           py310h06a4308_0  
nbconvert-core               7.16.6           py310h06a4308_0  
nbformat                     5.10.4           py310h06a4308_0  
ncurses                      6.4              h6a678d5_0  
nest-asyncio                 1.6.0            py310h06a4308_0  
networkx                     3.4.2            py310h06a4308_0  
notebook-shim                0.2.4            py310h06a4308_0  
numba                        0.61.2           py310h6a678d5_0  
numexpr                      2.10.1           py310h3c60e43_0  
numpy                        2.0.2            pypi_0           pypi  
numpy-base                   2.2.5            py310h06ae042_0  
openssl                      3.0.16           h5eee18b_0  
overrides                    7.4.0            py310h06a4308_0  
packaging                    24.2             py310h06a4308_0  
pandas                       2.2.3            py310h6a678d5_0  
pandocfilters                1.5.0            pyhd3eb1b0_0  
parso                        0.8.4            py310h06a4308_0  
partd                        1.4.2            pypi_0           pypi  
patsy                        1.0.1            pypi_0           pypi  
pexpect                      4.9.0            py310h06a4308_0  
pickleshare                  0.7.5            pypi_0           pypi  
pillow                       11.3.0           pypi_0           pypi  
pip                          25.1             pyhc872135_2  
pipreqs                      0.5.0            pypi_0           pypi  
platformdirs                 4.3.7            py310h06a4308_0  
prometheus_client            0.21.1           py310h06a4308_0  
prompt-toolkit               3.0.43           py310h06a4308_0  
prompt_toolkit               3.0.43           hd3eb1b0_0  
propcache                    0.3.2            pypi_0           pypi  
psutil                       5.9.0            py310h5eee18b_1  
pthread-stubs                0.3              h0ce48e5_1  
ptyprocess                   0.7.0            pyhd3eb1b0_2  
pure_eval                    0.2.2            pyhd3eb1b0_0  
pyarrow                      20.0.0           pypi_0           pypi  
pycparser                    2.21             pyhd3eb1b0_0  
pygments                     2.19.1           py310h06a4308_0  
pynndescent                  0.5.13           pypi_0           pypi  
pyparsing                    3.2.3            pypi_0           pypi  
pyscenic                     0.12.1           pypi_0           pypi  
pysocks                      1.7.1            py310h06a4308_0  
python                       3.10.18          h1a3bd86_0  
python-dateutil              2.9.0post0       py310h06a4308_2  
python-fastjsonschema        2.20.0           py310h06a4308_0  
python-json-logger           3.2.1            py310h06a4308_0  
python-tzdata                2025.2           pyhd3eb1b0_0  
python_abi                   3.10             2_cp310          cctbx202211  
pytz                         2024.1           py310h06a4308_0  
pyyaml                       6.0.2            py310h5eee18b_0  
pyzmq                        26.2.0           py310h6a678d5_0  
readline                     8.2              h5eee18b_0  
referencing                  0.30.2           py310h06a4308_0  
requests                     2.32.4           py310h06a4308_0  
rfc3339-validator            0.1.4            py310h06a4308_0  
rfc3986-validator            0.1.1            py310h06a4308_0  
rpds-py                      0.22.3           py310h4aa5aa6_0  
scanpy                       1.11.3           pypi_0           pypi  
scikit-learn                 1.7.0            pypi_0           pypi  
scipy                        1.15.3           py310h525edd1_0  
seaborn                      0.13.2           pypi_0           pypi  
send2trash                   1.8.2            py310h06a4308_1  
session-info2                0.1.2            pypi_0           pypi  
setuptools                   72.1.0           py310h06a4308_0  
six                          1.17.0           py310h06a4308_0  
sniffio                      1.3.0            py310h06a4308_0  
sortedcontainers             2.4.0            pypi_0           pypi  
soupsieve                    2.5              py310h06a4308_0  
sqlite                       3.45.3           h5eee18b_0  
stack_data                   0.2.0            pyhd3eb1b0_0  
statsmodels                  0.14.5           pypi_0           pypi  
tbb                          2021.8.0         hdb19cb5_0  
tblib                        3.1.0            pypi_0           pypi  
terminado                    0.17.1           py310h06a4308_0  
threadpoolctl                3.6.0            pypi_0           pypi  
tinycss2                     1.4.0            py310h06a4308_0  
tk                           8.6.14           h993c535_1  
tomli                        2.0.1            py310h06a4308_0  
toolz                        1.0.0            pypi_0           pypi  
tornado                      6.5.1            py310h5eee18b_0  
tqdm                         4.67.1           py310h2f386ee_0  
traitlets                    5.14.3           py310h06a4308_0  
typing                       3.10.0.0         py310h06a4308_0  
typing-extensions            4.12.2           py310h06a4308_0  
typing_extensions            4.12.2           py310h06a4308_0  
tzdata                       2025b            h04d1e81_0  
umap-learn                   0.5.9.post2      pypi_0           pypi  
urllib3                      2.5.0            py310h06a4308_0  
wcwidth                      0.2.13           py310h06a4308_0  
webencodings                 0.5.1            py310h06a4308_1  
websocket-client             1.8.0            py310h06a4308_0  
wheel                        0.45.1           py310h06a4308_0  
xorg-libx11                  1.8.12           h9b100fa_1  
xorg-libxau                  1.0.12           h9b100fa_0  
xorg-libxdmcp                1.1.5            h9b100fa_0  
xorg-xorgproto               2024.1           h5eee18b_1  
xyzservices                  2025.4.0         pypi_0           pypi  
xz                           5.6.4            h5eee18b_1  
yaml                         0.2.5            h7b6447c_0  
yarg                         0.1.9            pypi_0           pypi  
yarl                         1.20.1           pypi_0           pypi  
zeromq                       4.3.5            h6a678d5_0  
zict                         3.0.0            pypi_0           pypi  
zipp                         3.21.0           py310h06a4308_0  
zlib                         1.2.13           h5eee18b_1  
zstd                         1.5.6            hc292b87_0  
