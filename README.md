# qtl-viewer

The qtl-viewer is a web application for viewing QTL data.


## Requirements

python >= 2.7

### Python requirements:

Cython==0.23.4

Flask==0.10.1

Jinja2==2.8

MarkupSafe==0.23

Werkzeug==0.11.4

h5py==2.5.0

itsdangerous==0.24

numpy==1.10.4

scipy==0.17.0

six==1.10.0

wsgiref==0.1.2


## Preparation:

__For format of the HDF5 file, please see the docs subdirectory.__

The qtl-viewer needs a *settings.cfg* file for configuration.  Please see an example *settings.cfg* file in the docs sub directory.

Once the settings file has been configured, the data files used to generate the matrix grid will need to be generated.

This can be done by issuing the following command:

`python qtl_viewer\utils\data_utils.py [HDF5 file]`

Make sure the files are in the specified `DATA_BASE_DIR` from the *settings.cfg* file.


## Running

Once the data files have been generated, the application can be started with the following command:

`python run.py [settings.cfg]`


## Example

<img src="ftp://ftp.jax.org/mvincent/qtl-viewer/qtl-viewer-screenshot.png">

### Example Data

Example data and a settings file can be located at:

ftp://ftp.jax.org/mvincent/qtl-viewer/


