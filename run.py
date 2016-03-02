import os
import sys

from qtl_viewer.application import app
from qtl_viewer.application import CONF
from qtl_viewer.utils import data_utils
from qtl_viewer.utils import search_utils

if __name__ == "__main__":

    settings_file = ''

    if len(sys.argv) > 1:
        settings_file = sys.argv[1]
    else:
        print 'Please supply a configurations/settings file.'
        sys.exit();

    settings_file = os.path.abspath(settings_file)

    print "Configuring application from: {}".format(settings_file)
    app.config.from_pyfile(settings_file)

    host = app.config['HOST'] if 'HOST' in app.config else '0.0.0.0'
    port = app.config['PORT'] if 'PORT' in app.config else 80
    threaded = app.config['THREADED'] if 'THREADED' in app.config else True

    # let's attempt to make sure all files are where they are supposed to be
    if 'DATA_BASE_DIR' in app.config:
        CONF.DATA_BASE_DIR = os.path.abspath(app.config['DATA_BASE_DIR'])
    else:
        CONF.DATA_BASE_DIR = os.path.abspath('.')

    if not os.path.exists(CONF.DATA_BASE_DIR):
        print 'Either CONF.DATA_BASE_DIR is not set or does not exist, please check settings'
        print 'CONF.DATA_BASE_DIR={}'.format(CONF.DATA_BASE_DIR)
        sys.exit(1)

    CONF.DATA_CHROMOSOMES = os.path.join(CONF.DATA_BASE_DIR, app.config['DATA_CHROMOSOMES'])
    if not os.path.exists(CONF.DATA_CHROMOSOMES):
        print 'Either CONF.DATA_CHROMOSOMES is not set or does not exist, please check settings'
        sys.exit(1)

    CONF.DATA_SEARCH_DB = os.path.join(CONF.DATA_BASE_DIR, app.config['DATA_SEARCH_DB'])
    if not os.path.exists(CONF.DATA_SEARCH_DB):
        print 'Either CONF.DATA_SEARCH_DB is not set or does not exist, please check settings'
        sys.exit(1)

    CONF.DATA_HDF5 = os.path.join(CONF.DATA_BASE_DIR, app.config['DATA_HDF5'])
    if not os.path.exists(CONF.DATA_HDF5):
        print 'Either CONF.DATA_HDF5 is not set or does not exist, please check settings'
        sys.exit(1)

    search_utils.DATABASE = CONF.DATA_SEARCH_DB
    data_utils.HDF5_FILENAME = CONF.DATA_HDF5
    data_utils.init()

    CONF.URL_BASE = app.config['URL_BASE']
    CONF.URL_STATIC = app.config['URL_BASE_STATIC']

    if 'DEFAULT_DATASET' in app.config:
        CONF.DEFAULT_DATASET = app.config['DEFAULT_DATASET']
    else:
        datasets = data_utils.get_datasets()
        CONF.DEFAULT_DATASET = datasets.keys()[0]

    print 'Attempting to start application on: {}:{}'.format(host, port)
    print 'Use {}'.format(CONF.URL_BASE)

    # configure all WWW settings
    CONF.WWW_APP_HEADER = app.config['WWW_APP_HEADER']
    CONF.WWW_DATASET_TEXT = app.config['WWW_DATASET_TEXT']
    CONF.WWW_DATASET_DEFAULT = app.config['WWW_DATASET_DEFAULT']
    CONF.WWW_MATRIX_HEADER = app.config['WWW_MATRIX_HEADER']
    CONF.WWW_MATRIX_DATA_TAB_TEXT = app.config['WWW_MATRIX_DATA_TAB_TEXT']
    CONF.WWW_MATRIX_LIST_TAB_TEXT = app.config['WWW_MATRIX_LIST_TAB_TEXT']
    CONF.WWW_MATRIX_LOD = app.config['WWW_MATRIX_LOD']
    CONF.WWW_LOD_THRESHOLD_SLIDER_TEXT = app.config['WWW_LOD_THRESHOLD_SLIDER_TEXT']
    CONF.WWW_LOD_THRESHOLD_SLIDER_DEFAULT = app.config['WWW_LOD_THRESHOLD_SLIDER_DEFAULT']
    CONF.WWW_LOD_THRESHOLD_SLIDER_MIN = app.config['WWW_LOD_THRESHOLD_SLIDER_MIN']
    CONF.WWW_LOD_THRESHOLD_SLIDER_MAX = app.config['WWW_LOD_THRESHOLD_SLIDER_MAX']
    CONF.WWW_LOD_THRESHOLD_SLIDER_STEP = app.config['WWW_LOD_THRESHOLD_SLIDER_STEP']
    CONF.WWW_MATRIX_X_AXIS_TEXT = app.config['WWW_MATRIX_X_AXIS_TEXT']
    CONF.WWW_MATRIX_Y_AXIS_TEXT = app.config['WWW_MATRIX_Y_AXIS_TEXT']
    CONF.WWW_SEARCH_HEADER = app.config['WWW_SEARCH_HEADER']
    CONF.WWW_SEARCH_SEARCH_TAB_TEXT = app.config['WWW_SEARCH_SEARCH_TAB_TEXT']
    CONF.WWW_SEARCH_HISTORY_TAB_TEXT = app.config['WWW_SEARCH_HISTORY_TAB_TEXT']
    CONF.WWW_SEARCH_BUTTON_FIND_TEXT = app.config['WWW_SEARCH_BUTTON_FIND_TEXT']
    CONF.WWW_LOD_HEADER = app.config['WWW_LOD_HEADER']
    CONF.WWW_LOD_X_AXIS_TEXT = app.config['WWW_LOD_X_AXIS_TEXT']
    CONF.WWW_LOD_Y_AXIS_TEXT = app.config['WWW_LOD_Y_AXIS_TEXT']
    CONF.WWW_EFFECT_HEADER = app.config['WWW_EFFECT_HEADER']
    CONF.WWW_EFFECT_X_AXIS_TEXT = app.config['WWW_EFFECT_X_AXIS_TEXT']
    CONF.WWW_EFFECT_Y_AXIS_TEXT = app.config['WWW_EFFECT_Y_AXIS_TEXT']
    CONF.WWW_FACT_HEADER = app.config['WWW_FACT_HEADER']
    CONF.WWW_FACT_FACT_ORDER_HEADER = app.config['WWW_FACT_FACT_ORDER_HEADER']
    CONF.WWW_FACT_FACT_ORDER_TEXT = app.config['WWW_FACT_FACT_ORDER_TEXT']
    CONF.WWW_FACT_BUTTON_PLOT_TEXT = app.config['WWW_FACT_BUTTON_PLOT_TEXT']

    app.run(host=host, port=port, threaded=threaded)
