from flask import Blueprint, request, url_for, redirect, render_template, g
from qtl_viewer.utils import search_utils
from qtl_viewer.utils import data_utils

mod = Blueprint('general', __name__)


@mod.errorhandler(404)
def not_found(error):
    return render_template('general/404.html'), 404


@mod.before_request
def load_current_user():
    """
    only need the user for this portion
    """
    #g.CONF = CONF
    pass


@mod.route("/")
def index():
    # POST search_term = request.form.get("search_term")
    # GET  search_term = request.args.get("search_term")
    # BOTH search_term = request.values.get("search_term")

    print g.CONF

    search_term = request.values.get('search_term', '')
    lod_threshold = request.values.get('lod_threshold', g.CONF.WWW_LOD_THRESHOLD_SLIDER_DEFAULT)
    dataset = request.values.get('dataset', g.CONF.DEFAULT_DATASET)
    all_datasets = data_utils.get_datasets()

    render_effect_plot = not data_utils.has_samples(dataset)
    if render_effect_plot:
        strains = data_utils.get_strains(dataset)
        factors = None
    else:
        strains = None
        factors = data_utils.get_factors_web(dataset)

    pvalueToLodTable = {}
    if not g.CONF.WWW_MATRIX_LOD:
        for uk,v in all_datasets.iteritems():

            k = str(uk)
            pvalueToLodTable[k] = data_utils.pvalue2lodtable(k,
                                                             g.CONF.WWW_LOD_THRESHOLD_SLIDER_MIN,
                                                             g.CONF.WWW_LOD_THRESHOLD_SLIDER_MAX)

    return render_template('general/index.html', search_term=search_term, lod_threshold=lod_threshold,
                           all_datasets=all_datasets, dataset=dataset, strains=strains,
                           render_effect_plot=render_effect_plot, factors=factors,
                           pvalueToLodTable=pvalueToLodTable)
