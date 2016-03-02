import os

from flask import Blueprint, request, abort, jsonify, make_response, render_template, send_file, g
import numpy as np

from qtl_viewer.utils import search_utils
from qtl_viewer.utils import data_utils

mod = Blueprint('api', __name__, url_prefix='/api')


RESPONSE_OK = 200
RESPONSE_CREATED = 201
RESPONSE_BAD = 400
RESPONSE_NOT_FOUND = 404


def str2bool(v):
    return v.lower() in ("yes", "true", "t", "1")


@mod.errorhandler(404)
def not_found(error):
    return '', 404


@mod.route('/chromosomes', methods=['GET'])
def chromosomes():
    # POST search_term = request.form.get("search_term")
    # GET  search_term = request.args.get("search_term")
    # BOTH search_term = request.values.get("search_term")
    return send_file(g.CONF.DATA_CHROMOSOMES)


@mod.route('/matrix/<dataset>', methods=['GET'])
def matrix(dataset=None):
    if dataset is None:
        datasets = data_utils.get_datasets()
        dataset = datasets.keys()[0]
    filename = os.path.join(g.CONF.DATA_BASE_DIR, '{}.tsv'.format(dataset))
    print filename
    return send_file(filename)


@mod.route('/effect/<dataset>/<feature_id>', methods=['GET'])
def effect(dataset, feature_id):
    # POST search_term = request.form.get("search_term")
    # GET  search_term = request.args.get("search_term")
    # BOTH search_term = request.values.get("search_term")
    return data_utils.get_effect_data(dataset, feature_id)


@mod.route('/lod/<dataset>/<feature_id>', methods=['GET'])
def lod(dataset, feature_id):
    # POST search_term = request.form.get("search_term")
    # GET  search_term = request.args.get("search_term")
    # BOTH search_term = request.values.get("search_term")
    lod_data = data_utils.get_lod_data(dataset, feature_id)
    ret = ['marker_id\tchrom\tlocation\tlod']
    for d in lod_data:
        ret.append('\t'.join(map(str, d)))

    # TODO: check that value is not None or no Error thrown
    # TODO: pass back data, not string
    return '\n'.join(ret)


@mod.route('/factexp/<dataset>/<feature_id>/<marker_id>/<factors>', methods=['GET'])
def factexp(dataset, feature_id, marker_id, factors):
    # POST search_term = request.form.get("search_term")
    # GET  search_term = request.args.get("search_term")
    # BOTH search_term = request.values.get("search_term")

    # OR THIS (SAMPLES, FACT VIEWER)

    # dataset	/samples			table of samples, 1st col is unique identifier

    # group		/phenotypes
    # dataset		/phenotypes/factors		table of factors (sex, diet, tissue), 1st col is unique identifier
    # dataset		/phenotypes/phenotypes	table of phenotypes samples (rows) x factors (columns)

    # group		/genotypes
    # dataset		/genotypes/genotypes	table of genotypes markers (rows) x samples (columns)

    # group		/expression
    # dataset		/expression/expression	table of expression values features (rows) x samples (columns)

    factors_phenotypes = data_utils.get_factors_phenotypes(dataset, marker_id)
    samples = data_utils.get_samples(dataset)
    expression = data_utils.get_expression(dataset, feature_id)

    samples_arr = []

    for s in samples:
        if s['name'] is not None:
            samples_arr.append(s['name'])
        else:
            samples_arr.append(s['sample_id'])

    factors_arr = factors.split(',')

    #
    # check to make sure we actually have an expression value to plot
    #
    nan_indices = np.where(np.isnan(expression)==True)
    expression = list(np.delete(expression, nan_indices))
    samples_arr = list(np.delete(samples_arr, nan_indices))

    xvar = []

    for factor_id in factors_arr:
        factor_pheno = factors_phenotypes[factor_id]
        factor_pheno_values = list(np.delete(factor_pheno['values'], nan_indices))
        f = {"kind": "factor",
             "name": factor_pheno['name'],
             "levels": list(factor_pheno['levels']),
             "values": factor_pheno_values}
        xvar.append(f)

    yvar = [{
        "kind": "number",
        "name": feature_id,
        "values" : expression
    }]

    d = {
      "samples": samples_arr,
      "xvar": xvar,
      "yvar": yvar
    }

    return jsonify(d)


@mod.route("/pvalue2lod", methods=['GET','POST'])
def pvalue2lod():
    pvalue = float(request.values.get('pvalue', None))
    df = float(request.values.get('df', None))
    return str(data_utils.pvalue2lod(pvalue, df))


@mod.route("/search", methods=['GET','POST'])
def ws_search():
    term = search_utils.nvl(request.values.get('term', None), '')
    species = search_utils.nvl(request.values.get('species', None), '')
    exact = search_utils.str2bool(request.values.get('exact', 'True'))
    verbose = search_utils.str2bool(request.values.get('verbose', 'False'))
    callbackFunction = request.values.get('callback', None)

    valid_term = term.strip()
    result = None

    if len(valid_term) == 0:
        status = search_utils.get_status(True, 'No term specified')
    else:
        result, status = search_utils.search(term=valid_term, exact=exact, species_id=species, verbose=verbose)

    response = make_response(render_template('api/searchresults.json', callbackFunction=callbackFunction, result=result, status=status))
    response.headers['Content-Type'] = 'application/json'

    return response
