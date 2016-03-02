# -*- coding: utf-8 -*-

"""
Copyright (c) 2015 Matthew Vincent, The Churchill Lab

This software was developed by Gary Churchill's Lab.
(see http://research.jax.org/faculty/churchill)

This is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This software is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this software. If not, see <http://www.gnu.org/licenses/>.
"""

from flask import Flask
from flask import Response
from flask import abort
from flask import g
from flask import jsonify
from flask import make_response
from flask import render_template
from flask import request
from flask import send_file
from flask.ext.compress import Compress
from flask.ext.cors import CORS

from flask_debugtoolbar import DebugToolbarExtension

from utils import search_utils

import h5py

import rpy2.robjects as ro

import os
import sys

class Config:
    pass

app = Flask(__name__)

app.debug = True

#compress = Compress()
#compress.init_app(app)


#CORS(app)

CONF = Config()

HDF_DATA = None
R_DATA = None


from flask import after_this_request, request
from cStringIO import StringIO as IO
import gzip
import functools

def gzipped(f):
    @functools.wraps(f)
    def view_func(*args, **kwargs):
        @after_this_request
        def zipper(response):
            accept_encoding = request.headers.get('Accept-Encoding', '')

            if 'gzip' not in accept_encoding.lower():
                return response

            response.direct_passthrough = False

            if (response.status_code < 200 or
                response.status_code >= 300 or
                'Content-Encoding' in response.headers):
                return response
            gzip_buffer = IO()
            gzip_file = gzip.GzipFile(mode='wb',
                                      fileobj=gzip_buffer)
            gzip_file.write(response.data)
            gzip_file.close()

            response.data = gzip_buffer.getvalue()
            response.headers['Content-Encoding'] = 'gzip'
            response.headers['Vary'] = 'Accept-Encoding'
            response.headers['Content-Length'] = len(response.data)

            return response

        return f(*args, **kwargs)

    return view_func





@app.teardown_appcontext
def close_connection(exception):
    db = getattr(g, '_database', None)
    if db is not None:
        db.close()


@app.before_request
def before_request():
    request.URL_BASE = CONF.URL_BASE
    if 'HTTP_X_FORWARDED_PATH' in request.environ:
        protocol = 'http:'
        if 'wsgi.url_scheme' in request.environ:
            protocol = request.environ['wsgi.url_scheme'] + ':'

        request.URL_BASE = protocol + '//' + request.environ['HTTP_X_FORWARDED_HOST'] + request.environ['HTTP_X_FORWARDED_PATH']
        if request.URL_BASE[-1] == '/':
            request.URL_BASE = request.URL_BASE[:-1]


@app.route("/", methods=['GET', 'POST'])
def index():
    # POST search_term = request.form.get("search_term")
    # GET  search_term = request.args.get("search_term")
    # BOTH search_term = request.values.get("search_term")

    search_term = request.values.get('search_term', '')
    lod_threshold = request.values.get('lod_threshold', CONF.LOD_THRESHOLD)
    layout = request.values.get('layout', '1')

    print layout

    template_html = 'index.html'
    if layout == '2':
        template_html = 'index2.html'

    print template_html

    return render_template(template_html, CONF=CONF, search_term=search_term, lod_threshold=lod_threshold)


@app.route("/matrix", methods=['GET'])
def matrix():
    return send_file(CONF.DATA_MATRIX)


def get_lod_data(identifier):
    """
    Return tab delimited file

    rs  chromosome  position    lod

    :param identifier:
    :return:
    """
    data = []
    data.append('rs\tchromosome\tposition\tlod')

    idx = HDF_DATA['genes'][identifier]

    markers = HDF_DATA['data']['lod']['markers'][()]

    for i, l in enumerate(HDF_DATA['data']['lod']['lod'][idx]):
        line = list(markers[i])
        line.append(l)
        data.append('\t'.join(map(str, line)))

    return '\n'.join(data)


@app.route("/lod", methods=['GET'])
@gzipped
def lod():
    identifier = search_utils.nvl(request.args.get('identifier', None), '')
    if not identifier:
        return 'error'
    return get_lod_data(identifier)


def get_effect_data(identifier):
    idx = HDF_DATA['genes'][identifier]

    founders = list(HDF_DATA['data']['coef']['founders'][()])

    top = ['rs']
    top.extend(zip(*list(HDF_DATA['data']['lod']['markers'][()]))[0])

    data = []
    data.append(top)

    for i, c in enumerate(HDF_DATA['data']['coef']['coef'][idx]):
        line = [founders[i]]
        line.extend(c)
        data.append(line)
        #print line
        #print '\n\n\n'

    # data is now list of list of (founder, coef, coef, coef, ..., coef * number of markers)

    lines = []

    # transpose
    for d in zip(*data):
        lines.append("\t".join(map(str, d)))

    # lines is now list of (marker, founder val, founder val, ..., founder val * number of founders) * number markers

    return '\n'.join(lines)


@app.route("/effect", methods=['GET'])
@gzipped
def effect():
    identifier = search_utils.nvl(request.args.get('identifier', None), '')
    if not identifier:
        return 'error'
    return get_effect_data(identifier)


@app.route("/factexp", methods=["GET", "POST"])
def factexp():
    gene = request.values.get('gene', '')
    factors = request.values.get('factors', 'Sex,Age')

    x_variables = []

    for f in factors.split(','):
        r_factor = ro.r('pheno[, "{}"]'.format(f))
        j_factor = {"kind": "factor", "name": f}

        if isinstance(r_factor, ro.vectors.FactorVector):
            j_factor['levels'] = list(r_factor.levels)
            j_factor['values'] = list(r_factor.iter_labels())
        else:
            j_factor['levels'] = list(set(r_factor))
            j_factor['values'] = list(r_factor)

        x_variables.append(j_factor)

    y_variables = []
    j_factor = {"kind": "number", "name": gene, "values": list(ro.r('expr.mrna[,"{}"]'.format(gene)))}
    y_variables.append(j_factor)

    return jsonify(xvar=x_variables, yvar=y_variables, samples=len(j_factor["values"]))


@app.route("/ws/search", methods=['GET','POST'])
def ws_search():
    if request.method == 'POST':
        term = search_utils.nvl(request.form('term'), '')
        species = search_utils.nvl(request.form('species'), '')
        exact = search_utils.str2bool(request.form('exact'))
        verbose = search_utils.str2bool(request.form('verbose'))
        callbackFunction = request.form('callback')

    elif request.method == 'GET':
        term = search_utils.nvl(request.args.get('term', None), '')
        species = search_utils.nvl(request.args.get('species', None), '')
        exact = search_utils.str2bool(request.args.get('exact', 'True'))
        verbose = search_utils.str2bool(request.args.get('verbose', 'False'))
        callbackFunction = request.args.get('callback', None)
    else:
        abort(500)

    valid_term = term.strip()
    result = None

    if len(valid_term) == 0:
        status = search_utils.get_status(True, 'No term specified')
    else:
        result, status = search_utils.search(term=valid_term, exact=exact, species_id=species, verbose=verbose)

    response = make_response(render_template('json/searchresults.json', callbackFunction=callbackFunction, result=result, status=status))
    response.headers['Content-Type'] = 'application/json'

    return response


if __name__ == "__main__":
    """
    if len(sys.argv) == 1:
        app.config.from_pyfile('settings.cfg')
    else:
        app.config.from_pyfile(os.path.abspath(sys.argv[1]))

    #toolbar = DebugToolbarExtension(app)

    search_utils.DATABASE = app.config['SEARCH_DB']

    CONF.URL_BASE = app.config['URL_BASE']
    CONF.URL_STATIC = app.config['URL_BASE_STATIC']
    CONF.URL_CHROMOSOME = app.config['URL_CHROMOSOME']
    CONF.HEADER_APP = app.config['HEADER_APP']
    CONF.HEADER_MATRIX = app.config['HEADER_MATRIX']
    CONF.DATA_MATRIX = app.config['DATA_MATRIX']
    CONF.DATA_RDATA = app.config['DATA_RDATA']
    CONF.LOD_THRESHOLD = app.config['LOD_THRESHOLD']
    CONF.LOD_THRESHOLD_MIN = app.config['LOD_THRESHOLD_MIN']
    CONF.LOD_THRESHOLD_MAX = app.config['LOD_THRESHOLD_MAX']

    HDF_FILENAME = app.config['DATA_HDF']
    HDF_DATA = {
        'data': h5py.File(HDF_FILENAME),
        'genes': {}
    }

    # build a dictionary to look up the index of the gene by name
    for idx, gene in enumerate(HDF_DATA['data']['lod']['genes']):
        HDF_DATA['genes'][gene] = idx

    R_DATA = ro.r['load'](CONF.DATA_RDATA)

    print str(app.config)

    app.run(app.config['HOST'], app.config['PORT'], app.config['DEBUG'], threaded=True)
    """






