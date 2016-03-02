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

import sqlite3
import re

from flask import g

DATABASE = ''


def connect_to_database():
    return sqlite3.connect(DATABASE)

def get_db():
    db = getattr(g, '_database', None)
    if db is None:
        db = g._database = connect_to_database()
    db.row_factory = sqlite3.Row
    return db

REGEX_ENSEMBL_MOUSE_ID = re.compile("ENSMUS([EGTP])[0-9]{11}", re.IGNORECASE)
REGEX_ENSEMBL_HUMAN_ID = re.compile("ENS([EGTP])[0-9]{11}", re.IGNORECASE)
REGEX_MGI_ID = re.compile("MGI:[0-9]{1,}", re.IGNORECASE)
REGEX_LOCATION = re.compile("(CHR|)*\s*([0-9]{1,2}|X|Y|MT)\s*(-|:)?\s*(\d+)\s*(MB|M|K|)?\s*(-|:|)?\s*(\d+|)\s*(MB|M|K|)?", re.IGNORECASE)


SQL_TERM_EXACT = '''
SELECT MAX(s.text_score||'||'||s.score_description||'||'||e.text_value) AS match_description, e.ensembl_gene_id, g.*
  FROM ensembl_genes g,
       ensembl_genes_lookup e,
       text_type_score s
WHERE g.ensembl_gene_id = e.ensembl_gene_id
  AND e.text_type = s.text_type
  AND e.text_value = :term
GROUP BY e.ensembl_gene_id
ORDER BY match_description - length(match_description) DESC, g.symbol ASC
'''

SQL_TERM_EXACT_MM = '''
SELECT MAX(s.text_score||'||'||s.score_description||'||'||e.text_value) AS match_description, e.ensembl_gene_id, g.*
  FROM ensembl_genes g,
       ensembl_genes_lookup e,
       text_type_score s
WHERE g.ensembl_gene_id = e.ensembl_gene_id
  AND e.text_type = s.text_type
  AND e.text_value = :term
  AND e.species_id = 'Mm'
GROUP BY e.ensembl_gene_id
ORDER BY match_description - length(match_description) DESC, g.symbol ASC
'''

SQL_TERM_EXACT_HS = '''
SELECT MAX(s.text_score||'||'||s.score_description||'||'||e.text_value) AS match_description, e.ensembl_gene_id, g.*
  FROM ensembl_genes g,
       ensembl_genes_lookup e,
       text_type_score s
WHERE g.ensembl_gene_id = e.ensembl_gene_id
  AND e.text_type = s.text_type
  AND e.text_value = :term
  AND e.species_id = 'Hs'
GROUP BY e.ensembl_gene_id
ORDER BY match_description - length(match_description) DESC, g.symbol ASC
'''

SQL_TERM_LIKE = '''
SELECT MAX(s.text_score||'||'||s.score_description||'||'||e.text_value) AS match_description, e.ensembl_gene_id, g.*
  FROM ensembl_genes g,
       ensembl_genes_lookup e,
       text_type_score s
WHERE g.ensembl_gene_id = e.ensembl_gene_id
  AND e.text_type = s.text_type
  AND e.text_value like :term
GROUP BY e.ensembl_gene_id
ORDER BY match_description - length(match_description) DESC, g.symbol ASC
'''

SQL_TERM_LIKE_MM = '''
SELECT MAX(s.text_score||'||'||s.score_description||'||'||e.text_value) AS match_description, e.ensembl_gene_id, g.*
  FROM ensembl_genes g,
       ensembl_genes_lookup e,
       text_type_score s
WHERE g.ensembl_gene_id = e.ensembl_gene_id
  AND e.text_type = s.text_type
  AND e.text_value like :term
  AND e.species_id = 'Mm'
GROUP BY e.ensembl_gene_id
ORDER BY match_description - length(match_description) DESC, g.symbol ASC
'''

SQL_TERM_LIKE_HS = '''
SELECT MAX(s.text_score||'||'||s.score_description||'||'||e.text_value) AS match_description, e.ensembl_gene_id, g.*
  FROM ensembl_genes g,
       ensembl_genes_lookup e,
       text_type_score s
WHERE g.ensembl_gene_id = e.ensembl_gene_id
  AND e.text_type = s.text_type
  AND e.text_value like :term
  AND e.species_id = 'Hs'
GROUP BY e.ensembl_gene_id
ORDER BY match_description - length(match_description) DESC, g.symbol ASC
'''

SQL_ID = '''
SELECT MAX(s.text_score||'||'||s.score_description||'||'||e.text_value) AS match_description, e.ensembl_gene_id, g.*
  FROM ensembl_genes g,
       ensembl_genes_lookup e,
       text_type_score s
WHERE g.ensembl_gene_id = e.ensembl_gene_id
  AND e.text_type = s.text_type
  AND e.text_value = :term
  AND e.text_type in ('MI', 'EG', 'ET', 'EE', 'HI')
GROUP BY e.ensembl_gene_id
ORDER BY match_description - length(match_description) DESC, g.symbol ASC
'''

SQL_ID_MM = '''
SELECT MAX(s.text_score||'||'||s.score_description||'||'||e.text_value) AS match_description, e.ensembl_gene_id, g.*
  FROM ensembl_genes g,
       ensembl_genes_lookup e,
       text_type_score s
WHERE g.ensembl_gene_id = e.ensembl_gene_id
  AND e.text_type = s.text_type
  AND e.text_value = :term
  AND e.species_id = 'Mm'
  AND e.text_type in ('MI', 'EG', 'ET', 'EE', 'HI')
GROUP BY e.ensembl_gene_id
ORDER BY match_description - length(match_description) DESC, g.symbol ASC
'''

SQL_ID_HS = '''
SELECT MAX(s.text_score||'||'||s.score_description||'||'||e.text_value) AS match_description, e.ensembl_gene_id, g.*
  FROM ensembl_genes g,
       ensembl_genes_lookup e,
       text_type_score s
WHERE g.ensembl_gene_id = e.ensembl_gene_id
  AND e.text_type = s.text_type
  AND e.text_value = :term
  AND e.species_id = 'Hs'
  AND e.text_type in ('MI', 'EG', 'ET', 'EE', 'HI')
GROUP BY e.ensembl_gene_id
ORDER BY match_description - length(match_description) DESC, g.symbol ASC
'''

SQL_LOCATION = '''
SELECT *
  FROM ensembl_genes e
 WHERE e.chromosome = :chromosome
   AND e.start_position >= :start_position
   AND e.end_position <= :end_position
 ORDER BY cast(replace(replace(replace(e.chromosome, 'X', '50'), 'Y', '51'), 'MT', 51) AS int), e.start_position, e.end_position
'''

SQL_LOCATION_MM = '''
SELECT *
  FROM ensembl_genes e
 WHERE e.chromosome = :chromosome
   AND e.start_position >= :start_position
   AND e.end_position <= :end_position
   AND e.species_id = 'Mm'
 ORDER BY cast(replace(replace(replace(e.chromosome, 'X', '50'), 'Y', '51'), 'MT', 51) AS int), e.start_position, e.end_position
'''

SQL_LOCATION_HS = '''
SELECT *
  FROM ensembl_genes e
 WHERE e.chromosome = :chromosome
   AND e.start_position >= :start_position
   AND e.end_position <= :end_position
   AND e.species_id = 'Hs'
 ORDER BY cast(replace(replace(replace(e.chromosome, 'X', '50'), 'Y', '51'), 'MT', 51) AS int), e.start_position, e.end_position
'''


QUERIES = {}
QUERIES['SQL_TERM_EXACT']=SQL_TERM_EXACT
QUERIES['SQL_TERM_EXACT_MM']=SQL_TERM_EXACT_MM
QUERIES['SQL_TERM_EXACT_HS']=SQL_TERM_EXACT_HS
QUERIES['SQL_TERM_LIKE']=SQL_TERM_LIKE
QUERIES['SQL_TERM_LIKE_MM']=SQL_TERM_LIKE_MM
QUERIES['SQL_TERM_LIKE_HS']=SQL_TERM_LIKE_HS
QUERIES['SQL_ID']=SQL_ID
QUERIES['SQL_ID_MM']=SQL_ID_MM
QUERIES['SQL_ID_HS']=SQL_ID_HS
QUERIES['SQL_LOCATION']=SQL_LOCATION
QUERIES['SQL_LOCATION_MM']=SQL_LOCATION_MM
QUERIES['SQL_LOCATION_HS']=SQL_LOCATION_HS


class Status:
    pass


class Location:
    """ Encapsulates a genomic location

    """
    def __init__(self):
        self.chromosome = ''
        self.start_position = None
        self.end_position = None
    def __str__(self):
        return str(self.chromosome) + ':' + str(self.start_position) + '-' + str(self.end_position)
    def __repr__(self):
        location = {}
        location['chromosome'] = self.chromosome
        location['start_position'] = self.start_position
        location['end_position'] = self.end_position
        return location


class Query:
    """ Simple class to encapsulate query objects

    """
    def __init__(self, term=None, species_id=None, exact=False, verbose=False):
        self.term = term
        self.species_id = species_id
        self.exact = exact
        self.verbose = verbose
        self._location = None
        self.query = None

    def __str__(self):
        return 'QUERY: \n' + str(self.query) + '\nTERM: ' + str(self.term) + '\nSPECIES_ID: ' + str(self.species_id) + \
               '\nEXACT: ' + str(self.exact) + '\nVERBOSE: ' + str(self.verbose) + '\nLOCATION: ' + str(self._location) + \
               '\nPARAMS: ' + str(self.get_params())
    def get_params(self):
        if self._location:
            return {'chromosome':self._location.chromosome, 'start_position':self._location.start_position, 'end_position':self._location.end_position}
        return {'term': self.term}


class Match:
    """ Simple class to encapsulate a match

    """
    def __init__(self, ensembl_gene_id=None, external_id=None, symbol=None, name=None, description=None,
                 synonyms=None, species_id= None, chromosome=None, position_start=None, position_end=None, strand=None,
                 match_reason=None, match_value=None):
        self.ensembl_gene_id = ensembl_gene_id
        self.external_id = external_id
        self.species_id = species_id
        self.symbol = symbol
        self.name = name
        self.description = description
        self.synonyms = synonyms
        self.chromosome = chromosome
        self.position_start = position_start
        self.position_end = position_end
        self.strand = strand
        self.match_reason = match_reason
        self.match_value = match_value

    def __str__(self):
        return str(self.ensembl_gene_id)


class Result:
    """ Simple class to encapsulate a Query and matches

    """
    def __init__(self, query=None, matches=[]):
        self.query = query
        self.matches = matches


def str2bool(v):
    """ Checks to see if v is a string representing True

    Parameters:
        v: a string representing a boolean value
    Returns:
        True or False
    """
    if v:
        return v.lower() in ("yes", "true", "t", "1")
    return False


def nvl(value, default):
    """ Evaluates if value es empty or None, if so returns default

    Parameters:
        value: the evalue to evaluate
        default: the default value
    Returns:
        value or default
    """
    if value:
        return value
    return default


def nvli(value, default):
    ret = default
    if value:
        try:
            ret = int(value)
        except ValueError:
            pass
    return ret


def get_status(error = False, message=None):
    """ Create a status with no error and no message

    Parameters:
        error: True if there is an error, False otherwise
        message: The error message
    Returns:
        status: Status object
    """
    _status = Status()
    _status.error = error
    _status.message = nvl(message, '')
    return _status


def get_multiplier(factor):
    if factor.lower() in ['g', 'gb', 'gbp']:
        return 1000000000
    if factor.lower() in ['m', 'mb', 'mbp']:
        return 1000000
    elif factor.lower() in ['k', 'kb', 'kbp']:
        return 1000

    return 1


def str_to_location(location):
    status = get_status()
    if not location:
        return None, get_status(True, 'No location')

    valid_location = str(location).strip()

    if len(valid_location) <= 0:
        return None, get_status(True, 'Empty location')

    match = REGEX_LOCATION.match(valid_location)

    loc = None

    if match:
        loc = Location()
        loc.chromosome = match.group(2)
        loc.start_position = match.group(4)
        loc.end_position = match.group(7)
        multiplier_one = match.group(5)
        multiplier_two = match.group(8)

        if type(loc.start_position) is str:
            loc.start_position = long(loc.start_position)
        if type(loc.end_position) is str:
            loc.end_position = long(loc.end_position)

        if multiplier_one:
            loc.start_position = loc.start_position * get_multiplier(multiplier_one)

        if multiplier_two:
            loc.end_position = loc.end_position * get_multiplier(multiplier_two)
    else:
        status = get_status(True, 'Invalid location string')

    return loc, status


def convert(input):
    if isinstance(input, dict):
        return {convert(key): convert(value) for key, value in input.iteritems()}
    elif isinstance(input, list):
        return [convert(element) for element in input]
    elif isinstance(input, unicode):
        return input.encode('utf-8')
    else:
        return input


def validate_ensembl_id(ensembl_id):
    """ Validate an id to make sure it conforms to the convention.

    Parameters:
        ensembl_id: the ensembl identifer (string)
    Returns:
        valid_id: the valid ensembl id or None if not valid
        status: True if error, false otherwise and message

    """
    status = get_status()
    if not ensembl_id:
        return None, get_status(True, 'No Ensembl id')

    valid_id = ensembl_id.strip()

    if len(valid_id) <= 0:
        return None, get_status(True, 'Empty Ensembl id')

    if REGEX_ENSEMBL_HUMAN_ID.match(ensembl_id):
        return valid_id, status
    elif REGEX_ENSEMBL_MOUSE_ID.match(id):
        return valid_id, status

    return None, get_status(True, 'Invalid Ensembl ID')


def _get_query(term, species_id=None, exact=True, verbose=False):
    """ Get query based upon parameters

    Parameters:
        term: the query object
        species_id: either 'Hs', 'Mm', or None
        exact: True for exact matches
        verbose: True to give details
    Returns:
        query: the query
        status: True if error, false otherwise and message

    """
    status = get_status()
    if not term:
        return None, get_status(True, 'No term')

    valid_term = str(term).strip()

    if len(valid_term) <= 0:
        return None, get_status(True, 'Empty term')

    query = Query(term, species_id, exact, verbose)

    if species_id:
        species_id = '_' + species_id.upper()
    else:
        species_id = ''

    if REGEX_LOCATION.match(valid_term):
        location, status = str_to_location(valid_term)
        if not status.error:
            query.query = QUERIES['SQL_LOCATION' + species_id]
            query._location = location
        else:
            return None, status
    elif REGEX_ENSEMBL_MOUSE_ID.match(valid_term):
        query.query = QUERIES['SQL_ID' + species_id]
    elif REGEX_ENSEMBL_HUMAN_ID.match(valid_term):
        query.query = QUERIES['SQL_ID' + species_id]
    elif REGEX_MGI_ID.match(valid_term):
        query.query = QUERIES['SQL_ID' + species_id]
    else:
        if exact:
            query.query = QUERIES['SQL_TERM_EXACT' + species_id]
        else:
            query.query = QUERIES['SQL_TERM_LIKE' + species_id]

            if valid_term != '%':
                valid_term = '%' + valid_term
            if valid_term[-1] != '%':
                valid_term = valid_term + '%'

            query.term = valid_term

    return query, status


def _query(query):
    status = get_status()

    if not query:
        return None, get_status(True, 'No query')

    matches = []

    try:
        cursor = get_db().cursor()

        for row in cursor.execute(query.query, query.get_params()):
            match = Match()
            match.ensembl_gene_id = row['ensembl_gene_id']
            match.external_id = row['external_id']
            match.species_id = row['species_id']
            match.symbol = row['symbol']
            match.name = row['name']
            match.description = row['description']

            row_synonyms = row['synonyms']
            synonyms = []
            if row_synonyms:
                synonyms = row_synonyms.split('||')

            match.synonyms = synonyms
            match.chromosome = row['chromosome']
            match.position_start = row['start_position']
            match.position_end = row['end_position']
            match.strand = row['strand']

            if query._location:
                match.match_reason = 'Location'
                match.match_value = str(match.chromosome) + ':' + str(match.position_start) + '-' + str(match.position_end)
            else:
                row_match_description = row['match_description']
                if row_match_description:
                    desc = row_match_description.split('||')
                match.match_reason = desc[1]
                match.match_value = desc[2]


            matches.append(match)

        cursor.close()
    except sqlite3.Error, e:
        return None, get_status(True, 'Database Error: ' + str(e))

    return Result(query, matches), status


def search(term, species_id=None, exact=True, verbose=False):
    query, status = _get_query(term, species_id, exact, verbose)

    if verbose:
        print str(query)

    if status.error:
        return None, status

    result, status = _query(query)

    if status.error:
        return None, status

    return result, status




if __name__ == '__main__':
    #results, status = search("CHR1:1MB-1000MB", 'Mm', False, True)
    loc,status=str_to_location("CHR1:10MB-11MB")
    print loc






