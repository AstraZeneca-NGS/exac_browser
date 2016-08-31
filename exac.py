#!/usr/bin/env python2

import itertools
import json
import os
from os.path import basename

import pymongo
import pysam
import gzip
from parsing import *
import logging
import lookups
import random
import sys
from utils import *

from flask import Flask, request, session, g, redirect, url_for, abort, render_template, flash, jsonify, send_from_directory
from flask.ext.compress import Compress
from flask.ext.runner import Runner
#from flask.ext.script import Manager, Server
from flask_errormail import mail_on_500

from flask import Response
from collections import defaultdict, OrderedDict
from werkzeug.contrib.cache import SimpleCache

from multiprocessing import Process
import glob
import sqlite3
import traceback
import time

logging.getLogger().addHandler(logging.StreamHandler())
logging.getLogger().setLevel(logging.INFO)

ADMINISTRATORS = (
    'exac.browser.errors@gmail.com',
)

app = Flask(__name__)
mail_on_500(app, ADMINISTRATORS)
Compress(app)
app.config['COMPRESS_DEBUG'] = True
cache = SimpleCache(default_timeout=60*60*24)

EXAC_FILES_DIRECTORY = '../exac_data/'
REGION_LIMIT = 1E5
EXON_PADDING = 50

# Load default config and override config from an environment variable
app.config.update(dict(
    DB_HOST='localhost',
    DB_PORT=27017,
    DB_NAME='exac', 
    DEBUG=True,
    SECRET_KEY='development key',
    LOAD_DB_PARALLEL_PROCESSES = 1,  # contigs assigned to threads, so good to make this a factor of 24 (eg. 2,3,4,6,8)
    FEATURES_FILE=os.path.join(os.path.dirname(__file__), EXAC_FILES_DIRECTORY, '%s', 'all_features.bed.gz'),
    CANONICAL_TRANSCRIPT_FILE=os.path.join(os.path.dirname(__file__), EXAC_FILES_DIRECTORY, '%s', 'canonical_transcripts.txt.gz'),
    OMIM_FILE=os.path.join(os.path.dirname(__file__), EXAC_FILES_DIRECTORY, '%s', 'omim_info.txt.gz'),
    SITES_VCFS=os.path.join(os.path.dirname(__file__), EXAC_FILES_DIRECTORY, '%s', 'vardict', '%s.vcf.gz'),
    BASE_COVERAGE_FILES=os.path.join(os.path.dirname(__file__), EXAC_FILES_DIRECTORY, '%s', 'coverage', '%s', '*.txt.gz'),
    DBNSFP_FILE=os.path.join(os.path.dirname(__file__), EXAC_FILES_DIRECTORY, 'dbNSFP2.6_gene.gz'),
    CONSTRAINT_FILE=os.path.join(os.path.dirname(__file__), EXAC_FILES_DIRECTORY, 'forweb_cleaned_exac_r03_march16_z_data_pLI.txt.gz'),
    MNP_FILE=os.path.join(os.path.dirname(__file__), EXAC_FILES_DIRECTORY, 'MNPs_NotFiltered_ForBrowserRelease.txt.gz'),

    # How to get a dbsnp142.txt.bgz file:
    #   wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b142_GRCh37p13/database/organism_data/b142_SNPChrPosOnRef_105.bcp.gz
    #   zcat b142_SNPChrPosOnRef_105.bcp.gz | awk '$3 != ""' | perl -pi -e 's/ +/\t/g' | sort -k2,2 -k3,3n | bgzip -c > dbsnp142.txt.bgz
    #   tabix -s 2 -b 3 -e 3 dbsnp142.txt.bgz
    DBSNP_FILE=os.path.join(os.path.dirname(__file__), EXAC_FILES_DIRECTORY, 'dbsnp142.txt.bgz'),

    READ_VIZ_DIR=os.path.join(os.path.dirname(__file__), EXAC_FILES_DIRECTORY, '%s', 'combined_bams', '%s'),
    READ_VIZ_DIR_HTML=os.path.join('/%s', 'combined_bams', '%s'),
))

GENE_CACHE_DIR = os.path.join(os.path.dirname(__file__), 'gene_cache')
# GENES_TO_CACHE = {l.strip('\n') for l in open(os.path.join(os.path.dirname(__file__), 'genes_to_cache.txt'))}

def connect_db():
    """
    Connects to the specific database.
    """
    client = pymongo.MongoClient(host=app.config['DB_HOST'], port=app.config['DB_PORT'])
    return client[app.config['DB_NAME']]


def parse_tabix_file_subset(tabix_filenames, subset_i, subset_n, record_parser, canonical_transcripts=None):
    """
    Returns a generator of parsed record objects (as returned by record_parser) for the i'th out n subset of records
    across all the given tabix_file(s). The records are split by files and contigs within files, with 1/n of all contigs
    from all files being assigned to this the i'th subset.

    Args:
        tabix_filenames: a list of one or more tabix-indexed files. These will be opened using pysam.Tabixfile
        subset_i: zero-based number
        subset_n: total number of subsets
        record_parser: a function that takes a file-like object and returns a generator of parsed records
    """
    start_time = time.time()
    open_tabix_files = [pysam.Tabixfile(tabix_filename) for tabix_filename in tabix_filenames]
    tabix_file_contig_pairs = [(tabix_file, contig) for tabix_file in open_tabix_files for contig in tabix_file.contigs]
    tabix_file_contig_subset = tabix_file_contig_pairs[subset_i : : subset_n]  # get every n'th tabix_file/contig pair
    short_filenames = ", ".join(map(os.path.basename, tabix_filenames))
    num_file_contig_pairs = len(tabix_file_contig_subset)
    print(("Loading subset %(subset_i)s of %(subset_n)s total: %(num_file_contig_pairs)s contigs from "
           "%(short_filenames)s") % locals())
    counter = 0
    for tabix_file, contig in tabix_file_contig_subset:
        header_iterator = tabix_file.header
        records_iterator = tabix_file.fetch(contig, 0, 10**9, multiple_iterators=True)
        for parsed_record in record_parser(itertools.chain(header_iterator, records_iterator), canonical_transcripts):
            counter += 1
            yield parsed_record

            if counter % 100000 == 0:
                seconds_elapsed = int(time.time()-start_time)
                print(("Loaded %(counter)s records from subset %(subset_i)s of %(subset_n)s from %(short_filenames)s "
                       "(%(seconds_elapsed)s seconds)") % locals())

    print("Finished loading subset %(subset_i)s from  %(short_filenames)s (%(counter)s records)" % locals())


def load_base_coverage(project_name=None, genome=None):
    def load_coverage(coverage_files, i, n, db):
        coverage_generator = parse_tabix_file_subset(coverage_files, i, n, get_base_coverage_from_file)
        try:
            db.base_coverage.insert(coverage_generator, w=0)
        except pymongo.errors.InvalidOperation:
            pass  # handle error when coverage_generator is empty

    full_db = get_db()
    if project_name:
        projects = [get_one_project(full_db, project_name, genome)]
    else:
        projects = get_projects(full_db)
    procs = []
    for project_name, genome in projects:
        db = full_db[get_project_key(project_name, genome)]
        db.base_coverage.drop()
        print("Dropped db.base_coverage for " + project_name)
        # load coverage first; variant info will depend on coverage
        db.base_coverage.ensure_index('xpos')

        coverage_files = glob.glob(app.config['BASE_COVERAGE_FILES'] % (genome, project_name))
        num_procs = app.config['LOAD_DB_PARALLEL_PROCESSES']
        # random.shuffle(app.config['BASE_COVERAGE_FILES'])
        max_procs = max(1, num_procs / len(projects))
        for i in range(max_procs):
            p = Process(target=load_coverage, args=(coverage_files, i, num_procs, db))
            p.start()
            procs.append(p)
    return procs

    #print 'Done loading coverage. Took %s seconds' % int(time.time() - start_time)


def get_projects(full_db):
    projects = full_db.projects.find()
    projects = [(project['name'], project['genome']) for project in projects]
    return projects


def get_one_project(full_db, project_name, genome):
    if not genome:
        print('Error! Please specify project name and genome')
        sys.exit(2)
    project = (project_name, genome)
    if not full_db.projects.find({'name': project_name, 'genome': genome}):
        full_db.projects.insert({'name': project_name, 'genome': genome})
    return project


def load_variants_file(project_name=None, genome=None):
    def load_variants(sites_file, i, n, db, canonical_transcripts):
        variants_generator = parse_tabix_file_subset([sites_file], i, n, get_variants_from_sites_vcf, canonical_transcripts)
        try:
            db.variants.insert(variants_generator, w=0)
        except pymongo.errors.InvalidOperation:
            pass  # handle error when variant_generator is empty

    full_db = get_db()
    if project_name:
        projects = [get_one_project(full_db, project_name, genome)]
    else:
        full_db.projects.drop()
        for genome in 'hg19', 'hg38':
            sites_vcfs = glob.glob(app.config['SITES_VCFS'] % (genome, '*'))
            project_names = [basename(vcf).split('.')[0] for vcf in sites_vcfs]
            if project_names:
                full_db.projects.insert({'name': project_name, 'genome': genome} for project_name in project_names)
        projects = get_projects(full_db)
    procs = []
    for project_name, genome in projects:
        db = full_db[get_project_key(project_name, genome)]
        db.variants.drop()
        print("Dropped db.variants for " + project_name)

        # grab variants from sites VCF
        db.variants.ensure_index('xpos')
        db.variants.ensure_index('xstart')
        db.variants.ensure_index('xstop')
        db.variants.ensure_index('rsid')
        db.variants.ensure_index('genes')
        db.variants.ensure_index('transcripts')

        sites_vcfs = glob.glob(app.config['SITES_VCFS'] % (genome, project_name))
        if len(sites_vcfs) == 0:
            raise IOError("No vcf file found")
        elif len(sites_vcfs) > 1:
            raise Exception("More than one sites vcf file found: %s" % sites_vcfs)

        canonical_transcripts = defaultdict()
        with gzip.open(app.config['CANONICAL_TRANSCRIPT_FILE'] % genome) as canonical_transcript_file:
            for gene, transcript in get_canonical_transcripts(canonical_transcript_file):
                canonical_transcripts[gene] = transcript

        num_procs = app.config['LOAD_DB_PARALLEL_PROCESSES']
        max_procs = max(1, num_procs / len(projects))
        for i in range(max_procs):
            p = Process(target=load_variants, args=(sites_vcfs[0], i, num_procs, db, canonical_transcripts))
            p.start()
            procs.append(p)
    return procs

    #print 'Done loading variants. Took %s seconds' % int(time.time() - start_time)

'''
def load_constraint_information():
    db = get_db()

    db.constraint.drop()
    print 'Dropped db.constraint.'

    start_time = time.time()

    with gzip.open(app.config['CONSTRAINT_FILE']) as constraint_file:
        for transcript in get_constraint_information(constraint_file):
            db.constraint.insert(transcript, w=0)

    db.constraint.ensure_index('transcript')
    print 'Done loading constraint info. Took %s seconds' % int(time.time() - start_time)


def load_mnps():
    db = get_db()
    start_time = time.time()

    db.variants.ensure_index('has_mnp')
    print 'Done indexing.'
    while db.variants.find_and_modify({'has_mnp' : True}, {'$unset': {'has_mnp': '', 'mnps': ''}}):
        pass
    print 'Deleted MNP data.'

    with gzip.open(app.config['MNP_FILE']) as mnp_file:
        for mnp in get_mnp_data(mnp_file):
            variant = lookups.get_raw_variant(db, mnp['xpos'], mnp['ref'], mnp['alt'], True)
            db.variants.find_and_modify({'_id': variant['_id']}, {'$set': {'has_mnp': True}, '$push': {'mnps': mnp}}, w=0)

    db.variants.ensure_index('has_mnp')
    print 'Done loading MNP info. Took %s seconds' % int(time.time() - start_time)'''


def load_gene_models():
    full_db = get_db()
    for genome in 'hg19', 'hg38':
        db = full_db[genome]
        db.genes.drop()
        db.transcripts.drop()
        db.exons.drop()
        print 'Dropped db.genes, db.transcripts, and db.exons for ' + genome

        start_time = time.time()

        canonical_transcripts = {}
        with gzip.open(app.config['CANONICAL_TRANSCRIPT_FILE'] % genome) as canonical_transcript_file:
            for gene, transcript in get_canonical_transcripts(canonical_transcript_file):
                canonical_transcripts[gene] = transcript

        omim_annotations = {}
        '''with gzip.open(app.config['OMIM_FILE'] % genome) as omim_file:
            for fields in get_omim_associations(omim_file):
                if fields is None:
                    continue
                gene, transcript, accession, description = fields
                omim_annotations[gene] = (accession, description)

        dbnsfp_info = {}
        with gzip.open(app.config['DBNSFP_FILE'] % genome) as dbnsfp_file:
            for dbnsfp_gene in get_dbnsfp_info(dbnsfp_file):
                other_names = [other_name.upper() for other_name in dbnsfp_gene['gene_other_names']]
                dbnsfp_info[dbnsfp_gene['ensembl_gene']] = (dbnsfp_gene['gene_full_name'], other_names)'''

        print 'Done loading metadata. Took %s seconds' % int(time.time() - start_time)

        # grab genes from GTF
        start_time = time.time()
        with gzip.open(app.config['FEATURES_FILE'] % genome) as features_file:
            for gene in get_genes_from_features(features_file):
                gene_id = gene['gene_id']
                if gene_id in canonical_transcripts:
                    gene['canonical_transcript'] = canonical_transcripts[gene_id]
                if gene_id in omim_annotations:
                    gene['omim_accession'] = omim_annotations[gene_id][0]
                    gene['omim_description'] = omim_annotations[gene_id][1]
                '''if gene_id in dbnsfp_info:
                    gene['full_gene_name'] = dbnsfp_info[gene_id][0]
                    gene['other_names'] = dbnsfp_info[gene_id][1]'''
                db.genes.insert(gene, w=0)

        print 'Done loading genes. Took %s seconds' % int(time.time() - start_time)

        start_time = time.time()
        db.genes.ensure_index('gene_id')
        db.genes.ensure_index('gene_name_upper')
        db.genes.ensure_index('gene_name')
        db.genes.ensure_index('other_names')
        db.genes.ensure_index('xstart')
        db.genes.ensure_index('xstop')
        print 'Done indexing gene table. Took %s seconds' % int(time.time() - start_time)

        # and now transcripts
        start_time = time.time()
        with gzip.open(app.config['FEATURES_FILE'] % genome) as features_file:
            db.transcripts.insert((transcript for transcript in get_transcripts_from_features(features_file)), w=0)
        print 'Done loading transcripts. Took %s seconds' % int(time.time() - start_time)

        start_time = time.time()
        db.transcripts.ensure_index('transcript_id')
        db.transcripts.ensure_index('gene_id')
        print 'Done indexing transcript table. Took %s seconds' % int(time.time() - start_time)

        # Building up gene definitions
        start_time = time.time()
        with gzip.open(app.config['FEATURES_FILE'] % genome) as features_file:
            db.exons.insert((exon for exon in get_exons_from_features(features_file)), w=0)
        print 'Done loading exons. Took %s seconds' % int(time.time() - start_time)

        start_time = time.time()
        db.exons.ensure_index('exon_id')
        db.exons.ensure_index('transcript_id')
        db.exons.ensure_index('gene_id')
        print 'Done indexing exon table. Took %s seconds' % int(time.time() - start_time)

    return []


def load_dbsnp_file():
    db = get_db()

    def load_dbsnp(dbsnp_file, i, n, db):
        if os.path.isfile(dbsnp_file + ".tbi"):
            dbsnp_record_generator = parse_tabix_file_subset([dbsnp_file], i, n, get_snp_from_dbsnp_file)
            try:
                db.dbsnp.insert(dbsnp_record_generator, w=0)
            except pymongo.errors.InvalidOperation:
                pass  # handle error when coverage_generator is empty

        else:
            with gzip.open(dbsnp_file) as f:
                db.dbsnp.insert((snp for snp in get_snp_from_dbsnp_file(f)), w=0)

    db.dbsnp.drop()
    db.dbsnp.ensure_index('rsid')
    db.dbsnp.ensure_index('xpos')
    start_time = time.time()
    dbsnp_file = app.config['DBSNP_FILE']

    print "Loading dbsnp from %s" % dbsnp_file
    if os.path.isfile(dbsnp_file + ".tbi"):
        num_procs = app.config['LOAD_DB_PARALLEL_PROCESSES']
    else:
        # see if non-tabixed .gz version exists
        if os.path.isfile(dbsnp_file):
            print(("WARNING: %(dbsnp_file)s.tbi index file not found. Will use single thread to load dbsnp."
                "To create a tabix-indexed dbsnp file based on UCSC dbsnp, do: \n"
                "   wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/snp141.txt.gz \n"
                "   gzcat snp141.txt.gz | cut -f 1-5 | bgzip -c > snp141.txt.bgz \n"
                "   tabix -0 -s 2 -b 3 -e 4 snp141.txt.bgz") % locals())
            num_procs = 1
        else:
            raise Exception("dbsnp file %s(dbsnp_file)s not found." % locals())

    procs = []
    for i in range(num_procs):
        p = Process(target=load_dbsnp, args=(dbsnp_file, i, num_procs, db))
        p.start()
        procs.append(p)

    return procs
    #print 'Done loading dbSNP. Took %s seconds' % int(time.time() - start_time)

    #start_time = time.time()
    #db.dbsnp.ensure_index('rsid')
    #print 'Done indexing dbSNP table. Took %s seconds' % int(time.time() - start_time)


def load_db():
    """
    Load the database
    """
    # Initialize database
    # Don't need to explicitly create tables with mongo, just indices
    confirm = raw_input('This will drop the database and reload. Are you sure you want to continue? [no] ')
    if not confirm.startswith('y'):
        print('Exiting...')
        sys.exit(1)
    all_procs = []
    for load_function in [load_variants_file, load_base_coverage, load_gene_models]:
        procs = load_function()
        all_procs.extend(procs)
        print("Started %s processes to run %s" % (len(procs), load_function.__name__))

    [p.join() for p in all_procs]
    # print('Done! Loading MNPs...')
    # load_mnps()
    print('Done! Creating cache...')
    create_cache()
    print('Done!')


def delete_project(project_name, genome, silent=False):
    full_db = get_db()
    full_db[get_project_key(project_name, genome)].drop()
    full_db.projects.remove({'name': project_name, 'genome': genome})
    if not silent:
        print(project_name + ' was deleted from the database')
        print('Refreshing cache...')
    create_cache()
    if not silent:
        print('Done!')


def add_project(project_name, genome):
    delete_project(project_name, genome, silent=True)
    all_procs = []
    full_db = get_db()
    full_db.projects.insert({'name': project_name, 'genome': genome})
    print('Adding ' + project_name + ' to the database')
    for load_function in [load_variants_file, load_base_coverage]:
        procs = load_function(project_name, genome)
        all_procs.extend(procs)
        print("Started %s processes to run %s" % (len(procs), load_function.__name__))

    [p.join() for p in all_procs]
    print('Done! Creating cache...')
    create_cache()
    precalculate_metrics(project_name=project_name, genome=genome)
    print('Done!')


def save_autocomplete_data(autocomplete_strings, output_fpath):
    f = open(os.path.join(os.path.dirname(__file__), output_fpath), 'w')
    for s in sorted(autocomplete_strings):
        f.write(s+'\n')
    f.close()


def create_cache():
    """
    This is essentially a compile step that generates all cached resources.
    Creates files like autocomplete_entries.txt
    Should be run on every redeploy.
    """
    # create autocomplete_entries.txt
    autocomplete_strings = []
    full_db = get_db()
    projects = get_projects(full_db)
    for genome in 'hg19', 'hg38':
        for gene in full_db[genome].genes.find():
            autocomplete_strings.append(gene['gene_name'])
            if 'other_names' in gene:
                autocomplete_strings.extend(gene['other_names'])

        save_autocomplete_data(set(autocomplete_strings), genome + '_autocomplete_strings.txt')
    save_autocomplete_data([project_name for (project_name, project_genome) in projects], 'autocomplete_projects.txt')

    '''
    # create static gene pages for genes in
    if not os.path.exists(GENE_CACHE_DIR):
        os.makedirs(GENE_CACHE_DIR)

    # get list of genes ordered by num_variants
    for gene_id in GENES_TO_CACHE:
        try:
            page_content = get_gene_page_content(gene_id)
        except Exception as e:
            print e
            continue
        f = open(os.path.join(GENE_CACHE_DIR, '{}.html'.format(gene_id)), 'w')
        f.write(page_content)
        f.close()'''


def precalculate_metrics(project_name=None, genome=None):
    import numpy
    full_db = get_db()
    if project_name:
        projects = [get_one_project(full_db, project_name, genome)]
    else:
        projects = get_projects(full_db)
    for project_name, genome in projects:
        db = full_db[get_project_key(project_name, genome)]
        print 'Reading %s variants for %s...' % (db.variants.count(), project_name)
        metrics = defaultdict(list)
        binned_metrics = defaultdict(list)
        progress = 0
        start_time = time.time()
        for variant in db.variants.find(projection=['quality_metrics', 'site_quality', 'allele_num', 'allele_count']):
            for metric, value in variant['quality_metrics'].iteritems():
                metrics[metric].append(float(value))
            qual = float(variant['site_quality'])
            metrics['site_quality'].append(qual)
            if variant['allele_num'] == 0: continue
            if variant['allele_count'] == 1:
                binned_metrics['singleton'].append(qual)
            elif variant['allele_count'] == 2:
                binned_metrics['doubleton'].append(qual)
            else:
                for af in AF_BUCKETS:
                    if float(variant['allele_count'])/variant['allele_num'] < af:
                        binned_metrics[af].append(qual)
                        break
            progress += 1
            if not progress % 100000:
                print 'Read %s variants. Took %s seconds' % (progress, int(time.time() - start_time))
        print 'Done reading variants. Dropping metrics database... '
        db.metrics.drop()
        print 'Dropped metrics database. Calculating metrics...'
        for metric in metrics:
            bin_range = None
            data = map(numpy.log, metrics[metric]) if metric == 'DP' else metrics[metric]
            if metric == 'FS':
                bin_range = (0, 20)
            elif metric == 'VQSLOD':
                bin_range = (-20, 20)
            elif metric == 'InbreedingCoeff':
                bin_range = (0, 1)
            if bin_range is not None:
                data = [x if (x > bin_range[0]) else bin_range[0] for x in data]
                data = [x if (x < bin_range[1]) else bin_range[1] for x in data]
            hist = numpy.histogram(data, bins=40, range=bin_range)
            edges = hist[1]
            # mids = [(edges[i]+edges[i+1])/2 for i in range(len(edges)-1)]
            lefts = [edges[i] for i in range(len(edges)-1)]
            db.metrics.insert({
                'metric': metric,
                'mids': lefts,
                'hist': list(hist[0])
            })
        for metric in binned_metrics:
            hist = numpy.histogram(map(numpy.log, binned_metrics[metric]), bins=40)
            edges = hist[1]
            mids = [(edges[i]+edges[i+1])/2 for i in range(len(edges)-1)]
            db.metrics.insert({
                'metric': 'binned_%s' % metric,
                'mids': mids,
                'hist': list(hist[0])
            })
        db.metrics.ensure_index('metric')
    print 'Done pre-calculating metrics!'


def get_db():
    """
    Opens a new database connection if there is none yet for the
    current application context.
    """
    if not hasattr(g, 'db_conn'):
        g.db_conn = connect_db()
    return g.db_conn


# @app.teardown_appcontext
# def close_db(error):
#     """Closes the database again at the end of the request."""
#     if hasattr(g, 'db_conn'):
#         g.db_conn.close()


@app.route('/')
def homepage():
    cache_key = 't-homepage'
    t = cache.get(cache_key)
    if t is None:
        t = render_template('homepage.html')
        cache.set(cache_key, t)
    return t


@app.route('/<project_genome>/<project_name>/')
def project_page(project_name, project_genome):
    t = render_template(
        'project_page.html',
        project_name=project_name,
        genome=project_genome
    )
    return t


@app.route('/<project_genome>/<project_name>/autocomplete/<query>')
def awesome_autocomplete(query, project_name, project_genome):
    if not hasattr(g, 'autocomplete_strings'):
        g.autocomplete_strings = dict()
        for genome in 'hg19', 'hg38':
            g.autocomplete_strings[genome] = [s.strip() for s in open(os.path.join(os.path.dirname(__file__), genome + '_autocomplete_strings.txt'))]
    suggestions = lookups.get_awesomebar_suggestions(g.autocomplete_strings[project_genome], query)
    return Response(json.dumps([{'value': s} for s in suggestions]),  mimetype='application/json')


@app.route('/autocomplete_project/')
def show_all_projects():
    if not hasattr(g, 'autocomplete_projects'):
        g.autocomplete_projects = [s.strip() for s in open(os.path.join(os.path.dirname(__file__), 'autocomplete_projects.txt'))]
    suggestions = g.autocomplete_projects
    return Response(json.dumps([{'value': s} for s in suggestions]),  mimetype='application/json')


@app.route('/autocomplete_project/<query>')
def awesome_project_autocomplete(query):
    if not hasattr(g, 'autocomplete_projects'):
        g.autocomplete_projects = [s.strip() for s in open(os.path.join(os.path.dirname(__file__), 'autocomplete_projects.txt'))]
    suggestions = lookups.get_awesomebar_suggestions(g.autocomplete_projects, query)
    return Response(json.dumps([{'value': s} for s in suggestions]),  mimetype='application/json')


@app.route('/<project_genome>/<project_name>/awesome')
def awesome(project_name, project_genome):
    db = get_db()
    query = request.args.get('query')
    datatype, identifier = lookups.get_awesomebar_result(db, project_name, project_genome, query)

    print "Searched for %s: %s" % (datatype, identifier)
    if datatype == 'gene':
        return redirect('/{}/{}/gene/{}'.format(project_genome, project_name, identifier))
    elif datatype == 'transcript':
        return redirect('/{}/{}/transcript/{}'.format(project_genome, project_name, identifier))
    elif datatype == 'variant':
        return redirect('/{}/{}/variant/{}'.format(project_genome, project_name, identifier))
    elif datatype == 'region':
        return redirect('/{}/{}/region/{}'.format(project_genome, project_name, identifier))
    elif datatype == 'dbsnp_variant_set':
        return redirect('/{}/{}/dbsnp/{}'.format(project_genome, project_name, identifier))
    elif datatype == 'error':
        return redirect('/{}/{}/error/{}'.format(project_genome, project_name, identifier))
    elif datatype == 'not_found':
        return redirect('/{}/{}/not_found/{}'.format(project_genome, project_name, identifier))
    else:
        raise Exception


@app.route('/awesomeproject')
def awesome_project():
    db = get_db()
    query = request.args.get('query')
    project = lookups.get_project_by_project_name(db, query)
    if project:
        return redirect('/{}/{}'.format(project['genome'], project['name']))

    return redirect('/not_found/{}'.format(query))


@app.route('/<project_genome>/<project_name>/variant/<variant_str>')
def variant_page(project_name, project_genome, variant_str):
    db = get_db()
    try:
        chrom, pos, ref, alt = variant_str.split('-')
        pos = int(pos)
        # pos, ref, alt = get_minimal_representation(pos, ref, alt)
        xpos = get_xpos(chrom, pos)
        variant = lookups.get_variant(db, project_name, project_genome, xpos, ref, alt)

        if variant is None:
            variant = {
                'chrom': chrom,
                'pos': pos,
                'xpos': xpos,
                'ref': ref,
                'alt': alt
            }
        consequences = OrderedDict()
        if 'vep_annotations' in variant:
            add_consequence_to_variant(variant)
            variant['vep_annotations'] = remove_extraneous_vep_annotations(variant['vep_annotations'])
            variant['vep_annotations'] = order_vep_by_csq(variant['vep_annotations'])  # Adds major_consequence
            for annotation in variant['vep_annotations']:
                annotation['HGVS'] = get_proper_hgvs(annotation)
                consequences.setdefault(annotation['major_consequence'], {}).setdefault(annotation['Gene_Name'], []).append(annotation)
        base_coverage = lookups.get_coverage_for_bases(db, project_name, project_genome, xpos, xpos + len(ref) - 1)
        any_covered = any([x['has_coverage'] for x in base_coverage])
        metrics = lookups.get_metrics(db, project_name, project_genome, variant)

        igv_genes_path = ''
        if project_genome == 'hg19':
            igv_genes_path = '//igv.broadinstitute.org/annotations/hg19/genes/gencode.v18.collapsed.bed'
            igv_genes_index = igv_genes_path + '.idx'
        elif project_genome == 'hg38':
            igv_genes_path = '//igv.broadinstitute.org/annotations/hg38/genes/gencode.v24.annotation.sorted.gtf.gz'
            igv_genes_index = igv_genes_path + '.tbi'

        # check the appropriate sqlite db to get the *expected* number of
        # available bams and *actual* number of available bams for this variant
        sqlite_db_path = os.path.join(
            app.config["READ_VIZ_DIR"],
            "combined_bams",
            chrom,
            "combined_chr%s_%03d.db" % (chrom, pos % 1000))
        print(sqlite_db_path)
        try:
            read_viz_db = sqlite3.connect(sqlite_db_path)
            n_het = read_viz_db.execute("select n_expected_samples, n_available_samples from t "
                "where chrom=? and pos=? and ref=? and alt=? and het_or_hom=?", (chrom, pos, ref, alt, 'het')).fetchone()
            n_hom = read_viz_db.execute("select n_expected_samples, n_available_samples from t "
                "where chrom=? and pos=? and ref=? and alt=? and het_or_hom=?", (chrom, pos, ref, alt, 'hom')).fetchone()
            read_viz_db.close()
        except Exception, e:
            logging.debug("Error when accessing sqlite db: %s - %s", sqlite_db_path, e)
            n_het = n_hom = None

        read_viz_dict = {
            'het': {'n_expected': n_het[0] if n_het is not None and n_het[0] is not None else -1, 'n_available': n_het[1] if n_het and n_het[1] else 0,},
            'hom': {'n_expected': n_hom[0] if n_hom is not None and n_hom[0] is not None else -1, 'n_available': n_hom[1] if n_hom and n_hom[1] else 0,},
        }

        for het_or_hom in ('het', 'hom',):
            #read_viz_dict[het_or_hom]['some_samples_missing'] = (read_viz_dict[het_or_hom]['n_expected'] > 0)    and (read_viz_dict[het_or_hom]['n_expected'] - read_viz_dict[het_or_hom]['n_available'] > 0)
            read_viz_dict[het_or_hom]['all_samples_missing'] = (read_viz_dict[het_or_hom]['n_expected'] != 0) and (read_viz_dict[het_or_hom]['n_available'] == 0)
            read_viz_dict[het_or_hom]['readgroups'] = [
                '%(chrom)s-%(pos)s-%(ref)s-%(alt)s_%(het_or_hom)s%(i)s' % locals()
                    for i in range(read_viz_dict[het_or_hom]['n_available'])
            ]   #eg. '1-157768000-G-C_hom1',

            read_viz_dict[het_or_hom]['urls'] = [
                os.path.join('combined_bams', chrom, 'combined_chr%s_%03d.bam' % (chrom, pos % 1000))
                    for i in range(read_viz_dict[het_or_hom]['n_available'])
            ]
        chrom = chrom.replace('chr', '')
        bam_fpath = os.path.join(
            app.config["READ_VIZ_DIR_HTML"] % (project_genome, project_name),
            "%s-" % chrom)

        read_group = None
        sample_names = []
        if 'sample_names' in variant:
            sample_names = [sample_name.replace('-', '_') for idx, sample_name in enumerate(variant['sample_names'])
                            if variant['sample_data'][idx]]

        if 'transcripts' in variant and variant['transcripts']:
            for transcript_id in variant['transcripts']:
                transcript_exons = lookups.get_exons_in_transcript(db, project_genome, transcript_id)
                for idx, exon in enumerate(transcript_exons):
                    if exon['start'] <= pos <= exon['stop']:
                        start, end = exon['start'], exon['stop']
                        transcript = sorted(variant['transcripts'])[0]
                        read_group = '{chrom}-{transcript}-{idx}-'.format(**locals())
                        break
                if read_group:
                    break
        if not read_group:
            read_group = '%(chrom)s-%(pos)s-%(ref)s-%(alt)s-' % locals()
        print 'Rendering variant: %s' % variant_str
        return render_template(
            'variant.html',
            project_name=project_name,
            genome=project_genome,
            variant=variant,
            variant_json=JSONEncoder().encode(variant),
            base_coverage=base_coverage,
            base_coverage_json=JSONEncoder().encode(base_coverage),
            consequences=consequences,
            consequences_json=JSONEncoder().encode(consequences),
            any_covered=any_covered,
            any_covered_json=JSONEncoder().encode(any_covered),
            metrics=metrics,
            metrics_json=JSONEncoder().encode(metrics),
            igv_genes_url=igv_genes_path,
            igv_genes_index=igv_genes_index,
            read_viz=read_viz_dict,
            bam_url=bam_fpath,
            sample_names=sample_names,
            read_group=read_group
        )
    except Exception:
        print 'Failed on variant:', variant_str, ';Error=', traceback.format_exc()
        abort(404)


@app.route('/<project_genome>/<project_name>/gene/<gene_id>')
def gene_page(project_name, project_genome, gene_id):
    # if gene_id in GENES_TO_CACHE:
    #    return open(os.path.join(GENE_CACHE_DIR, '{}.html'.format(gene_id))).read()
    # else:
    return get_gene_page_content(project_name, project_genome, gene_id)


def get_gene_page_content(project_name, project_genome, gene_id):
    db = get_db()
    try:
        gene = lookups.get_gene(db, project_genome, gene_id)
        if gene is None:
            abort(404)
        cache_key = 't-gene-{}-{}-{}'.format(project_genome, project_name, gene_id)
        t = cache.get(cache_key)
        print 'Rendering %sgene: %s' % ('' if t is None else 'cached ', gene_id)
        if t is None:
            variants_in_gene = lookups.get_variants_in_gene(db, project_name, project_genome, gene_id)
            transcripts_in_gene = lookups.get_transcripts_in_gene(db, project_genome, gene_id)

            # Get some canonical transcript and corresponding info
            transcript_id = gene['canonical_transcript']
            transcript = lookups.get_transcript(db, project_genome, transcript_id)
            variants_in_transcript = lookups.get_variants_in_transcript(db, project_name, project_genome, transcript_id)
            coverage_stats = lookups.get_coverage_for_transcript(db, project_name, project_genome, transcript['xstart'] - EXON_PADDING, transcript['xstop'] + EXON_PADDING)
            add_transcript_coordinate_to_variants(db, project_genome, variants_in_transcript, transcript_id)
            constraint_info = lookups.get_constraint_for_transcript(db, project_genome, transcript_id)

            t = render_template(
                'gene.html',
                project_name=project_name,
                genome=project_genome,
                gene=gene,
                gene_json=JSONEncoder().encode(gene),
                transcript=transcript,
                transcript_json=JSONEncoder().encode(transcript),
                variants_in_gene=variants_in_gene,
                variants_in_gene_json=JSONEncoder().encode(variants_in_gene),
                variants_in_transcript=variants_in_transcript,
                variants_in_transcript_json=JSONEncoder().encode(variants_in_transcript),
                transcripts_in_gene=transcripts_in_gene,
                transcripts_in_gene_json=JSONEncoder().encode(transcripts_in_gene),
                coverage_stats=coverage_stats,
                coverage_stats_json=JSONEncoder().encode(coverage_stats),
                constraint=constraint_info,
                csq_order=csq_order,
                csq_order_json=JSONEncoder().encode(csq_order)
            )
            cache.set(cache_key, t)
        return t
    except Exception, e:
        print 'Failed on gene:', gene_id, ';Error=', traceback.format_exc()
        abort(404)


@app.route('/<project_genome>/<project_name>/transcript/<transcript_id>')
def transcript_page(project_name, project_genome, transcript_id):
    db = get_db()
    try:
        transcript = lookups.get_transcript(db, project_genome, transcript_id)

        cache_key = 't-transcript-{}-{}-{}'.format(project_genome, project_name, transcript_id)
        t = cache.get(cache_key)
        print 'Rendering %stranscript: %s' % ('' if t is None else 'cached ', transcript_id)
        if t is None:

            gene = lookups.get_gene(db, project_genome, transcript['gene_id'])
            gene['transcripts'] = lookups.get_transcripts_in_gene(db, project_genome, transcript['gene_id'])
            variants_in_transcript = lookups.get_variants_in_transcript(db, project_name, project_genome, transcript_id)

            coverage_stats = lookups.get_coverage_for_transcript(db, project_name, project_genome, transcript['xstart'] - EXON_PADDING, transcript['xstop'] + EXON_PADDING)

            add_transcript_coordinate_to_variants(db, project_genome, variants_in_transcript, transcript_id)

            t = render_template(
                'transcript.html',
                project_name=project_name,
                genome=project_genome,
                transcript=transcript,
                transcript_json=JSONEncoder().encode(transcript),
                variants_in_transcript=variants_in_transcript,
                variants_in_transcript_json=JSONEncoder().encode(variants_in_transcript),
                coverage_stats=coverage_stats,
                coverage_stats_json=JSONEncoder().encode(coverage_stats),
                gene=gene,
                gene_json=JSONEncoder().encode(gene),
                csq_order=csq_order,
                csq_order_json=JSONEncoder().encode(csq_order)
            )
            cache.set(cache_key, t)
        return t
    except Exception, e:
        print 'Failed on transcript:', transcript_id, ';Error=', traceback.format_exc()
        abort(404)


@app.route('/<project_genome>/<project_name>/region/<region_id>')
def region_page(project_name, project_genome, region_id):
    db = get_db()
    try:
        region = region_id.split('-')
        cache_key = 't-region-{}-{}-{}'.format(project_genome, project_name, region_id)
        t = cache.get(cache_key)
        print 'Rendering %sregion: %s' % ('' if t is None else 'cached ', region_id)
        if t is None:
            chrom = region[0]
            start = None
            stop = None
            if len(region) == 3:
                chrom, start, stop = region
                start = int(start)
                stop = int(stop)
            if start is None or stop - start > REGION_LIMIT or stop < start:
                return render_template(
                    'region.html',
                    project_name=project_name,
                    genome=project_genome,
                    genes_in_region=None,
                    variants_in_region=None,
                    chrom=chrom,
                    start=start,
                    stop=stop,
                    coverage=None,
                    csq_order=csq_order,
                    csq_order_json=JSONEncoder().encode(csq_order)
                )
            if start == stop:
                start -= 20
                stop += 20
            genes_in_region = lookups.get_genes_in_region(db, project_genome, chrom, start, stop)
            variants_in_region = lookups.get_variants_in_region(db, project_name, project_genome, chrom, start, stop)
            xstart = get_xpos(chrom, start)
            xstop = get_xpos(chrom, stop)
            coverage_array = lookups.get_coverage_for_bases(db, project_name, project_genome, xstart, xstop)
            t = render_template(
                'region.html',
                project_name=project_name,
                genome=project_genome,
                genes_in_region=genes_in_region,
                genes_in_region_json=JSONEncoder().encode(genes_in_region),
                variants_in_region=variants_in_region,
                variants_in_region_json=JSONEncoder().encode(variants_in_region),
                chrom=chrom,
                start=start,
                stop=stop,
                coverage=coverage_array,
                coverage_json=JSONEncoder().encode(coverage_array),
                csq_order=csq_order,
                csq_order_json=JSONEncoder().encode(csq_order)
            )
            cache.set(cache_key, t)
        return t
    except Exception, e:
        print 'Failed on region:', region_id, ';Error=', traceback.format_exc()
        abort(404)


@app.route('/<genome>/<project_name>/dbsnp/<rsid>')
def dbsnp_page(rsid, project_name, genome):
    db = get_db()
    try:
        variants = lookups.get_variants_by_rsid(db, project_name, genome, rsid)
        chrom = None
        start = None
        stop = None
        print 'Rendering rsid: %s' % rsid
        return render_template(
            'region.html',
            rsid=rsid,
            genome=genome,
            variants_in_region=variants,
            chrom=chrom,
            start=start,
            stop=stop,
            coverage=None,
            genes_in_region=None,
            csq_order=csq_order,
            csq_order_json=JSONEncoder().encode(csq_order)
        )
    except Exception, e:
        print 'Failed on rsid:', rsid, ';Error=', traceback.format_exc()
        abort(404)


@app.route('/not_found/<query>')
@app.errorhandler(404)
def not_found_project_page(query):
    return render_template(
        'not_found.html',
        query=query
    ), 404


@app.route('/<project_genome>/<project_name>/not_found/<query>')
@app.errorhandler(404)
def not_found_page(project_genome, project_name, query):
    return render_template(
        'not_found.html',
        genome=project_genome,
        project_name=project_name,
        query=query
    ), 404


@app.route('/error/<query>')
@app.errorhandler(404)
def error_page(query):
    return render_template(
        'error.html',
        query=query
    ), 404


@app.route('/downloads')
def downloads_page():
    return render_template('downloads.html')


@app.route('/about')
def about_page():
    return render_template('about.html')


@app.route('/participants')
def participants_page():
    return render_template('about.html')


@app.route('/terms')
def terms_page():
    return render_template('terms.html')


@app.route('/contact')
def contact_page():
    return render_template('contact.html')


@app.route('/faq')
def faq_page():
    return render_template('faq.html')


@app.route('/text')
def text_page():
    db = get_db()
    query = request.args.get('text')
    datatype, identifier = lookups.get_awesomebar_result(db, query)
    if datatype in ['gene', 'transcript']:
        gene = lookups.get_gene(db, identifier)
        link = "http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=chr%(chrom)s%%3A%(start)s-%(stop)s" % gene
        output = '''Searched for %s. Found %s.
%s; Canonical: %s.
%s''' % (query, identifier, gene['full_gene_name'], gene['canonical_transcript'], link)
        output += '' if 'omim_accession' not in gene else '''
In OMIM: %(omim_description)s
http://omim.org/entry/%(omim_accession)s''' % gene
        return output
    elif datatype == 'error' or datatype == 'not_found':
        return "Gene/transcript %s not found" % query
    else:
        return "Search types other than gene transcript not yet supported"


@app.route('/<project_genome>/combined_bams/<project_name>/<bam>')
def read_viz_files(project_genome, project_name, bam):
    full_path = os.path.abspath(os.path.join(app.config["READ_VIZ_DIR"] % (project_genome, project_name), bam))

    # security check - only files under READ_VIZ_DIR should be accsessible
    # if not full_path.startswith(app.config["READ_VIZ_DIR"]):
    #     return "Invalid path: %s" % path

    logging.info("path: " + full_path)

    # handle igv.js Range header which it uses to request a subset of a .bam
    range_header = request.headers.get('Range', None)
    if not range_header:
        return send_from_directory(app.config["READ_VIZ_DIR"] % (project_genome, project_name), bam)

    m = re.search('(\d+)-(\d*)', range_header)
    if not m:
        error_msg = "ERROR: unexpected range header syntax: %s" % range_header
        logging.error(error_msg)
        return error_msg

    size = os.path.getsize(full_path)
    offset = int(m.group(1))
    length = int(m.group(2) or size) - offset

    data = None
    with open(full_path, 'rb') as f:
        f.seek(offset)
        data = f.read(length)

    rv = Response(data, 206, mimetype="application/octet-stream", direct_passthrough=True)
    rv.headers.add('Content-Range', 'bytes {0}-{1}/{2}'.format(offset, offset + length - 1, size))

    logging.info("GET range request: %s-%s %s" % (m.group(1), m.group(2), full_path))
    return rv


@app.after_request
def apply_caching(response):
    # prevent click-jacking vulnerability identified by BITs
    response.headers["X-Frame-Options"] = "SAMEORIGIN"
    return response


if __name__ == "__main__":
    # adds Flask command line options for setting host, port, etc.
    app.run(host='0.0.0.0')
    # runner = Runner(app)
    # runner.run()
    # manager = Manager(app)
    # manager.run(host="0.0.0.0")
    # manager.add_command("runserver", Server(
    # use_debugger = True,
    # use_reloader = True,
    # host = '0.0.0.0') )
