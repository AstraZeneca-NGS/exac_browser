#!/usr/bin/env python2
import sys
from flask.ext.script import Manager
from exac import app
import exac

manager = Manager(app)


@manager.command
def hello():
    print "hello"


@manager.command
def load_db():
    exac.load_db()


@manager.command
def load_base_coverage(project_name):
    exac.load_base_coverage(project_name)


@manager.command
def load_variants_file(project_name):
    exac.load_variants_file(project_name)


@manager.command
def reload_variants(project_name):
    exac.load_variants_file(project_name)
    # exac.load_mnps()
    exac.precalculate_metrics(project_name)


@manager.command
def load_gene_models():
    exac.load_gene_models()


@manager.command
def load_dbsnp_file():
    exac.load_dbsnp_file()


@manager.command
def load_constraint_information():
    exac.load_constraint_information()


@manager.command
def load_mnps():
    exac.load_mnps()


@manager.command
def create_cache():
    exac.create_cache()


@manager.command
def precalculate_metrics(project_name=None):
    exac.precalculate_metrics(project_name)

if __name__ == "__main__":
    manager.run()

