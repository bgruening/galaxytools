#!/usr/bin/env python

__author__ = "Bjoern Gruening"
__version__ = "1.0"
__date__ = "2012"
__license__ = "GLP3+"

"""
    Creates a GenericMolecule Database, with tables and schema
"""
import argparse
import datetime
import sys

import cheminfolib
import generic as config
import psycopg2
from sqlalchemy import MetaData, create_engine
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import backref, relation
from sqlalchemy.schema import (Column, ForeignKeyConstraint,
                               PrimaryKeyConstraint)
from sqlalchemy.types import *


class molecule_type(UserDefinedType):
    """
    Define the new molecule Type in PG in sqlalchemy
    """

    def __init__(self):
        pass

    def get_col_spec(self):
        return "molecule"

    def bind_processor(self, dialect):
        def process(value):
            return value

        return process

    def result_processor(self, dialect, coltype):
        def process(value):
            return value

        return process


Base = declarative_base()

Base.metadata = MetaData(schema=config.schema)


class Molecule(Base):
    __tablename__ = "tbl_molecule"

    id = Column(
        Integer, primary_key=True
    )  # unique_id will be auto-generated, and linked to its chembl code in the synonyms table
    mol = Column(molecule_type)
    set_bits = Column(Integer, index=True)
    # maccs       = Column(VARBIT)
    # fp3         = Column(VARBIT)
    # fp4         = Column(VARBIT)
    inchi = Column(Text)
    inchi_key_first = Column(CHAR(14), index=True)
    inchi_key_last = Column(CHAR(14), index=True)
    can_smiles = Column(Text, index=True)
    spectrophore = Column(Text)
    molwt = Column(Numeric, index=True)
    logp = Column(Numeric, index=True)
    mr = Column(Numeric, index=True)
    hba = Column(Integer, index=True)
    hbd = Column(Integer, index=True)
    psa = Column(Numeric, index=True)
    rotbonds = Column(Integer, index=True)
    rings = Column(Integer, index=True)
    atoms = Column(Integer, index=True)

    def __init__(self, molecule):
        prop = cheminfolib.get_properties_ext(molecule)

        self.mol = molecule.write("mol")
        self.set_bits = len(molecule.calcfp().bits)
        """
        maccs = ['0']*166
        for i in molecule.calcfp('MACCS').bits:
            maccs[i-1] = '1'
        self.maccs = ''.join(maccs)

        fp3 = ['0']*1024
        for i in molecule.calcfp('FP3').bits:
            fp3[i-1] = '1'
        self.fp3 = ''.join(fp3)
        
        fp4 = ['0']*1024
        for i in molecule.calcfp('FP4').bits:
            fp4[i-1] = '1'
        self.fp4 = ''.join(fp4)
        """
        self.inchi = prop["inchi"]
        self.inchi_key_first = prop["inchi_key"][:14]
        self.inchi_key_last = prop["inchi_key"][15:]
        self.can_smiles = prop["can"]
        self.spectrophore = prop["spectrophore"]
        self.molwt = prop["molwt"]
        self.logp = prop["logp"]
        self.mr = prop["mr"]
        self.hbd = prop["donors"]
        self.hba = prop["acceptors"]
        self.psa = prop["psa"]
        self.rotbonds = prop["rotbonds"]
        self.rings = prop["rings"]
        self.atoms = prop["atoms"]

    def __repr__(self):
        return "Molecule(%s)" % (self.can_smiles)

    # __table_args__  = (
    #    {'schema': Base.metadata.schema}
    # )


class Synonym(Base):
    __tablename__ = "tbl_molecule_synonym"

    fk_molecule_id = Column(Integer, index=True)
    synonym = Column(String, nullable=False, index=True)
    # Store the IUPAC chemical name if available
    note = Column(CHAR(10), index=True)

    __table_args__ = (
        ForeignKeyConstraint(
            ["fk_molecule_id"],
            ["tbl_molecule.id"],
            onupdate="CASCADE",
            ondelete="CASCADE",
        ),
        PrimaryKeyConstraint("fk_molecule_id", "synonym"),
        {},  # {'schema': Base.metadata.schema}  #'autoload': True
    )

    molecule = relation(Molecule, backref=backref("synonyms", order_by=synonym))

    def __init__(self, synonym, note=False):
        """
        synonym is a synonym corresponding to one molecule that is
        identified through fk_molecule_id
        note is an additional note to one synonym, maybe you know
        that a synonym is a PubChem id or is a IUPAC name, than give as note
        CID or IUPAC ... note is not necessary and will default to null
        """
        self.synonym = synonym
        # True is a placeholder for nothing in the parser so ignore any insert
        if note == False:
            self.note = note
        else:
            self.note = None

    def __repr__(self):
        return "Synonyms(%s)" % (self.chembl_id)


class File(Base):
    """
    Store the parsed files for logging purpose.
    """

    __tablename__ = "tbl_files"

    filename = Column(String(200), primary_key=True)
    date = Column(DateTime(timezone=True), default=datetime.datetime.utcnow)

    # __table_args__  = (
    # {'schema': Base.metadata.schema}
    # )

    def __init__(self, filename):
        self.filename = filename


def init(con, schema):
    """
    initialize the database and return the db_engine and the Base-Class for further usage
    any already existing DB won't be overridden, you still get the handle (engine) to the DB
    """
    engine = create_engine(con, pool_recycle=900)
    connection = engine.connect()
    connection.execute("SET search_path TO %s" % schema)
    metadata = MetaData(schema=schema)
    metadata.bind = engine
    # print MetaData
    Base.metadata = metadata
    Base.metadata.create_all(engine)

    return engine


def create_tables(con, schema):
    """
    reset the whole DB
    """
    engine = create_engine(con, pool_recycle=900)
    try:
        # metadata = MetaData(schema = schema)
        # metadata.bind = engine
        # Base.metadata = metadata
        Base.metadata.drop_all(engine)
        Base.metadata.create_all(engine)
    except Exception:
        sys.stderr.write("Can't create table")
        sys.exit(1)


def create_schema(dbname, user, password, host, schema):
    try:
        conn = psycopg2.connect(
            "dbname='%s' user='%s' host='%s' password='%s'"
            % (dbname, user, host, password)
        )
    except Exception:
        sys.stderr.write("Unable to connect to the database.")
        sys.exit(1)

    cur = conn.cursor()
    cur.execute("""DROP SCHEMA IF EXISTS %s CASCADE;""" % schema)
    cur.execute("""CREATE schema %s""" % schema)
    conn.commit()
    conn.close()


def create_indices(con, schema):
    try:
        engine = create_engine(con, pool_recycle=900)
        connection = engine.connect()

        connection.execute(
            """CREATE INDEX idx_mol
                    ON %s.tbl_molecule
                    USING gist(mol);"""
            % schema
        )
        connection.close()
        engine.dispose()
    except Exception:
        sys.stderr.write("Error in creating the GIST index!\n")
        sys.exit(1)


def build_database():
    try:

        host = config.host
        schema = config.schema
        user = config.username
        password = config.password
        database = config.database
        con = "postgresql://%s:%s@%s/%s" % (user, password, host, database)
        create_schema(database, user, password, host, schema)

    except Exception:
        sys.stderr.write("Error: Can't access informations in the config file.")
        sys.exit(1)

    create_tables(con, schema)
    return (schema, con)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="Create, clear and deletes databases.")
    parser.add_argument("--version", action="version", version="%(prog)s 0.1")
    parser.add_argument(
        "--config",
        "-c",
        dest="config_file",
        help="config file with database informations",
    )
    options = parser.parse_args()

    build_database()
