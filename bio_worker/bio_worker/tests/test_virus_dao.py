"""Testy DAO do wstawiania wirusów, sekwencji i gospodarzy."""
from io import StringIO
from Bio import SeqIO
from psycopg2.extensions import connection
from bio_worker.dao.virus_dao import VirusDAO


FAKE_GENBANK = """
LOCUS       TEST00001             25 bp    DNA     linear   SYN 01-JAN-1980
DEFINITION  test virus.
ACCESSION   TEST00001
ORIGIN
        1 atgctgacgt agcatgcagt gctag
//
"""


def test_insert_full_virus(db_conn: connection, fake_accession: str) -> None:
    """
    Testuje wstawienie wirusa do bazy na podstawie minimalnego GenBank.

    Sprawdza brak błędów przy insertach.
    """
    record = SeqIO.read(StringIO(FAKE_GENBANK), "genbank")
    record.annotations["organism"] = "Testvirus infectus"
    record.name = fake_accession
    record.features = []

    VirusDAO.insert_full_virus(db_conn, fake_accession, record)
