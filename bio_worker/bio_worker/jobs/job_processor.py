"""Moduł odpowiedzialny za przetwarzanie zadań typu 'fetch_ncbi'."""
import psycopg2
from bio_worker.jobs.ncbi_adapter import NCBIAdapter
from bio_worker.dao.virus_dao import VirusDAO
from bio_worker.dao.job_dao import JobDAO


def process_virus_job(
    conn: psycopg2.extensions.connection, job_id: int, accession: str
) -> None:
    """
    Przetwarza zadanie 'fetch_ncbi': pobiera dane z NCBI,
     zapisuje wirusa do bazy,
    aktualizuje status zadania.

    Args:
        conn: połączenie do bazy danych.
        job_id: identyfikator zadania.
        accession: identyfikator sekwencji NCBI (np. NC_045512.2).
    """
    adapter = NCBIAdapter()
    record = adapter.fetch_genbank(accession)

    if isinstance(record, str):  # Błąd
        JobDAO.mark_error(conn, job_id, record)
        return

    VirusDAO.insert_full_virus(conn, accession, record)
    JobDAO.mark_done(conn, job_id)
