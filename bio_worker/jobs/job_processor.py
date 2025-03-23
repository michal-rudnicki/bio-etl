import psycopg2
from dao.job_dao import JobDAO
from dao.virus_dao import VirusDAO
from jobs.ncbi_adapter import NCBIAdapter

def process_virus_job(
    conn: psycopg2.extensions.connection, job_id: int, accession: str
) -> None:
    adapter = NCBIAdapter()
    record = adapter.fetch_genbank(accession)

    if isinstance(record, str):  # Błąd
        JobDAO.mark_error(conn, job_id, record)
        return

    VirusDAO.insert_full_virus(conn, accession, record)
    JobDAO.mark_done(conn, job_id)
