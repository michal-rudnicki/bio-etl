"""Testy jednostkowe DAO do zarządzania zadaniami (jobs)."""
from psycopg2.extensions import connection
from bio_worker.dao.job_dao import JobDAO


def test_job_lifecycle(db_conn: connection) -> None:
    """
    Testuje pełen cykl życia zadania: insert → mark_done → sprawdzenie statusu.
    """
    cur = db_conn.cursor()
    cur.execute(
        """
        INSERT INTO jobs (type, parameter, created_by)
        VALUES ('fetch_ncbi', 'NC_123456.1', 'tester')
        RETURNING id
    """
    )
    job_id = cur.fetchone()[0]
    db_conn.commit()

    JobDAO.mark_done(db_conn, job_id)

    cur.execute("SELECT status FROM jobs WHERE id = %s", (job_id,))
    status = cur.fetchone()[0]
    assert status == "done"
