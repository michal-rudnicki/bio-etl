from typing import Union

import psycopg2
from Bio.SeqRecord import SeqRecord


def process_virus_job(
    conn: psycopg2.extensions.connection,
    job_id: int,
    accession: str,
    record: Union[SeqRecord, str],
) -> None:
    cur = conn.cursor()

    if isinstance(record, str):  # error
        cur.execute(
            "UPDATE jobs "
            "SET status = 'error', result = %s, completed_at = now() "
            "WHERE id = %s ",
            (record, job_id),
        )
        conn.commit()
        cur.close()
        return

    # tutaj możesz użyć insert_virus_data()
    # jak wcześniej – lub przepisać go tutaj
    from etl.virus_insert import insert_virus_data

    insert_virus_data(conn, accession, record)

    cur.execute(
        "UPDATE jobs "
        "SET status = 'done', result = 'OK', completed_at = now() "
        "WHERE id = %s",
        (job_id,),
    )
    conn.commit()
    cur.close()
