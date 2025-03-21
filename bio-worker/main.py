import time

from etl.db import get_next_job, wait_for_postgres
from etl.genbank import fetch_genbank
from etl.virus_loader import process_virus_job


def main() -> None:
    conn = wait_for_postgres()
    while True:
        job = get_next_job(conn)
        if job:
            job_id, accession = job
            print(f"Processing job {accession}")
            record = fetch_genbank(accession)
            process_virus_job(conn, job_id, accession, record)
        else:
            time.sleep(5)


if __name__ == "__main__":
    main()
