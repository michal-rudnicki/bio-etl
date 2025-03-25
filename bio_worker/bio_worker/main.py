"""Moduł uruchamiający główną pętlę workera
do przetwarzania zadań z bazy danych."""

import time

from bio_worker.db import get_next_job, wait_for_postgres
from bio_worker.jobs.job_processor import process_virus_job


def main() -> None:
    """
    Główna pętla przetwarzająca zadania typu 'fetch_ncbi'.
    """
    conn = wait_for_postgres()

    while True:
        job = get_next_job(conn)
        if job:
            job_id, accession = job
            print(f"[Main] 🧬 Processing job #{job_id}: {accession}")
            process_virus_job(conn, job_id, accession)
        else:
            time.sleep(5)


if __name__ == "__main__":
    main()
