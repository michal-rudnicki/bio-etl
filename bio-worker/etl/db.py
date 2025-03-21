import os
import time
from typing import Optional, Tuple

import psycopg2


def wait_for_postgres() -> psycopg2.extensions.connection:
    for i in range(30):
        try:
            conn = psycopg2.connect(
                host=os.environ["DB_HOST"],
                dbname=os.environ["DB_NAME"],
                user=os.environ["DB_USER"],
                password=os.environ["DB_PASS"],
            )
            print("[Worker] Connected to PostgreSQL")
            return conn
        except psycopg2.OperationalError:
            print(f"[Worker] Waiting for database... ({i+1}/30)")
            time.sleep(2)
    raise Exception(" PostgreSQL not available")


def get_next_job(
    conn: psycopg2.extensions.connection,
) -> Optional[Tuple[int, str]]:
    cur = conn.cursor()
    cur.execute(
        "SELECT id, parameter FROM jobs "
        "WHERE status = 'pending' AND type = 'fetch_ncbi' "
        "LIMIT 1 FOR UPDATE SKIP LOCKED"
    )
    row = cur.fetchone()
    cur.close()
    return row if row else None
