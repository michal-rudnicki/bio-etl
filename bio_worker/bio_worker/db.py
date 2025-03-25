"""Moduł do obsługi połączenia z bazą danych PostgreSQL i pobierania zadań."""

import os
import time
from typing import Optional, Tuple

import psycopg2


def wait_for_postgres() -> psycopg2.extensions.connection:
    """
    Oczekuje na dostępność PostgreSQL i zwraca połączenie.

    Raises:
        Exception: jeśli nie uda się połączyć w ciągu 30 prób.
    """
    for i in range(30):
        try:
            conn = psycopg2.connect(
                host=os.environ.get("WORKER_DB_HOST", "bio-db"),
                dbname=os.environ.get("WORKER_DB_NAME", "bio-db"),
                user=os.environ.get("WORKER_DB_USER", "worker"),
                password=os.environ.get("WORKER_DB_PASS", "worker"),
            )
            print("[DB] Connected to PostgreSQL")
            return conn
        except psycopg2.OperationalError:
            print(f"[DB] Waiting for PostgreSQL... ({i+1}/30)")
            time.sleep(2)
    raise RuntimeError("Could not connect to PostgreSQL")


def get_next_job(
    conn: psycopg2.extensions.connection,
) -> Optional[Tuple[int, str]]:
    """
    Pobiera kolejne zadanie typu 'fetch_ncbi' o statusie 'pending'.

    Args:
        conn: Połączenie z bazą danych.

    Returns:
        Krotka (id, parameter) lub None, jeśli brak zadań.
    """
    cur = conn.cursor()
    cur.execute(
        """
        SELECT id, parameter FROM jobs
        WHERE status = 'pending' AND type = 'fetch_ncbi'
        LIMIT 1 FOR UPDATE SKIP LOCKED
    """
    )
    row = cur.fetchone()
    cur.close()
    return row if row else None
