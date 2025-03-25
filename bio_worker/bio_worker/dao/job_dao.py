"""Moduł DAO do obsługi tabeli jobs – aktualizacja statusów i błędów."""
import psycopg2


class JobDAO:
    """
    DAO do aktualizacji zadań w tabeli jobs.
    """
    @staticmethod
    def mark_done(conn: psycopg2.extensions.connection, job_id: int) -> None:
        """
        Ustawia status zadania na 'done' i zapisuje datę zakończenia.

        Args:
            conn: połączenie z bazą danych.
            job_id: identyfikator zadania.
        """
        cur = conn.cursor()
        cur.execute(
            """
            UPDATE jobs
            SET status = 'done', result = 'OK', completed_at = now()
            WHERE id = %s
        """,
            (job_id,),
        )
        conn.commit()
        cur.close()

    @staticmethod
    def mark_error(
        conn: psycopg2.extensions.connection, job_id: int, error_msg: str
    ) -> None:
        """
        Ustawia status zadania na 'error' i zapisuje treść błędu
        oraz datę zakończenia.

        Args:
            conn: połączenie z bazą danych.
            job_id: identyfikator zadania.
            error_msg: treść błędu do zapisania w kolumnie 'result'.
        """
        cur = conn.cursor()
        cur.execute(
            """
            UPDATE jobs
            SET status = 'error', result = %s, completed_at = now()
            WHERE id = %s
        """,
            (error_msg, job_id),
        )
        conn.commit()
        cur.close()
