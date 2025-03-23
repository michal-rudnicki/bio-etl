import psycopg2


class JobDAO:
    @staticmethod
    def mark_done(conn: psycopg2.extensions.connection, job_id: int) -> None:
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
