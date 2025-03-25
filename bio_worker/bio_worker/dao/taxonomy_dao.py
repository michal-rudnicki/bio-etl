"""Moduł DAO do zapisu taksonomii do tabeli taxonomy w bazie danych."""
from typing import List, Tuple, Optional
import psycopg2


class TaxonomyDAO:
    """
    DAO do obsługi zapisu drzewa taksonomii organizmów
    (np. wirusów, gospodarzy).
    """
    @staticmethod
    def insert(
        conn: psycopg2.extensions.connection, lineage: List[Tuple[str, str]]
    ) -> Optional[int]:
        """
        Wstawia kolejne poziomy drzewa taksonomii do tabeli taxonomy.

        Jeśli dane już istnieją (ten sam parent + name + rank),
        zostają pominięte.

        Args:
            conn: połączenie do bazy PostgreSQL.
            lineage: lista krotek (nazwa, ranga) reprezentujących
            linię taksonomiczną.

        Returns:
            ID ostatniego wstawionego lub znalezionego węzła (czyli gatunku).
        """
        cur = conn.cursor()
        parent_id = None

        for name, rank in lineage:
            cur.execute(
                """
                SELECT id FROM taxonomy
                WHERE name = %s
                AND rank = %s
                AND parent_id IS NOT DISTINCT FROM %s
            """,
                (name, rank, parent_id),
            )
            row = cur.fetchone()

            if row:
                tax_id = row[0]
            else:
                cur.execute(
                    """
                    INSERT INTO taxonomy (name, rank, parent_id)
                    VALUES (%s, %s, %s)
                    RETURNING id
                """,
                    (name, rank, parent_id),
                )
                tax_id = cur.fetchone()[0]

            parent_id = tax_id

        conn.commit()
        cur.close()
        return parent_id
