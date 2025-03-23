from typing import List, Tuple, Optional
import psycopg2


class TaxonomyDAO:
    @staticmethod
    def insert(conn: psycopg2.extensions.connection, lineage: List[Tuple[str, str]]) -> Optional[int]:
        cur = conn.cursor()
        parent_id = None

        for name, rank in lineage:
            cur.execute("""
                SELECT id FROM taxonomy
                WHERE name = %s AND rank = %s AND parent_id IS NOT DISTINCT FROM %s
            """, (name, rank, parent_id))
            row = cur.fetchone()

            if row:
                tax_id = row[0]
            else:
                cur.execute("""
                    INSERT INTO taxonomy (name, rank, parent_id)
                    VALUES (%s, %s, %s)
                    RETURNING id
                """, (name, rank, parent_id))
                tax_id = cur.fetchone()[0]

            parent_id = tax_id

        conn.commit()
        cur.close()
        return parent_id