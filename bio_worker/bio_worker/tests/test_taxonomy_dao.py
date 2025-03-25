"""Testy jednostkowe DAO do zapisu taksonomii."""
from typing import List, Tuple
from psycopg2.extensions import connection
from bio_worker.dao.taxonomy_dao import TaxonomyDAO


def test_insert_taxonomy(
    db_conn: connection, random_taxonomy: List[Tuple[str, str]]
) -> None:
    """
    Testuje, czy insert_taxonomy wstawia poprawnie dane do tabeli taxonomy.

    Sprawdza, czy zwracany ID to liczba całkowita.
    """
    taxonomy_id = TaxonomyDAO.insert(db_conn, random_taxonomy)
    assert isinstance(taxonomy_id, int)
