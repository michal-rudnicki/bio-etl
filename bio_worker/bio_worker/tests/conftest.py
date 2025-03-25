"""Fixture'y testowe: połączenie z bazą danych oraz dane testowe."""
from typing import Generator, List, Tuple
import os
import pytest
import psycopg2
from faker import Faker
from psycopg2.extensions import connection

fake = Faker()


@pytest.fixture(scope="session")
def db_conn() -> Generator[connection, None, None]:
    """
    Inicjalizuje połączenie do bazy PostgreSQL dla testów.

    Yields:
        connection: obiekt połączenia psycopg2.
    """
    conn = psycopg2.connect(
        host=os.environ.get("TESTER_DB_HOST", "localhost"),
        dbname=os.environ.get("TESTER_DB_NAME", "bio"),
        user=os.environ.get("TESTER_DB_USER", "postgres"),
        password=os.environ.get("TESTER_DB_PASS", "password"),
    )
    yield conn
    conn.rollback()  # rollback wszystkiego po sesji
    conn.close()


@pytest.fixture
def random_taxonomy() -> List[Tuple[str, str]]:
    """
    Generuje losową linię taksonomiczną.

    Returns:
        Lista krotek (nazwa, ranga).
    """
    lineage = [
        ("Viruses", "superkingdom"),
        (fake.domain_word(), "family"),
        (fake.domain_word(), "genus"),
        (fake.domain_word(), "species"),
    ]
    return lineage


@pytest.fixture
def fake_accession() -> str:
    """
    Generuje losowy identyfikator GenBank (np. NC_123456.7).

    Returns:
        Łańcuch tekstowy z identyfikatorem.
    """
    return fake.bothify(text="NC_######.##")
