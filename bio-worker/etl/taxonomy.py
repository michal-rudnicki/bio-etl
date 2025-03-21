from typing import List, Optional, Tuple

import psycopg2
from Bio import Entrez


def fetch_taxonomy_lineage(name: str) -> List[Tuple[str, str]]:
    try:
        handle = Entrez.esearch(
            db="taxonomy", term=name + "[Organism]", retmode="xml"
        )
        record = Entrez.read(handle)
        handle.close()
        if not record["IdList"]:
            print(f"[Worker]️ No taxonomy found for: {name}")
            return []
        tax_id = record["IdList"][0]
        handle = Entrez.efetch(db="taxonomy", id=tax_id, retmode="xml")
        data = Entrez.read(handle)
        handle.close()
        lineage = data[0].get("LineageEx", [])
        lineage_data = [
            (item["ScientificName"], item["Rank"]) for item in lineage
        ]
        lineage_data.append((data[0]["ScientificName"], data[0]["Rank"]))
        return lineage_data
    except Exception as e:
        print(f"[Worker] Error fetching taxonomy for {name}: {e}")
        return []


def insert_taxonomy(
    conn: psycopg2.extensions.connection, lineage_data: List[Tuple[str, str]]
) -> Optional[int]:
    cursor = conn.cursor()
    parent_id = None
    for name, rank in lineage_data:
        cursor.execute(
            "SELECT id FROM taxonomy "
            "WHERE name = %s AND rank = %s "
            "AND parent_id IS NOT DISTINCT FROM %s ",
            (name, rank, parent_id),
        )
        row = cursor.fetchone()
        if row:
            tax_id = row[0]
        else:
            cursor.execute(
                "INSERT INTO taxonomy (name, rank, parent_id) "
                "VALUES (%s, %s, %s) "
                "RETURNING id ",
                (name, rank, parent_id),
            )
            tax_id = cursor.fetchone()[0]
        parent_id = tax_id
    conn.commit()
    cursor.close()
    return parent_id
