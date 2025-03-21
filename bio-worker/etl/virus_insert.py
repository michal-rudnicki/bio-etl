import psycopg2
from Bio.SeqRecord import SeqRecord
from etl.taxonomy import fetch_taxonomy_lineage, insert_taxonomy


def insert_virus_data(
    conn: psycopg2.extensions.connection, accession: str, record: SeqRecord
) -> None:
    cursor = conn.cursor()

    # Taksonomia wirusa
    organism = record.annotations.get("organism", "Unknown")
    virus_lineage = fetch_taxonomy_lineage(organism)
    taxonomy_id = (
        insert_taxonomy(conn, virus_lineage) if virus_lineage else None
    )

    # Typ genomu
    genome_type = None
    hosts = set()

    for feature in record.features:
        if feature.type == "source":
            if "host" in feature.qualifiers:
                hosts.update(feature.qualifiers["host"])
            if not genome_type:
                genome_type = feature.qualifiers.get("mol_type", [None])[0]

    if not genome_type:
        genome_type = record.annotations.get("molecule_type", "unknown")

    # Wstaw wirusa
    cursor.execute(
        "INSERT INTO viruses "
        "(ncbi_id, name, genome_type, species, taxonomy_id) "
        "VALUES (%s, %s, %s, %s, %s) "
        "ON CONFLICT (ncbi_id) DO NOTHING "
        "RETURNING id ",
        (accession, record.name, genome_type, organism, taxonomy_id),
    )

    virus_id = cursor.fetchone()
    if not virus_id:
        cursor.execute(
            "SELECT id FROM viruses WHERE ncbi_id = %s ", (accession,)
        )
        virus_id = cursor.fetchone()
    virus_id = virus_id[0]

    # 📥 Wstaw sekwencję
    cursor.execute(
        "INSERT INTO sequences (virus_id, type, sequence) "
        "VALUES (%s, %s, %s)",
        (virus_id, "genome", str(record.seq)),
    )

    # Relacje gospodarzowe
    for host in hosts:
        host_lineage = fetch_taxonomy_lineage(host)
        if host_lineage:
            host_tax_id = insert_taxonomy(conn, host_lineage)
            cursor.execute(
                "INSERT INTO virus_hosts (virus_id, host_taxonomy_id) "
                "VALUES (%s, %s) "
                "ON CONFLICT DO NOTHING ",
                (virus_id, host_tax_id),
            )

    conn.commit()
    cursor.close()
