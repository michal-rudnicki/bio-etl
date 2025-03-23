import psycopg2
from Bio.SeqRecord import SeqRecord
from dao.taxonomy_dao import TaxonomyDAO


class VirusDAO:
    @staticmethod
    def insert_full_virus(
        conn: psycopg2.extensions.connection, accession: str, record: SeqRecord
    ) -> None:
        cur = conn.cursor()

        species = record.annotations.get("organism", "Unknown")
        taxonomy = TaxonomyDAO.fetch_lineage(species)
        taxonomy_id = TaxonomyDAO.insert(conn, taxonomy) if taxonomy else None

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
        cur.execute(
            """
            INSERT INTO viruses
            (ncbi_id, name, genome_type, species, taxonomy_id)
            VALUES (%s, %s, %s, %s, %s)
            ON CONFLICT (ncbi_id) DO NOTHING
            RETURNING id
        """,
            (accession, record.name, genome_type, species, taxonomy_id),
        )

        virus_id = cur.fetchone()
        if not virus_id:
            cur.execute(
                "SELECT id FROM viruses WHERE ncbi_id = %s", (accession,)
            )
            virus_id = cur.fetchone()
        virus_id = virus_id[0]

        # Sekwencja
        cur.execute(
            """
            INSERT INTO sequences (virus_id, type, sequence)
            VALUES (%s, %s, %s)
        """,
            (virus_id, "genome", str(record.seq)),
        )

        # Gospodarze
        for host in hosts:
            host_lineage = TaxonomyDAO.fetch_lineage(host)
            host_tax_id = (
                TaxonomyDAO.insert(conn, host_lineage)
                if host_lineage
                else None
            )
            if host_tax_id:
                cur.execute(
                    """
                    INSERT INTO virus_hosts (virus_id, host_taxonomy_id)
                    VALUES (%s, %s)
                    ON CONFLICT DO NOTHING
                """,
                    (virus_id, host_tax_id),
                )

        conn.commit()
        cur.close()
