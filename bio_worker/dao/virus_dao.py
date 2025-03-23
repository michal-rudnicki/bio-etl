from typing import Set
from Bio.SeqRecord import SeqRecord
import psycopg2
from dao.taxonomy_dao import TaxonomyDAO
from jobs.ncbi_adapter import NCBIAdapter


class VirusDAO:
    @staticmethod
    def insert_full_virus(conn: psycopg2.extensions.connection, accession: str, record: SeqRecord) -> None:
        cursor = conn.cursor()
        adapter = NCBIAdapter()

        # 🧬 Pobierz taksonomię wirusa
        species = record.annotations.get("organism", "Unknown")
        virus_lineage = adapter.fetch_taxonomy(species)
        taxonomy_id = TaxonomyDAO.insert(conn, virus_lineage) if virus_lineage else None

        # 🧪 Zidentyfikuj typ genomu i gospodarzy
        genome_type = None
        hosts: Set[str] = set()

        for feature in record.features:
            if feature.type == "source":
                if "host" in feature.qualifiers:
                    hosts.update(feature.qualifiers["host"])
                if not genome_type:
                    genome_type = feature.qualifiers.get("mol_type", [None])[0]

        if not genome_type:
            genome_type = record.annotations.get("molecule_type", "unknown")

        # 💾 Wstaw wirusa
        cursor.execute("""
            INSERT INTO viruses (ncbi_id, name, genome_type, species, taxonomy_id)
            VALUES (%s, %s, %s, %s, %s)
            ON CONFLICT (ncbi_id) DO NOTHING
            RETURNING id
        """, (accession, record.name, genome_type, species, taxonomy_id))

        virus_id = cursor.fetchone()
        if not virus_id:
            cursor.execute("SELECT id FROM viruses WHERE ncbi_id = %s", (accession,))
            virus_id = cursor.fetchone()

        virus_id = virus_id[0]

        # 💾 Wstaw sekwencję
        cursor.execute("""
            INSERT INTO sequences (virus_id, type, sequence)
            VALUES (%s, %s, %s)
        """, (
            virus_id,
            "genome",
            str(record.seq)
        ))

        # 🔗 Wstaw gospodarzy
        for host in hosts:
            host_lineage = adapter.fetch_taxonomy(host)
            host_tax_id = TaxonomyDAO.insert(conn, host_lineage) if host_lineage else None

            if host_tax_id:
                cursor.execute("""
                    INSERT INTO virus_hosts (virus_id, host_taxonomy_id)
                    VALUES (%s, %s)
                    ON CONFLICT DO NOTHING
                """, (virus_id, host_tax_id))

        conn.commit()
        cursor.close()