import time
import psycopg2
from Bio import Entrez, SeqIO
import os

Entrez.email = os.environ.get("EMAIL", "default@example.com")

def wait_for_postgres():
    for i in range(30):
        try:
            conn = psycopg2.connect(
                host=os.environ['DB_HOST'],
                dbname=os.environ['DB_NAME'],
                user=os.environ['DB_USER'],
                password=os.environ['DB_PASS']
            )
            print("[Worker] Connected to PostgreSQL")
            return conn
        except psycopg2.OperationalError:
            print(f"[Worker] Waiting for database... ({i+1}/30)")
            time.sleep(2)
    raise Exception(" PostgreSQL not available")

def fetch_genbank(accession):
    try:
        print(f"[Worker] Fetching {accession} from NCBI")
        handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        handle.close()
        return record
    except Exception as e:
        return str(e)

def fetch_taxonomy_lineage(name):
    try:
        handle = Entrez.esearch(db="taxonomy", term=name + "[Organism]", retmode="xml")
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
        lineage_data = [(item["ScientificName"], item["Rank"]) for item in lineage]
        lineage_data.append((data[0]["ScientificName"], data[0]["Rank"]))
        return lineage_data
    except Exception as e:
        print(f"[Worker] Error fetching taxonomy for {name}: {e}")
        return []

def insert_taxonomy(conn, lineage_data):
    cursor = conn.cursor()
    parent_id = None
    for name, rank in lineage_data:
        cursor.execute("""
            SELECT id FROM taxonomy
            WHERE name = %s AND rank = %s AND parent_id IS NOT DISTINCT FROM %s
        """, (name, rank, parent_id))
        row = cursor.fetchone()
        if row:
            tax_id = row[0]
        else:
            cursor.execute("""
                INSERT INTO taxonomy (name, rank, parent_id)
                VALUES (%s, %s, %s)
                RETURNING id
            """, (name, rank, parent_id))
            tax_id = cursor.fetchone()[0]
        parent_id = tax_id
    conn.commit()
    cursor.close()
    return parent_id

def insert_virus_data(conn, accession, record):
    cursor = conn.cursor()

    # Wirus – taksonomia
    organism = record.annotations.get('organism', 'Unknown')
    virus_lineage = fetch_taxonomy_lineage(organism)
    taxonomy_id = insert_taxonomy(conn, virus_lineage) if virus_lineage else None

    # Gospodarze – mogą być wielokrotni
    hosts = set()
    genome_type = None
    for feature in record.features:
        if feature.type == 'source':
            if 'host' in feature.qualifiers:
                hosts.update(feature.qualifiers['host'])
            if not genome_type:
                genome_type = feature.qualifiers.get('mol_type', [None])[0]

    if not genome_type:
        genome_type = record.annotations.get('molecule_type', 'unknown')

    # Wstaw wirusa
    cursor.execute("""
        INSERT INTO viruses (ncbi_id, name, genome_type, species, taxonomy_id)
        VALUES (%s, %s, %s, %s, %s)
        ON CONFLICT (ncbi_id) DO NOTHING
        RETURNING id
    """, (
        accession,
        record.name,
        genome_type,
        organism,
        taxonomy_id
    ))

    virus_id = cursor.fetchone()
    if not virus_id:
        cursor.execute("SELECT id FROM viruses WHERE ncbi_id = %s", (accession,))
        virus_id = cursor.fetchone()
    virus_id = virus_id[0]

    # Wstaw sekwencję
    cursor.execute("""
        INSERT INTO sequences (virus_id, type, sequence)
        VALUES (%s, %s, %s)
    """, (
        virus_id,
        "genome",
        str(record.seq)
    ))

    # Wstaw relacje gospodarzy
    for host in hosts:
        host_lineage = fetch_taxonomy_lineage(host)
        if host_lineage:
            host_tax_id = insert_taxonomy(conn, host_lineage)
            cursor.execute("""
                INSERT INTO virus_hosts (virus_id, host_taxonomy_id)
                VALUES (%s, %s)
                ON CONFLICT DO NOTHING
            """, (virus_id, host_tax_id))

    conn.commit()
    cursor.close()

def main():
    conn = wait_for_postgres()
    while True:
        cursor = conn.cursor()
        cursor.execute("""
            SELECT id, parameter FROM jobs
            WHERE status = 'pending' AND type = 'fetch_ncbi'
            LIMIT 1 FOR UPDATE SKIP LOCKED
        """)
        row = cursor.fetchone()
        if row:
            job_id, accession = row
            print(f"[Worker] Processing job #{job_id}: {accession}")
            result = fetch_genbank(accession)
            if isinstance(result, str):
                cursor.execute("""
                    UPDATE jobs
                    SET status = 'error', result = %s, completed_at = now()
                    WHERE id = %s
                """, (result, job_id))
            else:
                insert_virus_data(conn, accession, result)
                cursor.execute("""
                    UPDATE jobs
                    SET status = 'done', result = 'OK', completed_at = now()
                    WHERE id = %s
                """, (job_id,))
            conn.commit()
        else:
            time.sleep(5)

main()