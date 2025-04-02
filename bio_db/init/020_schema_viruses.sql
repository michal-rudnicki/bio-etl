-- 020_schema_viruses.sql
DO $$ BEGIN RAISE NOTICE 'Tworzę tabelę viruses'; END; $$;
CREATE TABLE IF NOT EXISTS viruses (
    accession TEXT PRIMARY KEY,
    name TEXT,
    taxonomy_id INTEGER REFERENCES taxonomy(tax_id),
    host TEXT,
    segment TEXT,
    strain TEXT,
    genotype TEXT,
    serotype TEXT,
    job_id UUID REFERENCES jobs(job_id) ON DELETE SET NULL
);

DO $$ BEGIN RAISE NOTICE 'Tworzę indeks idx_viruses_taxonomy_id'; END; $$;
CREATE INDEX IF NOT EXISTS idx_viruses_taxonomy_id ON viruses(taxonomy_id);