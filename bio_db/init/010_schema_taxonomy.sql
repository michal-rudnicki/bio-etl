-- 010_schema_taxonomy.sql
DO $$ BEGIN RAISE NOTICE 'Tworzę tabelę taxonomy'; END; $$;
CREATE TABLE IF NOT EXISTS taxonomy (
    tax_id INTEGER PRIMARY KEY,
    parent_id INTEGER,
    rank TEXT,
    name_txt TEXT,
    unique_name TEXT,
    name_class TEXT,
    job_id UUID REFERENCES jobs(job_id) ON DELETE SET NULL,
    UNIQUE(name_txt, name_class)
);

DO $$ BEGIN RAISE NOTICE 'Tworzę indeks idx_taxonomy_parent_id'; END; $$;
CREATE INDEX IF NOT EXISTS idx_taxonomy_parent_id ON taxonomy(parent_id);