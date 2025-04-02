-- 030_schema_sequences.sql
DO $$ BEGIN RAISE NOTICE 'Tworzę tabelę sequences'; END; $$;
CREATE TABLE IF NOT EXISTS sequences (
    id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    accession TEXT UNIQUE,
    sequence TEXT,
    sequence_type TEXT,
    organism INTEGER REFERENCES taxonomy(tax_id),
    length INTEGER,
    source TEXT,
    date DATE,
    job_id UUID REFERENCES jobs(job_id) ON DELETE SET NULL
);

DO $$ BEGIN RAISE NOTICE 'Tworzę indeksy dla tabeli sequences'; END; $$;
CREATE INDEX IF NOT EXISTS idx_sequences_organism ON sequences(organism);
CREATE INDEX IF NOT EXISTS idx_sequences_accession ON sequences(accession);
CREATE INDEX IF NOT EXISTS idx_sequences_date ON sequences(date);