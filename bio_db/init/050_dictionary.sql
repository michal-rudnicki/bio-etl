-- 050_schema_dictionary.sql
DO $$ BEGIN RAISE NOTICE 'Tworzę tabelę dictionary'; END; $$;
CREATE TABLE IF NOT EXISTS dictionary (
    id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    category TEXT NOT NULL,            -- np. 'error', 'warning', 'info', 'job_type', 'sequence_type'
    code TEXT NOT NULL,                -- np. 'SEQ_TOO_SHORT', 'VIRUS_UNKNOWN'
    message TEXT NOT NULL,             -- komunikat dla użytkownika
    description TEXT,                  -- opcjonalny opis techniczny
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    UNIQUE(category, code)
);

DO $$ BEGIN RAISE NOTICE 'Tworzę indeksy dla tabeli dictionary'; END; $$;
CREATE INDEX IF NOT EXISTS idx_dictionary_code ON dictionary(code);
CREATE INDEX IF NOT EXISTS idx_dictionary_category ON dictionary(category);