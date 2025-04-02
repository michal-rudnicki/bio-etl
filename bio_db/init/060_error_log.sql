-- 060_schema_error_log.sql
DO $$ BEGIN RAISE NOTICE 'Tworzę tabelę error_log'; END; $$;
CREATE TABLE IF NOT EXISTS error_log (
    id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    job_id UUID REFERENCES jobs(job_id) ON DELETE SET NULL,
    error_code TEXT,
    message TEXT,
    logged_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);
DO $$ BEGIN RAISE NOTICE 'Tworzę indeksy dla tabeli error_log'; END; $$;
CREATE INDEX IF NOT EXISTS idx_error_log_job_id ON error_log(job_id);
CREATE INDEX IF NOT EXISTS idx_error_log_logged_at ON error_log(logged_at);