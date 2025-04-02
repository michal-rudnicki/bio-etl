-- 005_schema_jobs.sql
DO $$ BEGIN RAISE NOTICE 'Tworzę tabelę jobs'; END; $$;
CREATE TABLE IF NOT EXISTS jobs (
    job_id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    job_type TEXT,
    status TEXT,
    input_params JSONB,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    finished_at TIMESTAMP,
    status_message TEXT
);
DO $$ BEGIN RAISE NOTICE 'Tworzę indeks idx_jobs_created_at'; END; $$;
CREATE INDEX IF NOT EXISTS idx_jobs_created_at ON jobs(created_at);

DO $$ BEGIN RAISE NOTICE 'Tworzę funkcje set_created_by'; END; $$;
CREATE FUNCTION set_created_by()
RETURNS TRIGGER AS $$
BEGIN
  IF NEW.created_by IS NULL THEN
    NEW.created_by := current_user;
  END IF;
  RETURN NEW;
END;
$$ LANGUAGE plpgsql;

DO $$ BEGIN RAISE NOTICE 'Tworzę trigger trg_job_created_by'; END; $$;
CREATE TRIGGER trg_jobs_created_by
BEFORE INSERT ON jobs
FOR EACH ROW EXECUTE FUNCTION set_created_by();