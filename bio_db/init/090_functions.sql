-- 090_functions.sql
DO $$ BEGIN RAISE NOTICE 'Tworzę funkcję gc_content(seq)'; END; $$;
CREATE OR REPLACE FUNCTION gc_content(seq TEXT)
RETURNS FLOAT AS $$
    gc = seq.count('G') + seq.count('C')
    return round(100.0 * gc / len(seq), 2) if len(seq) > 0 else 0.0
$$ LANGUAGE plpython3u IMMUTABLE;

DO $$ BEGIN RAISE NOTICE 'Tworzę funkcję generate_uuid()'; END; $$;
CREATE OR REPLACE FUNCTION generate_uuid()
RETURNS uuid AS $$
SELECT uuid_generate_v4();
$$ LANGUAGE SQL;

DO $$ BEGIN RAISE NOTICE 'Tworzę funkcję log_error(job_id, error_code, message)'; END; $$;
CREATE OR REPLACE FUNCTION log_error(
    p_job_id UUID,
    p_error_code TEXT,
    p_message TEXT
) RETURNS void AS $$
BEGIN
    INSERT INTO error_log (job_id, error_code, message)
    VALUES (p_job_id, p_error_code, p_message);
END;
$$ LANGUAGE plpgsql;