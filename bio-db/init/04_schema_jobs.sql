CREATE TABLE jobs (
    id SERIAL PRIMARY KEY,
    type TEXT,
    parameter TEXT,
    status TEXT DEFAULT 'pending',
    result TEXT,
    created_by TEXT,
    created_at TIMESTAMP DEFAULT now(),
    completed_at TIMESTAMP
);

CREATE FUNCTION set_created_by()
RETURNS TRIGGER AS $$
BEGIN
  IF NEW.created_by IS NULL THEN
    NEW.created_by := current_user;
  END IF;
  RETURN NEW;
END;
$$ LANGUAGE plpgsql;

CREATE TRIGGER trg_jobs_created_by
BEFORE INSERT ON jobs
FOR EACH ROW EXECUTE FUNCTION set_created_by();