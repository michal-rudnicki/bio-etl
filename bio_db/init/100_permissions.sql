DO $$ BEGIN RAISE NOTICE 'Tworzę użytkownika worker'; END; $$;
CREATE USER worker WITH PASSWORD 'worker';
GRANT CONNECT ON DATABASE "bio-db" TO worker;
GRANT USAGE ON SCHEMA public TO worker;
GRANT SELECT, INSERT, UPDATE ON ALL TABLES IN SCHEMA public TO worker;
GRANT USAGE, SELECT ON ALL SEQUENCES IN SCHEMA public TO worker;

GRANT EXECUTE ON FUNCTION gc_content(TEXT) TO worker;
GRANT EXECUTE ON FUNCTION set_created_by() TO worker;
GRANT EXECUTE ON FUNCTION generate_uuid() TO worker;
GRANT EXECUTE ON FUNCTION log_error(UUID, TEXT, TEXT) TO worker;

DO $$ BEGIN RAISE NOTICE 'Tworzę użytkownika tester'; END; $$;
CREATE USER tester WITH PASSWORD 'tester';
GRANT CONNECT ON DATABASE "bio-db" TO tester;
GRANT USAGE ON SCHEMA public TO tester;
GRANT SELECT, INSERT, UPDATE, DELETE ON ALL TABLES IN SCHEMA public TO tester;
GRANT USAGE, SELECT ON ALL SEQUENCES IN SCHEMA public TO tester;

GRANT EXECUTE ON FUNCTION gc_content(TEXT) TO tester;
GRANT EXECUTE ON FUNCTION set_created_by() TO tester;
GRANT EXECUTE ON FUNCTION generate_uuid() TO tester;
GRANT EXECUTE ON FUNCTION log_error(UUID, TEXT, TEXT) TO tester;