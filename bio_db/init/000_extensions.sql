-- 00_extensions.sql
DO $$ BEGIN RAISE NOTICE 'Instaluję rozszerzenia: plpython3u'; END; $$;
CREATE EXTENSION IF NOT EXISTS "plpython3u";

DO $$ BEGIN RAISE NOTICE 'Instaluję rozszerzenia: uuid-ossp'; END; $$;
CREATE EXTENSION IF NOT EXISTS "uuid-ossp";