CREATE TABLE viruses (
    id SERIAL PRIMARY KEY,
    ncbi_id TEXT UNIQUE,
    name TEXT,
    genome_type TEXT,
    species TEXT,
    taxonomy_id INTEGER REFERENCES taxonomy(id),
    created_at TIMESTAMP DEFAULT now()
);

CREATE TABLE virus_hosts (
    virus_id INTEGER REFERENCES viruses(id) ON DELETE CASCADE,
    host_taxonomy_id INTEGER REFERENCES taxonomy(id) ON DELETE CASCADE,
    PRIMARY KEY (virus_id, host_taxonomy_id)
);