CREATE TABLE taxonomy (
    id SERIAL PRIMARY KEY,
    name TEXT NOT NULL,
    rank TEXT NOT NULL,
    parent_id INTEGER REFERENCES taxonomy(id),
    UNIQUE (name, rank, parent_id)
);