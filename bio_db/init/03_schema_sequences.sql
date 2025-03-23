CREATE TABLE sequences (
    id SERIAL PRIMARY KEY,
    virus_id INTEGER REFERENCES viruses(id),
    type TEXT,
    sequence TEXT,
    length INTEGER GENERATED ALWAYS AS (char_length(sequence)) STORED
);