CREATE OR REPLACE FUNCTION gc_content(seq TEXT)
RETURNS FLOAT AS $$
    gc = seq.count('G') + seq.count('C')
    return round(100.0 * gc / len(seq), 2) if len(seq) > 0 else 0.0
$$ LANGUAGE plpython3u IMMUTABLE;