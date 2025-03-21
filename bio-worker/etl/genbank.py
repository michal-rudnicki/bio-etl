import os
from typing import Union

from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord


def fetch_genbank(accession: str) -> Union[SeqRecord, str]:
    try:
        Entrez.email = os.environ.get("EMAIL", "default@example.com")
        print(f"[Worker] Fetching {accession} from NCBI")
        handle = Entrez.efetch(
            db="nucleotide", id=accession, rettype="gb", retmode="text"
        )
        record = SeqIO.read(handle, "genbank")
        handle.close()
        return record
    except Exception as e:
        return str(e)
