import os
from typing import Union

from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord

Entrez.email = os.environ.get("EMAIL", "default@example.com")


class GenBankAdapter:
    def fetch(self, accession: str) -> Union[SeqRecord, str]:
        try:
            handle = Entrez.efetch(
                db="nucleotide", id=accession, rettype="gb", retmode="text"
            )
            record = SeqIO.read(handle, "genbank")
            handle.close()
            return record
        except Exception as e:
            return str(e)
