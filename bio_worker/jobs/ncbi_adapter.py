import os
from typing import Union, List, Tuple
from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord

Entrez.email = os.environ.get("EMAIL", "default@example.com")


class NCBIAdapter:
    def fetch_genbank(self, accession: str) -> Union[SeqRecord, str]:
        try:
            handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
            record = SeqIO.read(handle, "genbank")
            handle.close()
            return record
        except Exception as e:
            return str(e)

    def fetch_taxonomy(self, name: str) -> List[Tuple[str, str]]:
        try:
            handle = Entrez.esearch(db="taxonomy", term=name + "[Organism]", retmode="xml")
            record = Entrez.read(handle)
            handle.close()

            if not record["IdList"]:
                print(f"[NCBI] ❗ No taxonomy for: {name}")
                return []

            tax_id = record["IdList"][0]
            handle = Entrez.efetch(db="taxonomy", id=tax_id, retmode="xml")
            data = Entrez.read(handle)
            handle.close()

            lineage = data[0].get("LineageEx", [])
            lineage_data = [(item["ScientificName"], item["Rank"]) for item in lineage]
            lineage_data.append((data[0]["ScientificName"], data[0]["Rank"]))

            return lineage_data
        except Exception as e:
            print(f"[NCBI] ❌ Taxonomy fetch error: {e}")
            return []