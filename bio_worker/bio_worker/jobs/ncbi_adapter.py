"""Adapter do komunikacji z serwisem NCBI (Entrez, GenBank, Taxonomy)."""
import os
from typing import Union, List, Tuple
from Bio import Entrez, SeqIO
from Bio.Entrez.Parser import ValidationError
from Bio.SeqRecord import SeqRecord


class NCBIAdapter:
    """
    Adapter do pobierania danych z NCBI: sekwencje GenBank i taksonomia.
    """

    def __init__(self) -> None:
        """
        Inicjalizuje adapter, ustawia e-mail wymagany przez NCBI.
        """
        Entrez.email = os.environ.get("EMAIL", "default@example.com")

    def fetch_genbank(self, accession: str) -> Union[SeqRecord, str]:
        """
        Pobiera rekord GenBank o podanym numerze dostępu.

        Args:
            accession: identyfikator GenBank (np. NC_045512.2)

        Returns:
            SeqRecord lub komunikat błędu.
        """
        try:
            handle = Entrez.efetch(
                db="nucleotide", id=accession, rettype="gb", retmode="text"
            )
            record = SeqIO.read(handle, "genbank")
            handle.close()
            return record
        except (ValidationError, ValueError, TypeError) as error:
            return str(error)

    def fetch_taxonomy(self, name: str) -> List[Tuple[str, str]]:
        """
        Pobiera linię taksonomiczną organizmu na podstawie jego nazwy.

        Args:
            name: pełna nazwa organizmu.

        Returns:
            Lista krotek (nazwa, ranga) reprezentujących drzewo taksonomiczne.
        """
        try:
            handle = Entrez.esearch(
                db="taxonomy", term=name + "[Organism]", retmode="xml"
            )
            record = Entrez.read(handle)
            handle.close()

            if not record["IdList"]:
                print(f"[NCBI] No taxonomy for: {name}")
                return []

            tax_id = record["IdList"][0]
            handle = Entrez.efetch(db="taxonomy", id=tax_id, retmode="xml")
            data = Entrez.read(handle)
            handle.close()

            lineage = data[0].get("LineageEx", [])
            lineage_data = [
                (item["ScientificName"], item["Rank"]) for item in lineage
            ]
            lineage_data.append((data[0]["ScientificName"], data[0]["Rank"]))

            return lineage_data
        except (ValidationError, ValueError, TypeError) as error:
            print(f"[NCBI] Taxonomy fetch error: {error}")
            return []
