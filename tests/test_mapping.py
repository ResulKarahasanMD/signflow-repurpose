import pytest

from signflow.mapping import GeneRecord, MappingError, strip_ensembl_version, validate_gene_record


def test_strip_ensembl_version():
    assert strip_ensembl_version("ENSG00000240520.3") == "ENSG00000240520"


def test_validate_gene_record_rejects_uox_like_non_protein_coding():
    with pytest.raises(MappingError):
        validate_gene_record(
            GeneRecord(
                symbol="UOX",
                ensembl_id="ENSG00000240520.3",
                biotype="transcribed_unitary_pseudogene",
            )
        )
