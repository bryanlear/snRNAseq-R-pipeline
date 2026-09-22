from pathlib import Path
from geneformer import TranscriptomeTokenizer

project = Path("/Users/bry_lee/snrnaseq_r_pipeline")
input_dir = project / "data_preparation_geneformer/h5ad_RNA"
output_dir = project / "geneformer_output/tokenized"

output_dir.mkdir(parents=True, exist_ok=True)

tokenizer = TranscriptomeTokenizer(custom_attr_name_dict={"sample_id": "sample_id", "condition": "condition", "batch": "batch", "cell.type": "cell_type",}, nproc=4, model_version="V2",)
tokenizer.tokenize_data(data_directory=str(input_dir), output_directory=str(output_dir), output_prefix="gse174367_rna_v2", file_format="h5ad",)
