from __future__ import annotations
import argparse
import json
import torch
from pathlib import Path
from datasets import load_from_disk
from geneformer import EmbExtractor

PROJECT_DIR = Path(__file__).resolve().parent
MODEL_DIR = PROJECT_DIR / "geneformer_v2"
GENEFORMER_OUTPUT_DIR = PROJECT_DIR / "geneformer_output"
TOKENIZED_DIR = GENEFORMER_OUTPUT_DIR / "tokenized"
DATASET_DIR = TOKENIZED_DIR / "gse174367_rna_v2.dataset"
EMBEDDING_OUTPUT_DIR = GENEFORMER_OUTPUT_DIR / "embeddings"
REQUIRED_MODEL_FILES = ("config.json", "model.safetensors")
EMBEDDING_LABELS = ("sample_id", "condition", "batch", "cell_type")

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Extract mean-pooled Geneformer V2 cell embeddings.")
    parser.add_argument("--model-dir", type=Path, default=MODEL_DIR, help=f"Geneformer model directory. Default: {MODEL_DIR}")
    parser.add_argument("--dataset-dir", type=Path, default=DATASET_DIR, help=f"Tokenized Hugging Face dataset. Default: {DATASET_DIR}")
    parser.add_argument("--output-dir", type=Path, default=EMBEDDING_OUTPUT_DIR, help=f"Embedding output directory. Default: {EMBEDDING_OUTPUT_DIR}")
    parser.add_argument("--output-prefix", default="gse174367_rna_v2_cell_embeddings", help="Prefix for the CSV and PyTorch output files.")
    parser.add_argument("--batch-size", type=int, default=1, help="GPU forward-pass batch size. Start with 1 for 4096-token cells.")
    parser.add_argument("--nproc", type=int, default=4, help="Number of CPU processes used for dataset operations.")
    parser.add_argument("--max-cells", type=int, default=None, help="Optional cell limit for a test run. The default processes all cells.")
    return parser.parse_args()

def validate_inputs(args: argparse.Namespace) -> tuple[int, int]:
    missing_model_files = [args.model_dir / name for name in REQUIRED_MODEL_FILES if not (args.model_dir / name).is_file()]
    if missing_model_files:
        missing = "\n".join(str(path) for path in missing_model_files)
        raise FileNotFoundError(f"Missing model files:\n{missing}")
    if not args.dataset_dir.is_dir():
        raise FileNotFoundError(f"Tokenized dataset not found: {args.dataset_dir}")
    if args.batch_size < 1:
        raise ValueError("--batch-size must be at least 1.")
    if args.nproc < 1:
        raise ValueError("--nproc must be at least 1.")
    if args.max_cells is not None and args.max_cells < 1:
        raise ValueError("--max-cells must be at least 1.")
    if not torch.cuda.is_available():
        raise RuntimeError("CUDA is not available. Run this script on the GPU instance with a CUDA-enabled PyTorch installation.")

    dataset = load_from_disk(str(args.dataset_dir))
    required_columns = {"input_ids", "length", *EMBEDDING_LABELS}
    missing_columns = sorted(required_columns.difference(dataset.column_names))
    if missing_columns:
        raise ValueError("Tokenized dataset is missing columns: " + ", ".join(missing_columns))
    with (args.model_dir / "config.json").open() as handle:
        model_config = json.load(handle)

    model_input_size = int(model_config["max_position_embeddings"])
    max_sequence_length = max(dataset["length"])
    if max_sequence_length > model_input_size:
        raise ValueError(f"Dataset sequence length {max_sequence_length} exceeds model limit {model_input_size}.")
    return len(dataset), max_sequence_length

def main() -> None:
    args = parse_args()
    args.model_dir = args.model_dir.expanduser().resolve()
    args.dataset_dir = args.dataset_dir.expanduser().resolve()
    args.output_dir = args.output_dir.expanduser().resolve()

    cell_count, max_sequence_length = validate_inputs(args)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    csv_path = args.output_dir / f"{args.output_prefix}.csv"
    tensor_path = args.output_dir / f"{args.output_prefix}.pt"
    if csv_path.exists() or tensor_path.exists():
        raise FileExistsError(f"Output already exists. Use a different --output-prefix or move the existing files.\nCSV: {csv_path}\nTensor: {tensor_path}")
    device_index = torch.cuda.current_device()
    device_name = torch.cuda.get_device_name(device_index)
    device_memory_gib = torch.cuda.get_device_properties(device_index).total_memory / 1024**3

    selected_cell_count = min(cell_count, args.max_cells) if args.max_cells is not None else cell_count
    print(f"Project: {PROJECT_DIR}")
    print(f"Model: {args.model_dir}")
    print(f"Dataset: {args.dataset_dir}")
    print(f"Cells selected: {selected_cell_count:,} of {cell_count:,}")
    print(f"Maximum sequence length: {max_sequence_length:,}")
    print(f"CUDA device: {device_name} ({device_memory_gib:.1f} GiB)")
    print(f"Forward batch size: {args.batch_size}")

    extractor = EmbExtractor(model_type="Pretrained", num_classes=0, emb_mode="cell", cell_emb_style="mean_pool", filter_data=None, max_ncells=args.max_cells, emb_layer=-1, emb_label=list(EMBEDDING_LABELS), labels_to_plot=None, forward_batch_size=args.batch_size, nproc=args.nproc, model_version="V2")
    embedding_table, embedding_tensor = extractor.extract_embs(model_directory=str(args.model_dir), input_data_file=str(args.dataset_dir), output_directory=str(args.output_dir), output_prefix=args.output_prefix, output_torch_embs=True)
    embedding_tensor = embedding_tensor.detach().cpu()
    torch.save(embedding_tensor, tensor_path)

    expected_rows = selected_cell_count
    if len(embedding_table) != expected_rows:
        raise RuntimeError(f"Expected {expected_rows:,} embeddings, but received {len(embedding_table):,}.")
    if embedding_tensor.shape[0] != expected_rows:
        raise RuntimeError(f"Expected {expected_rows:,} tensor rows, but received {embedding_tensor.shape[0]:,}.")

    print(f"Embedding table shape: {embedding_table.shape}")
    print(f"Embedding tensor shape: {tuple(embedding_tensor.shape)}")
    print(f"Saved CSV: {csv_path}")
    print(f"Saved tensor: {tensor_path}")

if __name__ == "__main__":
    main()
