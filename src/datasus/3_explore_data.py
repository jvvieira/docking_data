import shutil
from pathlib import Path

import pandas as pd
from dbc_reader import DbcReader

INPUT_DIR = Path(__file__).resolve().parent.parent.parent / "input_datasus/PA"
CACHE_DIR = Path(__file__).resolve().parent.parent.parent / "cache_datasus/parquet"
CHUNK_SIZE = 100_000


def peek_columns(dbc_file):
    """Return column names from the first record (no full load)."""
    rec = next(iter(DbcReader(str(dbc_file))))
    return list(rec.keys())


def _cache_path(dbc_file):
    return CACHE_DIR / dbc_file.stem


def dbc_to_parquet(dbc_file, columns=None, chunk_size=CHUNK_SIZE, force=False):
    """Stream a .dbc file to on-disk parquet shards (memory-safe)."""
    cache = _cache_path(dbc_file)
    shards = list(cache.glob("batch_*.parquet"))
    if shards and not force:
        print(f"Using cached parquet at {cache} ({len(shards)} shards)")
        return cache

    if cache.exists():
        shutil.rmtree(cache)
    cache.mkdir(parents=True)

    chunk = []
    batch_idx = 0
    row_count = 0

    for rec in DbcReader(str(dbc_file)):
        if columns:
            chunk.append({col: rec.get(col) for col in columns})
        else:
            chunk.append(rec)
        row_count += 1

        if len(chunk) >= chunk_size:
            pd.DataFrame(chunk).to_parquet(cache / f"batch_{batch_idx:04d}.parquet")
            print(f"  wrote shard {batch_idx} ({row_count:,} rows so far)")
            batch_idx += 1
            chunk = []

    if chunk:
        pd.DataFrame(chunk).to_parquet(cache / f"batch_{batch_idx:04d}.parquet")
        batch_idx += 1

    print(f"Cached {row_count:,} rows in {batch_idx} shards at {cache}")
    return cache


def load_dataframe(dbc_file, columns=None, chunk_size=CHUNK_SIZE, max_rows=None, force=False):
    """Load a .dbc file as a DataFrame without holding all rows in RAM at once."""
    cache = dbc_to_parquet(dbc_file, columns=columns, chunk_size=chunk_size, force=force)

    parts = []
    total = 0
    for shard in sorted(cache.glob("batch_*.parquet")):
        part = pd.read_parquet(shard, columns=columns)
        if max_rows is not None:
            remaining = max_rows - total
            if remaining <= 0:
                break
            part = part.head(remaining)
        parts.append(part)
        total += len(part)
        print(f"  read {shard.name} ({total:,} rows)")

    df = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()
    print(f"Loaded {len(df):,} rows, {len(df.columns)} columns")
    return df


if __name__ == "__main__":
    file_to_process = INPUT_DIR / "PASP2505a.dbc"

    print("Columns:", peek_columns(file_to_process))

    # First run streams the .dbc to parquet shards (~2–3 min for SP files).
    # Later runs reuse cache_datasus/parquet/<stem>/ automatically.
    #
    # Options to control memory:
    #   columns=[...]  — only keep columns you need (smaller cache & DataFrame)
    #   max_rows=N     — load only the first N rows
    #   force=True     — rebuild parquet cache
    my_df = load_dataframe(file_to_process)

    print(my_df.columns)
    print(my_df.head())
