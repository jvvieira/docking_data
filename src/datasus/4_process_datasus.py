import pandas as pd
from dbc_reader import DbcReader
from pathlib import Path
import tempfile

INPUT_DIR = Path(__file__).resolve().parent.parent.parent / "input_datasus/PA"
NEEDED_COLS = [
    # 'PA_CODUNI', 
    # 'PA_GESTAO', 
    # 'PA_CONDIC', 
    'PA_UFMUN', 
    # 'PA_REGCT', 
    # 'PA_INCOUT',
    # 'PA_INCURG', 
    # 'PA_TPUPS', 
    # 'PA_TIPPRE', 
    # 'PA_MN_IND', 
    # 'PA_CNPJCPF', 
    # 'PA_CNPJMNT', 
    # 'PA_CNPJ_CC', 
    # 'PA_MVMR', 
    # 'PA_CMPEF', 
    # 'PA_PROC_ID', 
    # 'PA_TPFIN', 
    # 'PA_SUBFIN', 
    # 'PA_NIVCPL', 
    # 'PA_DOCORIG', 
    # 'PA_AUTORIZ', 
    # 'PA_CNSMED', 
    # 'PA_CBOCOD', 
    # 'PA_MOTSAI', 
    # 'PA_OBITO', 
    # 'PA_ENCERR', 
    # 'PA_PERMAN', 
    # 'PA_ALTA', 
    # 'PA_TRANSF', 
    'PA_CIDPRI', 
    'PA_CIDSEC', 
    'PA_CIDCAS', 
    # 'PA_CATEND', 
    'PA_IDADE', 
    # 'IDADEMIN', 
    # 'IDADEMAX', 
    # 'PA_FLIDADE', 
    'PA_SEXO'
    # 'PA_RACACOR', 
    # 'PA_MUNPCN', 
    # 'PA_QTDPRO', 
    # 'PA_QTDAPR', 
    # 'PA_VALPRO', 
    # 'PA_VALAPR', 
    # 'PA_UFDIF', 
    # 'PA_MNDIF', 
    # 'PA_DIF_VAL', 
    # 'NU_VPA_TOT', 
    # 'NU_PA_TOT', 
    # 'PA_INDICA', 
    # 'PA_CODOCO', 
    # 'PA_FLQT', 
    # 'PA_FLER', 
    # 'PA_ETNIA', 
    # 'PA_VL_CF', 
    # 'PA_VL_CL', 
    # 'PA_VL_INC', 
    # 'PA_SRV_C', 
    # 'PA_INENE', 
    # 'PA_NAT_JUR'
]

CHUNK_SIZE = 100_000

tmp_parquet = Path(tempfile.mkdtemp())

CID_PRIORITARIO = "L40"

#clear results folder
# output_dir = Path(__file__).resolve().parent.parent.parent / "outputs"
# if output_dir.exists():
#     for file in output_dir.glob("processed_*.csv"):
#         file.unlink()
# else:    output_dir.mkdir()

def process_file(dbc_file):
    chunk = []
    batch_idx = 0
    for rec in DbcReader(str(dbc_file)):
        chunk.append({col: rec.get(col) for col in NEEDED_COLS})
        if len(chunk) >= CHUNK_SIZE:
            pd.DataFrame(chunk).to_parquet(
                tmp_parquet / f"batch_{batch_idx:04d}.parquet"
            )
            batch_idx += 1
            chunk = []
    if chunk:
        pd.DataFrame(chunk).to_parquet(
            tmp_parquet / f"batch_{batch_idx:04d}.parquet"
        )
        batch_idx += 1

    df = pd.read_parquet(tmp_parquet)
    
    # print(f"Loaded {len(df)} records with columns {list(df.columns)}")
    # print(df[df["PA_CIDPRI"].str.contains(CID_PRIORITARIO, na=False)]["PA_CIDPRI"].value_counts())

    df_final = pd.DataFrame({
        "PA_SEXO": df["PA_SEXO"],
        "PA_IDADA": df["PA_IDADA"],
        "PA_CIDPRI": df["PA_CIDPRI"],
        "PA_CMP": df["PA_CMP"]})
    
    df_final = df_final[df_final["PA_CIDPRI"].str.contains(CID_PRIORITARIO, na=False)]
    df_final.to_csv(f"./outputs/processed_{dbc_file.stem}.csv", index=False)


for dbc_file in sorted(INPUT_DIR.glob("*.dbc")):
    # Ignore if already processed
    output_file = Path(f"./outputs/processed_{dbc_file.stem}.csv")
    if output_file.exists():
        print(f"Skipping {dbc_file.name} (already processed)")
        continue
    print(f"Processing {dbc_file.name}...")
    process_file(dbc_file)
                