## Download from FTP folder of Datasus
## ftp://ftp.datasus.gov.br/dissemin/publicos/SIASUS/200801_/Dados/PASP2502d.dbc


import ftplib
from pathlib import Path
import re

FTP_HOST = "ftp.datasus.gov.br"
FTP_DIR = "/dissemin/publicos/SIASUS/200801_/Dados/"
OUTPUT_DIR = Path(__file__).resolve().parent.parent / "input_datasus/PA"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

BASE = 'PA' # Produćão Ambulatorial
UF = {'AC',
    'AL',
    'AP',
    'AM',
    'BA',
    'CE',
    'DF',
    'ES',
    'GO',
    'MA',
    'MS',
    'MT',
    'MG',
    'PA',
    'PB',
    'PR',
    'PE',
    'PI',
    'RJ',
    'RN',
    'RS',
    'RO',
    'RR',
    'SC',
    'SP',
    'SE',
    'TO'
}

ANOS = {'20', '21', '22', '23', '24', '25', '26'} 

# clear the folder before downloading
# for file in OUTPUT_DIR.glob("*.dbc"):
#     file.unlink()

with ftplib.FTP(FTP_HOST) as ftp:
    ftp.login()
    ftp.cwd(FTP_DIR)
    for filename in ftp.nlst("*.dbc"):
        #donwload only files that starts with PA and the year is greater than 2020
        if not filename.startswith(BASE) or filename[2:4] not in UF or filename[4:6] not in ANOS:
            continue
        
        local_path = OUTPUT_DIR / filename
        if not local_path.exists():
            print(f"Downloading {filename}...")
            with open(local_path, "wb") as f:
                ftp.retrbinary(f"RETR {filename}", f.write)
        else:
            print(f"{filename} already exists, skipping download.")