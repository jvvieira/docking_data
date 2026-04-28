import os
import subprocess
import csv
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor

# =================================================================
# -------------- PAINEL DE CONTROLE --------------
# =================================================================
MODO_BUSCA = "GLOBAL"  # "GLOBAL" para blind docking ou "LOCAL" selecionar a grid
                       # IMPORTANTE: o ligante deve estar no bolsão no complexo inicial antes de ativar o LOCAL
                       # Dessa forma precisa tirar a função de formar o complexo
N_ESTRUTURAS = 100     # tamanho da amostra. 5-10 → teste rápido, 50 → ok, 100–200 → confiável 500+ → publicação
MAX_NUCLEOS = 4        # Quantos núcleos do M1 usar (4 é seguro para o Air)      
   
# Caminhos Rosetta
BASE_ROSETTA = os.path.expanduser("~/rosetta/rosetta-main")
ROSETTA_BIN_DOCK = f"{BASE_ROSETTA}/source/bin/docking_protocol.default.macosclangrelease"
ROSETTA_DB = f"{BASE_ROSETTA}/database"

# Diretórios
WORK_DIR = Path(os.path.expanduser("~/rosetta"))
PROTEINS_DIR = WORK_DIR / "proteins"
LIGANDS_DIR  = WORK_DIR / "ligands"
COMPLEX_DIR  = WORK_DIR / "complexes"
RESULTS_DIR  = WORK_DIR / "results"

# Cria as pastas necessárias
for d in [COMPLEX_DIR, RESULTS_DIR]: 
    d.mkdir(exist_ok=True)

# =================================================================
# ------------- FUNÇÕES DE SUPORTE -------------
# =================================================================

def preparar_complexo(prot_path, lig_path):
    """Une as duas proteínas garantindo Cadeia A para o receptor e B para o parceiro."""
    out_path = COMPLEX_DIR / f"{prot_path.stem}_{lig_path.stem}_complex.pdb"
    
    with open(out_path, "w") as out:
        # Escreve Receptor como Cadeia A
        with open(prot_path, "r") as f:
            for line in f:
                if line.startswith(("ATOM", "HETATM")):
                    out.write(line[:21] + "A" + line[22:])
        out.write("TER\n")
        
        # Escreve Parceiro como Cadeia B
        with open(lig_path, "r") as f:
            for line in f:
                if line.startswith(("ATOM", "HETATM")):
                    out.write(line[:21] + "B" + line[22:])
        out.write("END\n")
    return out_path

def gerar_resumo_csv():
    """Extrai os scores de interface (I_sc) e RMSD para análise de proteína-proteína."""
    resumo_path = RESULTS_DIR / "resumo_final_docking.csv"
    todas_linhas = []
    
    for score_file in RESULTS_DIR.rglob("*_scores.sc"):
        with open(score_file, "r") as f:
            header = None
            for line in f:
                parts = line.split()
                if not parts or parts[0] != "SCORE:": continue
                if "total_score" in parts:
                    header = parts
                    continue
                if header:
                    d = dict(zip(header, parts))
                    # I_sc é a métrica vital para docking de proteínas
                    i_score = d.get("I_sc", d.get("interface_delta_g", "N/A"))
                    rms_val = d.get("Irms", d.get("rms", "N/A"))

                    todas_linhas.append({
                        "complexo": score_file.parent.name,
                        "total_score": d.get("total_score", "N/A"),
                        "I_sc": i_score,
                        "rms": rms_val,
                        "description": d.get("description", "N/A")
                    })
    
    if todas_linhas:
        with open(resumo_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=["complexo", "total_score", "I_sc", "rms", "description"])
            writer.writeheader()
            writer.writerows(todas_linhas)
        print(f"\nResumo gerado com sucesso em: {resumo_path}")

def executar_docking(args):
    """Executa o comando via subprocesso e salva log."""
    cmd, out_dir, tag = args
    log_path = out_dir / f"{tag}.log"
    try:
        with open(log_path, "w") as log:
            subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, check=True)
        return f"{tag}: finalizado"
    except Exception as e:
        return f"{tag}: Erro (verifique o log)"

# =================================================================
# ------------- LOOP PRINCIPAL -------------
# =================================================================

def main():
    fila = []
    proteinas = list(PROTEINS_DIR.glob("*.pdb"))
    ligantes = list(LIGANDS_DIR.glob("*.pdb"))

    if not proteinas or not ligantes:
        print("Erro: Verifique se existem arquivos PDB nas pastas 'proteins' e 'ligands'.")
        return

    print(f"Preparando {len(proteinas) * len(ligantes)} combinações...")

    for p in proteinas:
        for l in ligantes:
            tag = f"{p.stem}_{l.stem}"
            out_dir = RESULTS_DIR / tag
            out_dir.mkdir(exist_ok=True)

            # 1. Prepara o complexo PDB
            complexo_pdb = preparar_complexo(p, l)

            # 2. Define o comando RosettaDock
            cmd = [ROSETTA_BIN_DOCK, "-database", ROSETTA_DB]
            
            if MODO_BUSCA == "GLOBAL":
                # Busca 'cega' com rotação exaustiva do parceiro
                cmd += ["-randomize1", "-randomize2", "-spin"]
            else:
                # Refinamento local com perturbação gaussiana
                cmd += ["-dock_pert", "3", "8"]

            # 3. Flags universais de docking proteína-proteína
            cmd += ["-ignore_zero_occupancy", "false", "-ignore_unrecognized_res", "true"]
            cmd += ["-ex1", "-ex2aro", "-use_input_sc"]
            
            # 4. Flag de saída
            cmd += ["-s", str(complexo_pdb), 
                    "-nstruct", str(N_ESTRUTURAS), 
                    "-out:path:all", str(out_dir), 
                    "-scorefile", f"{tag}_scores.sc", 
                    "-overwrite",
                    "-scoring:config_pw_whitelist", "true", # não roda com FlexPep
                    "-pack_separated", "true",              # Força o cálculo ao separar os parceiros; não roda com FlexPep
                    "-compute_punctuated_callbacks", "true" # Garante que os filtros de interface rodem; não roda com FlexPep
            ]
            
            fila.append((cmd, out_dir, tag))

    # Execução Paralela
    print(f"Rodando docking em {MAX_NUCLEOS} núcleos...")
    with ProcessPoolExecutor(max_workers=MAX_NUCLEOS) as executor:
        resultados = list(executor.map(executar_docking, fila))
    
    for r in resultados:
        print(r)

    # Finalização
    gerar_resumo_csv()

if __name__ == "__main__":
    main()