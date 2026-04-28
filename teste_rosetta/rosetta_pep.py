import os
import subprocess
import csv
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor

# =================================================================
# -------------- PAINEL DE CONTROLE (PEPTÍDEOS) --------------
# =================================================================
MODO_BUSCA = "GLOBAL"  # "GLOBAL" para blind docking ou "LOCAL" selecionar a grid
                       # IMPORTANTE: o ligante deve estar no bolsão no complexo inicial antes de ativar o LOCAL
                       # Dessa forma precisa tirar a função de formar o complexo
N_ESTRUTURAS = 100     # tamanho da amostra. 5-10 → teste rápido, 50 → ok, 100–200 → confiável 500+ → publicação
MAX_NUCLEOS = 4        # Quantos núcleos do M1 usar (4 é seguro para o Air)      
LIMITE_PEP = 95
         
# Caminhos Rosetta 
BASE_ROSETTA = os.path.expanduser("~/rosetta/rosetta-main")
ROSETTA_BIN_PEP = f"{BASE_ROSETTA}/source/bin/FlexPepDocking.default.macosclangrelease"
ROSETTA_DB = f"{BASE_ROSETTA}/database"

# Diretórios de Trabalho
WORK_DIR = Path(os.path.expanduser("~/rosetta"))
PROTEINS_DIR = WORK_DIR / "proteins"
LIGANDS_DIR  = WORK_DIR / "ligands"
COMPLEX_DIR  = WORK_DIR / "complexes"
RESULTS_DIR  = WORK_DIR / "results"

# Cria as pastas se elas não existirem
for d in [COMPLEX_DIR, RESULTS_DIR]: d.mkdir(exist_ok=True)

# =================================================================
# ------------- FUNÇÕES DE SUPORTE -------------
# =================================================================

def validar_peptideo(pep_path):
    # Conta resíduos (Carbonos Alfa) para garantir o limite.
    count = 0
    with open(pep_path, "r") as f:
        for line in f:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                count += 1
    return count

def preparar_complexo(prot_path, pep_path):
    # Une Receptor(A) e Peptídeo(B). Essencial para o FlexPepDocking.
    out_path = COMPLEX_DIR / f"{prot_path.stem}_{pep_path.stem}_prep.pdb"
    with open(out_path, "w") as out:
        # Receptor -> Cadeia A
        with open(prot_path, "r") as f:
            for line in f:
                if line.startswith(("ATOM", "HETATM")):
                    out.write(line[:21] + "A" + line[22:])
        out.write("TER\n")
        # Peptídeo -> Cadeia B
        with open(pep_path, "r") as f:
            for line in f:
                if line.startswith(("ATOM", "HETATM")):
                    out.write(line[:21] + "B" + line[22:])
        out.write("END\n")
    return out_path

def gerar_resumo_csv():
    # Extrai scores de interface (I_sc) e RMSD, salvando um ranking ordenado.
    resumo_path = RESULTS_DIR / "ranking_final_peptideos.csv"
    todas_linhas = []
    
    # Busca todos os arquivos de score na pasta de resultados
    for score_file in RESULTS_DIR.rglob("*_scores.sc"):
        with open(score_file, "r") as f:
            header = None
            for line in f:
                parts = line.split()
                if not parts or parts[0] != "SCORE:": continue
                
                # Identifica a linha do cabeçalho
                if "total_score" in parts:
                    header = parts
                    continue
                
                # Processa as linhas de dados
                if header:
                    d = dict(zip(header, parts))
                    
                    # Tenta capturar I_sc (essencial para peptídeos)
                    # O FlexPepDocking usa 'I_sc' ou 'interface_delta_g'
                    i_score = d.get("I_sc", d.get("interface_delta_g", "N/A"))
                    
                    # Tenta capturar o RMSD (pode variar o nome da coluna)
                    rms_val = d.get("Irms", d.get("rmsbb", d.get("rms", "N/A")))

                    todas_linhas.append({
                        "complexo": score_file.parent.name,
                        "total_score": d.get("total_score", "0.0"),
                        "I_sc": i_score,
                        "rms": rms_val,
                        "description": d.get("description", "N/A")
                    })
    
    if todas_linhas:
        # ORDENAÇÃO: Coloca os melhores I_sc (mais negativos) no topo
        try:
            todas_linhas.sort(key=lambda x: float(x["I_sc"]) if x["I_sc"] != "N/A" else 999.0)
        except ValueError:
            pass # Mantém a ordem original se houver erro de conversão

        with open(resumo_path, "w", newline="") as f:
            colunas = ["complexo", "I_sc", "rms", "total_score", "description"]
            writer = csv.DictWriter(f, fieldnames=colunas)
            writer.writeheader()
            writer.writerows(todas_linhas)
            
        print(f"\nRanking gerado: {resumo_path}")
        print(f"Melhor pose encontrada: {todas_linhas[0]['description']} com I_sc = {todas_linhas[0]['I_sc']}")


# =================================================================
# ------------- LOOP PREINCIPAL -------------
# =================================================================

def executar_docking(args):
    cmd, out_dir, tag = args
    log_file = out_dir / f"{tag}_exec.log"
    try:
        with open(log_file, "w") as log:
            subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, check=True)
        return f"{tag}: concluído"
    except Exception:
        return f"{tag}: Erro (Ver {log_file.name})"

def main():
    fila = []
    proteinas = list(PROTEINS_DIR.glob("*.pdb"))
    peptideos = list(LIGANDS_DIR.glob("*.pdb"))
    
    if not proteinas or not peptideos:
        print("Aviso: Certifique-se de que há arquivos .pdb nas pastas 'proteins' e 'ligands'.")
        return

    for r in proteinas:
        for p in peptideos:
            # Validação de tamanho
            n_res = validar_peptideo(p)
               
            if n_res > LIMITE_PEP:
                print(f"Pulando {p.name}: {n_res} resíduos excedem o limite.")
                continue

            tag = f"{r.stem}_{p.stem}"
            out_dir = RESULTS_DIR / tag
            out_dir.mkdir(exist_ok=True)
                
            # 1. Gera o complexo
            complexo_pdb = preparar_complexo(r, p)

            # 2. Construção do comando Rosetta
            cmd = [ROSETTA_BIN_PEP, "-database", ROSETTA_DB]
            
                
            # Flags de Protocolo (Conforme Documentação)
            if MODO_BUSCA == "GLOBAL":
            # Busca ab-initio para quando não se conhece a pose exata
                cmd += ["-lowres_abinitio", "-pep_refine"]
            else:
                # Refinamento de alta resolução LOCAL
                cmd += ["-pep_refine"]

            # 3. Flags de Correção e de Amostragem e Estabilidade
            cmd += [
                "-ex1", 
                "-ex2aro", 
                "-use_input_sc",
                "-ignore_zero_occupancy", "false",
                "-ignore_unrecognized_res", "true"
            ]        
                

            # 4. Flag de saída
            cmd += ["-run:preserve_header", "true",
                "-s", str(complexo_pdb),
                "-nstruct", str(N_ESTRUTURAS),
                "-out:path:all", str(out_dir),
                "-scorefile", f"{tag}_scores.sc",
                "-overwrite"
            ]
                
            fila.append((cmd, out_dir, tag))            

    # Execução Paralela
    print(f"Rodando docking em {MAX_NUCLEOS} núcleos...")
    with ProcessPoolExecutor(max_workers=MAX_NUCLEOS) as executor:
        for res in executor.map(executar_docking, fila):
            print(res)

    # Finalização
    gerar_resumo_csv()        

if __name__ == "__main__":
    main()

