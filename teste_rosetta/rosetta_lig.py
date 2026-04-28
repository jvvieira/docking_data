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
N_ESTRUTURAS = 200     # tamanho da amostra. 5-10 → teste rápido, 50 → ok, 100–200 → confiável 500+ → publicação
MAX_NUCLEOS = 4        # Quantos núcleos do M1 usar (4 é seguro para o Air)      
        

# Caminhos Rosetta
BASE_ROSETTA = os.path.expanduser("~/rosetta/rosetta-main")
ROSETTA_BIN_LIGAND = f"{BASE_ROSETTA}/source/bin/ligand_dock.default.macosclangrelease"
ROSETTA_DB = f"{BASE_ROSETTA}/database"
# Script de conversão (verifiCAR se a versão usa o sufixo _py3)
MOLFILE_TO_PARAMS = f"{BASE_ROSETTA}/source/scripts/python/public/molfile_to_params.py"

# Diretórios
WORK_DIR = Path(os.path.expanduser("~/rosetta"))
PROTEINS_DIR = WORK_DIR / "proteins"
LIGANDS_DIR  = WORK_DIR / "ligands"    
PARAMS_DIR   = WORK_DIR / "params"
COMPLEX_DIR  = WORK_DIR / "complexes"   
RESULTS_DIR  = WORK_DIR / "results"

for d in [PARAMS_DIR, RESULTS_DIR]: d.mkdir(exist_ok=True)

# =================================================================
# ------------- FUNÇÕES DE PREPARAÇÃO QUÍMICA -------------
# =================================================================

def preparar_ligante_params(lig_path):
    #Converte mol2 para .params e gera o PDB compatível com Rosetta.
    tag = lig_path.stem
    out_params = PARAMS_DIR / f"{tag}.params"
    out_pdb = PARAMS_DIR / f"{tag}_0001.pdb"
    
    if not out_params.exists():
        print ("Parametrizando")
        # Usa LG1 como código de 3 letras padrão para o ligante
        cmd = [
            "python3", MOLFILE_TO_PARAMS, 
            "-n", "LG1", 
            "-p", str(PARAMS_DIR / tag), # O script adiciona .params automaticamente
            str(lig_path)
        ]
        result = subprocess.run(cmd, cwd=PARAMS_DIR, capture_output=True, text=True)
    
        if result.returncode != 0:
            print(f"Erro para {lig_path.name}")
            print(f"Verifique se o arquivo {lig_path.name} começa com @<TRIPOS>MOLECULE e tem hidrogênios.")
            raise Exception(result.stderr)

    return out_params, out_pdb

def preparar_complexo(prot_path, lig_pdb_path):
    # Une Receptor(A) e Ligante(X).
    out_complex = WORK_DIR / "complexes" / f"{prot_path.stem}_{lig_pdb_path.stem}_ini.pdb"
    out_complex.parent.mkdir(exist_ok=True)
    
    with open(out_complex, "w") as out:
        # Escreve Proteína
        with open(prot_path, "r") as f:
            for line in f:
                if line.startswith("ATOM"): out.write(line[:21] + "A" + line[22:])
        out.write("TER\n")
        # Escreve Ligante (PDB gerado pelo molfile_to_params)
        with open(lig_pdb_path, "r") as f:
            for line in f:
                if line.startswith(("ATOM", "HETATM")):
                    out.write(line[:21] + "X" + line[22:])
        out.write("END\n")
    return out_complex

def gerar_resumo_csv_ligantes():
    resumo_path = RESULTS_DIR / "ranking_ligantes_mol2.csv"
    todas_linhas = []
    
    # Procura todos os ficheiros de score na pasta de resultados
    for score_file in RESULTS_DIR.rglob("*.sc"):
        with open(score_file, "r") as f:
            header = None
            for line in f:
                parts = line.split()
                if not parts or parts[0] != "SCORE:": continue
                
                # Identifica o cabeçalho (a linha que contém os nomes das colunas)
                if "total_score" in parts:
                    header = parts
                    continue
                
                # Processa os dados
                if header:
                    # Cria um dicionário mapeando o cabeçalho aos valores da linha
                    d = dict(zip(header, parts))
                    
                    # MÉTRICAS CRÍTICAS PARA LIGANTES:
                    # 1. interface_delta_g: A afinidade teórica (quanto menor, melhor)
                    # 2. ligand_rms_no_super: RMSD do ligante sem sobreposição (precisão da pose)
                    
                    i_dw = d.get("interface_delta_g", "N/A")
                    rms_lig = d.get("ligand_rms_no_super", d.get("rms", "N/A"))

                    todas_linhas.append({
                        "complexo": score_file.parent.name,
                        "interface_delta_g": i_dw,
                        "rms_ligante": rms_lig,
                        "total_score": d.get("total_score", "N/A"),
                        "pose_id": d.get("description", "N/A")
                    })
    
    if todas_linhas:
        # Ordena pelo interface_delta_g (Melhor afinidade primeiro)
        try:
            todas_linhas.sort(key=lambda x: float(x["interface_delta_g"]) if x["interface_delta_g"] != "N/A" else 999.0)
        except ValueError:
            pass

        with open(resumo_path, "w", newline="") as f:
            colunas = ["complexo", "interface_delta_g", "rms_ligante", "total_score", "pose_id"]
            writer = csv.DictWriter(f, fieldnames=colunas)
            writer.writeheader()
            writer.writerows(todas_linhas)
            
        print(f"\nResumo gerado: {resumo_path}")
        print(f"Melhor afinidade: {todas_linhas[0]['interface_delta_g']} REU ({todas_linhas[0]['pose_id']})")    

# =================================================================
# ------------- LOOP PRINCIPAL -------------
# =================================================================

def executar_docking(args):
    cmd, out_dir, tag = args
    log = out_dir / f"{tag}_dock.log"
    try:
        with open(log, "w") as f:
            subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, check=True)
        return f"{tag}: finalizado"
    except:
        return f"{tag}: Erro (ver log)"

def main():
    print(f"Iniciando Docking de Moléculas Pequenas | Modo: {MODO_BUSCA}")
    
    proteinas = list(PROTEINS_DIR.glob("*.pdb"))
    ligantes_mol2 = list(LIGANDS_DIR.glob("*.mol2"))
    fila = []

    for mol2 in ligantes_mol2:
        try:
            params_file, lig_pdb = preparar_ligante_params(mol2)
            
            for prot in proteinas:
                tag = f"{prot.stem}_{mol2.stem}"
                out_dir = RESULTS_DIR / tag
                out_dir.mkdir(exist_ok=True)
                
                # 1. Prepara o complexo PDB
                complexo_pdb = preparar_complexo(prot, lig_pdb)
                

                # Define que usar o modo ligante
                cmd = [ROSETTA_BIN_LIGAND, "-database", ROSETTA_DB]
                cmd += ["-s", str(complexo_pdb)]
                cmd += ["-extra_res_fa", str(params_file)] # OBRIGATÓRIO para moléculas pequenas
                cmd += ["-nstruct", str(N_ESTRUTURAS)]
                cmd += ["-out:path:all", str(out_dir)]
                cmd += ["-overwrite"]


                # Configuração GLOBAL vs LOCAL
                # abbrev2 é o protocolo padrão para ligantes
                cmd += ["-protocol", "abbrev2"]
                
                if MODO_BUSCA == "GLOBAL":
                    # Exploração: translação uniforme de 5A e amostragem de conformeros
                    cmd += ["-docking:uniform_trans", "5.0"] # translação uniforme
                    cmd += ["-docking:ligand:random_conformer"]
                else:
                    # Refinamento: sem translação inicial, apenas minimização local
                    cmd += ["-docking:ligand:dock_pert", "2", "8"] # Perturbação leve de 2A e 8 graus
                
                # 3. Flags de flexibilidade
                cmd += ["-run:preserve_header", "true"]
                cmd += ["-ignore_zero_occupancy", "false"]
                cmd += ["-ex1", "-ex2aro"]
                
                # 4. Saída de dados
                cmd += [
                    "-s", str(complexo_pdb), 
                    "-nstruct", str(N_ESTRUTURAS), 
                    "-out:path:all", str(out_dir), 
                    "-out:file:scorefile", f"{tag}_scores.sc",  # Adicionado o prefixo 'out:file:'
                    "-out:pdb", "true",
                    "-overwrite"
                ]

                fila.append((cmd, out_dir, tag))
                
        except Exception as e:
            print(f"Erro ao processar ligante {mol2.name}: {e}")

    # Execução Parlela
    print(f"Rodando docking em {MAX_NUCLEOS} núcleos...")        
    with ProcessPoolExecutor(max_workers=MAX_NUCLEOS) as ex:
        for res in ex.map(executar_docking, fila):
            print(res)

    gerar_resumo_csv_ligantes()        

if __name__ == "__main__":
    main()