import argparse
import os

def parse_nanostats(nanostats_path):
    n50 = None
    with open(nanostats_path, 'r', encoding='utf-8') as f:
        for line in f:
            if "Read length N50" in line:
                parts = line.split(':')
                if len(parts) > 1:
                    raw_value = parts[1].strip().replace(',', '')
                    try:
                        n50 = float(raw_value)
                        break
                    except ValueError:
                        continue
    
    if n50 is None:
        raise ValueError("Impossible de trouver la valeur 'Read length N50' dans le fichier nanostats.")
    return n50

def count_total_reads(fastq_path):
    print("Calcul du nombre total de reads dans le fichier fastq...")
    with open(fastq_path, 'r', encoding='utf-8') as f:
        line_count = sum(1 for _ in f)
    
    if line_count % 4 != 0:
        print("ATTENTION : Le fichier fastq semble corrompu.")
    return line_count // 4

def get_fastq_record(f):
    while True:
        header = f.readline()
        if not header: break
        seq = f.readline()
        plus = f.readline()
        qual = f.readline()
        yield header + seq + plus + qual

def main():
    parser = argparse.ArgumentParser(description="Script Python qui permet de réaliser un split de fichier FASTQ en fonction de la taille du génome à reconstruire et de la profondeur de lecture souhaitée (depth coverage).")
    
    parser.add_argument('--input_fastq', required=True, help="Chemin du fichier fastq")
    parser.add_argument('--output_folder', required=True, help="Dossier de sortie")
    parser.add_argument('--genome_size', required=True, type=int, help="Taille du génome (bp)")
    parser.add_argument('--nanostats', required=True, help="Fichier nanostats contenant le N50")
    parser.add_argument('--X', required=True, type=float, help="Couverture (Depth coverage) souhaitée par fichier")

    args = parser.parse_args()

    # 1. Création du dossier de sortie
    os.makedirs(args.output_folder, exist_ok=True)
    print(f"Dossier prêt : {args.output_folder}")

    # 2. Récupération du N50
    try:
        n50 = parse_nanostats(args.nanostats)
        print(f"N50 trouvé : {n50}")
    except Exception as e:
        print(f"Erreur lors de la lecture du nanostats : {e}")
        return

    # 3. Calcul du nombre de reads par split
    reads_per_split = int((args.genome_size * args.X) / n50)
    print(f"Objectif : {reads_per_split} reads par fichier (Basé sur X={args.X} et Genome={args.genome_size})")

    # 4. Compter les reads totaux pour gérer le dernier fichier
    total_reads_source = count_total_reads(args.input_fastq)
    print(f"Total reads source : {total_reads_source}")

    # 5. Boucle de split
    file_index = 1
    reads_processed_total = 0
    current_file_reads = 0
    
    # Ouverture du premier fichier de sortie
    current_output_path = os.path.join(args.output_folder, f"split_{file_index}.fastq")
    out_handle = open(current_output_path, 'w', encoding='utf-8')
    
    print("Début du sous-échantillonnage des données...")

    with open(args.input_fastq, 'r', encoding='utf-8') as in_handle:
        for record in get_fastq_record(in_handle):
            out_handle.write(record)
            current_file_reads += 1
            reads_processed_total += 1

            reads_remaining = total_reads_source - reads_processed_total
            
            if current_file_reads >= reads_per_split:
                if reads_remaining > (reads_per_split * 0.5):
                    out_handle.close()
                    file_index += 1
                    current_output_path = os.path.join(args.output_folder, f"split_{file_index}.fastq")
                    out_handle = open(current_output_path, 'w', encoding='utf-8')
                    current_file_reads = 0

    out_handle.close()
    print("Fractionnement terminé.")

    # 6. Vérification finale
    print("--- Vérification ---")
    if reads_processed_total == total_reads_source:
        print(f"SUCCÈS : Tous les reads ont été traités ({reads_processed_total}/{total_reads_source}).")
    else:
        print(f"ERREUR : Incohérence dans le nombre de reads traités ({reads_processed_total} vs {total_reads_source}).")

if __name__ == "__main__":
    main()