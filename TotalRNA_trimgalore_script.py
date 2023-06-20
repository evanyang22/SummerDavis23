import os
import subprocess
import argparse
from concurrent.futures import ThreadPoolExecutor

def run_trimgalore(file_path, output_dir, quality):
    command = ["trim_galore", "-a", "--illumina", "--fastqc", "--quality", quality, "--output_dir", output_dir, file_path]
    subprocess.run(command, check=True)

def main():
    parser = argparse.ArgumentParser(description="Run Trim Galore on a directory of FastQ files.")
    parser.add_argument("input_dir", help="Directory containing FastQ files.")
    parser.add_argument("output_dir", help="Directory to write output files.")
    parser.add_argument("quality", help="Quality cutoff for Trim Galore.")

    args = parser.parse_args()

    files = [os.path.join(args.input_dir, file) for file in os.listdir(args.input_dir) 
            if file.endswith('_R1_001.fastq') or file.endswith('_R1.fastq')]

    with ThreadPoolExecutor(max_workers=os.cpu_count()) as executor:
        executor.map(run_trimgalore, files, [args.output_dir]*len(files), 
                     [args.adapter_sequence]*len(files), [args.quality]*len(files))

if __name__ == "__main__":
    main()




            
