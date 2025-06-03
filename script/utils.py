from pathlib import Path
import logging

def setup_directories(directories):
    for dir in directories.values():
        Path(dir).mkdir(parents=True, exist_ok=True)
        
def validate_input(input_dir):
    if not Path(input_dir).exists():
        raise FileNotFoundError(f"Le dossier input {input_dir} n'existe pas")
        
def setup_logging(log_dir):
    logging.basicConfig(
        filename=Path(log_dir)/'pipeline.log',
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
