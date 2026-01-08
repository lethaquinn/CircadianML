import os
import pandas as pd
import requests
from typing import Optional


def download_geo_dataset(
    geo_id: str,
    output_dir: str = "data/raw",
    force_download: bool = False
) -> str:

    os.makedirs(output_dir, exist_ok=True)

    output_file = os.path.join(output_dir, f"{geo_id}_series_matrix.txt.gz")

    if os.path.exists(output_file) and not force_download:
        print(f"File already exists: {output_file}")
        return output_file


    base_url = "https://ftp.ncbi.nlm.nih.gov/geo/series"


    series_stub = geo_id[:-3] + "nnn"
    url = f"{base_url}/{series_stub}/{geo_id}/matrix/{geo_id}_series_matrix.txt.gz"

    print(f"Downloading {geo_id} from GEO...")
    print(f"URL: {url}")

    try:
        response = requests.get(url, stream=True, timeout=60)
        response.raise_for_status()

        with open(output_file, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)

        print(f"Download complete: {output_file}")
        return output_file

    except Exception as e:
        print(f"Error downloading {geo_id}: {e}")
        print("You may need to download the data manually from:")
        print(f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={geo_id}")
        return None


def parse_geo_matrix(
    matrix_file: str,
    output_file: Optional[str] = None
) -> pd.DataFrame:

    print(f"Parsing {matrix_file}...")

    with open(matrix_file if not matrix_file.endswith('.gz') else matrix_file, 'r') as f:
        lines = f.readlines()

    data_start = None
    for i, line in enumerate(lines):
        if '!series_matrix_table_begin' in line:
            data_start = i + 1
            break

    if data_start is None:
        raise ValueError("Could not find data table in matrix file")

    data_end = None
    for i in range(data_start, len(lines)):
        if '!series_matrix_table_end' in lines[i]:
            data_end = i
            break

    data_lines = lines[data_start:data_end]

    data = []
    for line in data_lines:
        data.append(line.strip().split('\t'))

    df = pd.DataFrame(data[1:], columns=data[0])

    if '"ID_REF"' in df.columns:
        df.set_index('"ID_REF"', inplace=True)

    print(f"Parsed matrix: {df.shape[0]} genes × {df.shape[1]} samples")

    if output_file:
        df.to_csv(output_file)
        print(f"Saved to {output_file}")

    return df


def list_available_datasets() -> pd.DataFrame:
    
    datasets = [
        {
            'GEO_ID': 'GSE11923',
            'Title': 'Circadian time series of mouse liver',
            'Organism': 'Mus musculus',
            'Timepoints': 48,
            'Description': 'Gene expression in mouse liver over 48h'
        },
        {
            'GEO_ID': 'GSE54650',
            'Title': 'Circadian rhythms in mouse tissues',
            'Organism': 'Mus musculus',
            'Timepoints': 24,
            'Description': 'Multiple tissues across 24h'
        },
    ]

    return pd.DataFrame(datasets)
