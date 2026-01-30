from pathlib import Path
import pandas as pd
import requests
from icepop.logging_config import logger
from time import time


class HomologyData:
    def __init__(self, data_dir='./cache', sp='10090'):
        self.data_dir = data_dir
        self.sp = sp

        organism_to_taxid = {
            "dpseudoobscura": "7070",
            "dvirilis": "7091",
            "dmelanogaster": "7227",
            "cbriggsae": "7460",
            "drerio": "7955",
            "ggallus": "9031",
            "hsapiens": "9606",
            "clupus": "9615",
            "fcatus": "9685",
            "btaurus": "9913",
            "mmusculus": "10090",
            "rnorvegicus": "10116",
            "tguttata": "36329",
            "acarolinensis": "185431",
            "olatipes": "237561",
            "xtropicalis": "330879",
            "cintestinalis": "347515",
            "gmorhua": "508771",
            "spurpuratus": "559292",
            "lgigantea": "660027",
            "aqueenslandica": "1340429",
        }

        if self.sp not in organism_to_taxid.keys():
            orgs = list(organism_to_taxid.keys())
            orgs = ', '.join(orgs)
            raise ValueError(f'Species must one of {orgs}')
        self.sp = organism_to_taxid[self.sp]

        self.ortho_file = f'{self.data_dir}/gene_orthologs.gz'
        if not Path(self.ortho_file).exists():
            self.download()

    def download(self):
        Path(self.data_dir).mkdir(parents=True, exist_ok=True)
        self.url = 'https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_orthologs.gz'
        response = requests.get(self.url)
        if response.status_code == 200:
            with open(self.ortho_file, 'wb') as f:
                f.write(response.content)
        else:
            raise RuntimeError(f"Failed to download file: {self.url}")

    def load(self):
        '''get one to one and many to one ortholog from a model species to human'''
        df = pd.read_csv(self.ortho_file, header=0, index_col=None, sep='\t', dtype=str)

        df = df[(df['#tax_id'] == '9606') & (df['Other_tax_id'] == self.sp)]

        human_to_sp_count = df.groupby("GeneID")["Other_GeneID"].nunique()
        sp_to_human_count = df.groupby("Other_GeneID")["GeneID"].nunique()

        # Identify one-to-one and many-to-one mappings
        filtered_df = df[
            (human_to_sp_count == 1).reindex(df["GeneID"]).values &
            (sp_to_human_count >= 1).reindex(df["Other_GeneID"]).values
        ]

        ortho_map = {}
        for human_gene, sp_gene in zip(filtered_df['GeneID'], filtered_df['Other_GeneID']):
            ortho_map.setdefault(human_gene, [])
            ortho_map[human_gene].append(sp_gene)
        return ortho_map
