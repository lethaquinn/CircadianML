import numpy as np
import pandas as pd

class CircadianDataLoader:
    def __init__(self, data_dir='data'):
        self.data_dir = data_dir
        
    def generate_synthetic_data(self, n_genes=100, n_timepoints=12):
        np.random.seed(42)
        time_points = np.arange(0, 48, 48/n_timepoints)
        data = []
        for i in range(n_genes):
            if i < 30:
                period = 24
                phase = np.random.uniform(0, 2*np.pi)
                expression = np.cos(2*np.pi*time_points/period + phase)
            else:
                expression = np.random.normal(0, 1, n_timepoints)
            data.append(expression)
        return pd.DataFrame(data)
    
    def preprocess(self, data):
        return (data - data.mean()) / data.std()
