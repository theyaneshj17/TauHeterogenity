import torch
from torch.utils.data import Dataset
import pandas as pd
import numpy as np



class PetDataset(Dataset):
    def __init__(self, datatype):
       
        file_path = datatype
        pet_data = pd.read_csv(file_path)
        pet_data = pet_data.dropna()  # Drop rows with missing values
        

        self.subject_ids = pet_data['SubjectID']  
        vertex_data = pet_data.drop(columns=['SubjectID'])  

      
        normalized_data = vertex_data

        
        self.x_data = torch.from_numpy(normalized_data.values).float()  

    def __getitem__(self, index):
        return self.x_data[index]

    def __len__(self):
        return len(self.x_data)

# Example usage
file_path =  '/N/slate/thjaya/Final/pet.csv'
pet_dataset = PetDataset(file_path)
