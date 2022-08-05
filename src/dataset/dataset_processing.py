import pandas as pd
import numpy as np
import os.path as osp
import torch
from torch_geometric.data import Data, Dataset
from torch_geometric.loader import NeighborLoader
import torch_geometric.transforms as T


class CFD_dataset(Dataset):
    def __init__(self, root, local_neighborhood = 0.01, transform=None, pre_filter=None, train_frac=0.8):
        self.pre_transform=T.RadiusGraph(local_neighborhood)
        self.root = root
        self.train_frac = train_frac
        super().__init__(root, transform, self.pre_transform, pre_filter)

    @property
    def raw_file_names(self):
        return ['data_air.csv']

    @property
    def processed_file_names(self):
        return ['data_air1.pt']

    def process(self):
        torch.manual_seed(0)
        for idx, raw_path in enumerate(self.raw_paths):
            data = raw_to_graph()
            # set training mask
            data.train_mask = torch.rand(data.num_nodes) < self.train_frac
            torch.save(data, osp.join(self.processed_dir, f'data_{idx}.pt'))
            idx += 1

    def len(self):
        return len(self.processed_file_names)

    def get(self, idx):
        data = torch.load(osp.join(self.processed_dir, f'data_{idx}.pt'))
        return data
    
    
def raw_to_graph():
    knn_func = T.KNNGraph(k=6)
    df_all = pd.read_csv('./data/raw/data_air.csv')
    df_air = df_all[df_all['Block Name']=='internalMesh'].iloc[::2]
    df_air.head()

    df_tube_all = pd.read_csv('./data/raw/data_tube.csv')
    df_tube = df_tube_all[df_tube_all['Block Name']=='internalMesh'].iloc[::2]
    df_tube.head()

    df_all_co2 = pd.read_csv('./data/raw/data_co2.csv')
    df_co2 = df_all_co2[df_all_co2['Block Name']=='internalMesh'].iloc[::2]
    df_co2.head()
    ##########################################
    df_air['Block Name'] = 0
    df_air['T_start'] = np.min(df_air['T'])
    df_air['P_start'] = 1e5
    ##########################################
    df_co2['Block Name'] = 1
    df_co2['T_start'] = np.max(df_co2['T'])
    df_co2['P_start'] = np.max(df_co2['p'])
    ##########################################
    df_tube['Block Name'] = 2
    df_tube['T_start'] = np.max(df_tube['T'])  #todo
    df_tube['U_0'] = 0 # give tube 0 velocity
    df_tube['U_1'] = 0
    df_tube['U_2'] = 0
    df_tube['P_start'] = 1e5
    n_points  = len(df_air)
    x_air = torch.Tensor(df_air[['U_0','U_1','U_2','T_start','Block Name','P_start']][0:n_points].to_numpy())
    y_air = torch.Tensor(df_air[['T']][0:n_points].to_numpy())
    pos_air = torch.Tensor(df_air[['Points_0','Points_1','Points_2']][0:n_points].to_numpy())
    x_co2 = torch.Tensor(df_co2[['U_0','U_1','U_2','T_start','Block Name','P_start']][0:n_points].to_numpy())
    y_co2 = torch.Tensor(df_co2[['T']][0:n_points].to_numpy())
    pos_co2 = torch.Tensor(df_co2[['Points_0','Points_1','Points_2']][0:n_points].to_numpy())
    x_tube = torch.Tensor(df_tube[['U_0','U_1','U_2','T_start','Block Name','P_start']][0:n_points].to_numpy())
    y_tube = torch.Tensor(df_tube[['T']][0:n_points].to_numpy())
    pos_tube = torch.Tensor(df_tube[['Points_0','Points_1','Points_2']][0:n_points].to_numpy())


    x_air_tube = torch.concat([x_air,x_tube],dim=0)
    pos_air_tube = torch.concat([pos_air,pos_tube],dim=0)
    y_air_tube = torch.concat([y_air,y_tube],dim=0)
    data_air_tube = Data(x = x_air_tube, edge_index = None, pos = pos_air_tube, y = y_air_tube)

    data_air_tube_preprocessed = knn_func(data_air_tube)


    x_tube_co2 = torch.concat([x_tube,x_co2],dim=0)
    pos_tube_co2 = torch.concat([pos_tube,pos_co2],dim=0)
    y_tube_co2 = torch.concat([y_tube,y_co2],dim=0)
    data_tube_co2= Data(x = x_tube_co2, edge_index = None, pos = pos_tube_co2, y = y_tube_co2)
    data_tube_co2_preprocessed = knn_func(data_tube_co2)

    edge_index_two_graphs = torch.concat([data_air_tube_preprocessed.edge_index, data_tube_co2_preprocessed.edge_index+x_air.shape[0] ],dim=1)
    x_air_tube_co2 = torch.concat([x_air,x_tube,x_co2],dim=0)
    pos_air_tube_co2 = torch.concat([pos_air,pos_tube,pos_co2],dim=0)
    y_air_tube_co2 = torch.concat([y_air,y_tube,y_co2],dim=0)
    data_air_tube_co2 = Data(x = x_air_tube_co2, edge_index = edge_index_two_graphs, pos = pos_air_tube_co2, y = y_air_tube_co2)
    return data_air_tube_co2