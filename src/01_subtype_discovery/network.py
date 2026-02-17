import torch
import torch.nn as nn
import torch.nn.functional as F


# number of input_node
num_nodes = [ 64,32,10]


class ImprovedMLP(nn.Module):
    def __init__(self, features, num_classes, dropout_rate=0.3, l2_lambda=1e-4, decoding=False):
        super(ImprovedMLP, self).__init__()
        self.features = features
        self.latent_dim = 10
        self.top_layer = nn.Linear(self.latent_dim, num_classes)
        self.dropout = nn.Dropout(dropout_rate)
        self.l2_lambda = l2_lambda
        self.decoding = decoding
        
        # Add decoder layers 
        if decoding:
            decoder_nodes = [ 10, 32] 
            self.decoder = self._make_decoder_layers(decoder_nodes, final_dim=features[0].in_features)
        
        self._initialize_weights()

    def forward(self, x, return_decoded=False):
        # Encoder pathway
        encoded = self.features(x)
        encoded = self.dropout(encoded)
        
        if self.top_layer:
            classification = self.top_layer(encoded)
            
            if self.decoding:
                # Decoder pathway
                decoded = self.decoder(encoded)
                if return_decoded:
                    return classification, decoded
                return classification  
            return classification
        
        return encoded
    
    def _make_decoder_layers(self, nodes, final_dim):
        layers = []
        input_dim = self.latent_dim
        
        for v in nodes:
            layers.extend([
                nn.Linear(input_dim, v),
                nn.BatchNorm1d(v),
                nn.ReLU(inplace=True),
                nn.Dropout(0.8)
            ])
            input_dim = v
            
        # Final layer to original dimension
        layers.extend([
            nn.Linear(input_dim, final_dim),
            nn.BatchNorm1d(final_dim)
        ])
        
        return nn.Sequential(*layers)

    def _initialize_weights(self):
        for m in self.modules():
            if isinstance(m, nn.BatchNorm1d):
                m.weight.data.fill_(1)
                m.bias.data.zero_()
            elif isinstance(m, nn.Linear):
                nn.init.kaiming_normal_(m.weight, mode='fan_in', nonlinearity='relu')
                m.bias.data.zero_()
    
    def get_l2_regularization(self):
        l2_reg = torch.tensor(0., requires_grad=True)
        for param in self.parameters():
            l2_reg = l2_reg + torch.norm(param, 2)
        return self.l2_lambda * l2_reg


    def get_training_outputs(self, x):
        encoded = self.features(x)
        encoded = self.dropout(encoded)
        classification = self.top_layer(encoded)
        decoded = self.decoder(encoded) if self.decoding else None
        return classification, decoded
	
	
def make_layers_features(num_node, input_dim, bn, dropout_rate=0.2):
    layers = []
    for i, v in enumerate(num_node):
        layers.extend([
            nn.Linear(input_dim, v),
            nn.BatchNorm1d(v) if bn else nn.Identity(),
            nn.ReLU(inplace=True),
            nn.Dropout(dropout_rate if i < len(num_node) - 1 else 0)  
        ])
        input_dim = v
    return nn.Sequential(*layers)

def improved_mlp(input_dim=1193, bn=True, output_dim=4, dropout_rate=0.3, l2_lambda=1e-4, decoding=True):
    model = ImprovedMLP(
        make_layers_features(num_nodes, input_dim, bn=bn, dropout_rate=dropout_rate),
        output_dim,
        dropout_rate=dropout_rate,
        l2_lambda=l2_lambda,
        decoding=decoding
    )
    return model