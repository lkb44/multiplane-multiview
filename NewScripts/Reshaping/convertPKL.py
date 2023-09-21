import pickle

# Load the .pkl array
pkl_file = '/Volumes/Crucial X6/MP_COLLAGE/Axial_Fat_NPY_PKL/Feats_Col_Patient-060_ax_ls_ws_11.pkl'
with open(pkl_file, 'rb') as f:
    pkl_data = pickle.load(f)

# Print the dimensions of the data
print("Dimensions of the data:", len(pkl_data))

# Print the contents of the array
print("Contents of the array:")
print(pkl_data)