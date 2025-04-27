# app.py

import streamlit as st
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdmolops
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.ensemble import RandomForestClassifier
from sklearn.multioutput import MultiOutputClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report
from scipy.sparse import hstack
from tabulate import tabulate
import streamlit as st  # Make sure this is at the top of your file!

file_path = "DATA SET 55.xlsx"

# Add a loading spinner
with st.spinner('Loading data, extracting features, and training model... 🧠💻 Please wait...'):

    df = pd.read_excel(file_path, sheet_name="Sheet1")






# === Functions ===

def get_graph_features_from_name(name):
    try:
        mol = Chem.MolFromSmiles(Chem.MolToSmiles(Chem.MolFromName(name)))
        if not mol:
            return [0, 0, 0]
        adj = rdmolops.GetAdjacencyMatrix(mol)
        degrees = np.sum(adj, axis=0)
        return [np.mean(degrees), np.max(degrees), int(np.sum(degrees) / 2)]
    except:
        return [0, 0, 0]


def create_adjacency_matrix(elements, bonds):
    size = len(elements)
    matrix = [[0] * size for _ in range(size)]
    for bond in bonds:
        i, j = bond
        matrix[i][j] = 1
        matrix[j][i] = 1
    return matrix


def display_matrix(matrix, elements):
    display = "    " + "  ".join(elements) + "\n"
    for label, row in zip(elements, matrix):
        display += label + " " + " ".join(map(str, row)) + "\n"
    return display


def predict_auto_graph(compound_name):
    graph_feat = get_graph_features_from_name(compound_name)
    vect = vectorizer.transform([compound_name])
    X_input = hstack([vect, np.array(graph_feat).reshape(1, -1)])
    prediction = model.predict(X_input)[0]
    return dict(zip(y.columns, prediction))


# === Load Data and Train Model ===

st.title("Chemical Compound Analyzer 🌟")


@st.cache_resource
def load_and_train():
    file_path = "DATA SET 55.xlsx"
    df = pd.read_excel(file_path, sheet_name="Sheet1")

    graph_feats = []
    for compound in df['Chemical Compound']:
        features = get_graph_features_from_name(compound)
        graph_feats.append(features)

    graph_df = pd.DataFrame(graph_feats, columns=["mean_degree", "max_degree", "total_bonds"])
    vectorizer_local = CountVectorizer()
    X_text = vectorizer_local.fit_transform(df['Chemical Compound'])
    X_full = hstack([X_text, graph_df])

    y_local = df[['Chemical Structure', 'Classification', 'Physical Properties', 'Molecular Weight', 'Melting Point',
                  'Boiling Point']]

    X_train, X_test, y_train, y_test = train_test_split(X_full, y_local, test_size=0.2, random_state=42)
    model_local = MultiOutputClassifier(RandomForestClassifier(n_estimators=100, random_state=42))
    model_local.fit(X_train, y_train)

    return df, vectorizer_local, model_local, y_local


df, vectorizer, model, y = load_and_train()

st.subheader("1. Predict Chemical Properties")

compound_name = st.text_input("Enter a chemical compound name:")
if st.button("Predict"):
    if compound_name.strip() != "":
        prediction = predict_auto_graph(compound_name)
        st.success("Prediction:")
        st.table(prediction)
    else:
        st.error("Please enter a valid compound name.")

st.markdown("---")

st.subheader("2. Create an Adjacency Matrix Manually")

elements = st.text_input("Enter element symbols (space-separated):").split()
num_bonds = st.number_input("Enter number of bonds:", min_value=0, step=1)

bonds = []
if num_bonds > 0:
    st.info("Enter bonds as two space-separated indices (starting from 0)")
    for i in range(num_bonds):
        bond_input = st.text_input(f"Bond {i + 1} (example: 0 1):")
        if bond_input:
            bond = tuple(map(int, bond_input.strip().split()))
            bonds.append(bond)

if st.button("Generate Adjacency Matrix"):
    if elements and bonds:
        matrix = create_adjacency_matrix(elements, bonds)
        st.code(display_matrix(matrix, elements))
    else:
        st.error("Please provide elements and bonds properly!")

st.markdown("---")

st.caption(" RUPESH PROJECT")
