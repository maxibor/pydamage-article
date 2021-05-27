import streamlit as st
import pydamage
import statsmodels
import pkg_resources
import pandas as pd
import numpy as np
import pickle
import gzip

@st.cache()
def load_model():
    """Returns the gml model"""
    model_path = "../models/accuracy_model_v2_python.pickle.gz"
    with gzip.open(model_path, "rb") as mod:
        return pickle.load(mod)


def prepare_df(coverage, contiglength, damage):
    """Prepare dataframe from input parameters

    Args:
        coverage (float): Mean coverage of contig
        contiglength (int): Length of contig
        damage (float): Damage on 5' end of contig

    Returns:
        [pandas DataFrame]: parameters as df
    """
    var_dict = {
        "coverage": [float(coverage)],
        "contiglength": [float(contiglength)],
        "damage": [float(damage)],
    }
    pd_df = pd.DataFrame(var_dict)
    return pd_df


def fit_model(df, model):
    """Fit GLM model to data
    Args:
        df (pandas DataFrame): prepared pydamage results
        model (pypmml model): GLM accuracy model
    """
    return model.predict(df)[0]


def interactive_predict(model):

    st.markdown("""
# [PyDamage](https://github.com/maxibor/pydamage) Accuracy prediction 
See [preprint](https://www.biorxiv.org/content/10.1101/2021.03.24.436838v1) for details
    """)

    coverage = st.number_input(label='Contig Mean Coverage', min_value=0.0, max_value=99999999999.0, value=25.3, step=0.01)
    contiglength = st.number_input(label='Contig Length (bp)', min_value=0, max_value=99999999999, value=10000, step=1)
    damage = st.number_input(label="Damage on 5' end", min_value=0.0, max_value=0.99, value=0.3, step=0.01)

    df = prepare_df(coverage, contiglength, damage)

    pred_acc = fit_model(df, model)

    st.markdown(f"**Predicted accuracy:** `{round(pred_acc,2)}`")

if __name__ == '__main__':
    model = load_model()
    interactive_predict(model)