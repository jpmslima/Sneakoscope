import streamlit as st
import pandas as pd
import requests
import matplotlib.pyplot as plt
import seaborn as sns
import io
import re

st.set_page_config(page_title="Protein Variant Visualizer", layout="wide")
st.title("🧬 Protein Variant & AlphaMissense Visualizer")
st.markdown("Enter a UniProt ID to fetch data from UniProt and AlphaFold DB and generate pathogenicity plots.")

# User Input
uniprot_id = st.text_input("Enter UniProt Accession (e.g., P63000)", value="P63000").strip().upper()

def get_am_data(uid):
    """Fetches AlphaMissense CSV from AlphaFold DB."""
    # Try the two common URL patterns for AlphaMissense data
    urls = [
        f"https://alphafold.ebi.ac.uk/files/AF-{uid}-F1-aa_substitutions.csv",
        f"https://alphafold.ebi.ac.uk/files/AF-{uid}-F1-hg38.csv"
    ]
    for url in urls:
        response = requests.get(url)
        if response.status_code == 200:
            return pd.read_csv(io.StringIO(response.text))
    return None

def get_uniprot_variants(uid):
    """Fetches variants from UniProt Proteins Variation API."""
    url = f"https://www.ebi.ac.uk/proteins/api/variation/{uid}"
    response = requests.get(url, headers={"Accept": "application/json"})
    if response.status_code == 200:
        data = response.json()
        features = data.get('features', [])
        variant_records = []
        for f in features:
            # We look for single amino acid changes
            wild = f.get('wildType')
            mut = f.get('alternativeSequence')
            pos = f.get('begin')
            if wild and mut and pos and len(wild) == 1 and len(mut) == 1:
                variant_records.append({
                    'protein_variant': f"{wild}{pos}{mut}",
                    'position': int(pos),
                    'source': f.get('sourceType', 'Unknown')
                })
        return pd.DataFrame(variant_records)
    return None

if st.button("Generate Plots"):
    with st.spinner(f"Fetching data for {uniprot_id}..."):
        am_df = get_am_data(uniprot_id)
        var_df = get_uniprot_variants(uniprot_id)

    if am_df is not None and var_df is not None:
        st.success(f"Data retrieved: {len(var_df)} variants found in UniProt.")
        
        # Merge datasets
        merged_df = pd.merge(var_df, am_df, on='protein_variant', how='inner')
        
        # Plot 1: Variant Frequency
        st.subheader("1. Frequency of Variants per Residue")
        fig1, ax1 = plt.subplots(figsize=(12, 5))
        counts = var_df['position'].value_counts().sort_index()
        ax1.bar(counts.index, counts.values, color='skyblue', edgecolor='navy')
        ax1.set_xlabel("Residue Position")
        ax1.set_ylabel("Number of Variants")
        ax1.set_title(f"Variant Distribution for {uniprot_id}")
        st.pyplot(fig1)

        # Plot 2: AlphaMissense Annotations
        st.subheader("2. Pathogenicity of Observed Variants")
        
        # Extract residue number for X-axis
        def get_pos(s):
            m = re.search(r'(\d+)', s)
            return int(m.group(1)) if m else None
        
        am_df['res_num'] = am_df['protein_variant'].apply(get_pos)
        
        col1, col2 = st.columns(2)
        
        with col1:
            st.write("**All Possible Substitutions (Background)**")
            fig2, ax2 = plt.subplots(figsize=(10, 6))
            sns.scatterplot(data=am_df, x='res_num', y='am_pathogenicity', hue='am_class', 
                            palette={'likely_benign': 'green', 'ambiguous': 'grey', 'likely_pathogenic': 'red'},
                            alpha=0.3, s=15, ax=ax2)
            ax2.set_title("Full AlphaMissense Landscape")
            st.pyplot(fig2)

        with col2:
            st.write("**Only UniProt Observed Variants**")
            fig3, ax3 = plt.subplots(figsize=(10, 6))
            sns.scatterplot(data=merged_df, x='position', y='am_pathogenicity', hue='am_class',
                            palette={'likely_benign': 'green', 'ambiguous': 'grey', 'likely_pathogenic': 'red'},
                            s=60, edgecolor='black', ax=ax3)
            ax3.axhline(y=0.34, color='green', linestyle='--', alpha=0.3)
            ax3.axhline(y=0.564, color='red', linestyle='--', alpha=0.3)
            ax3.set_title("Pathogenicity of Clinical/Population Variants")
            st.pyplot(fig3)
            
        # Data Table
        st.subheader("Data Summary")
        st.dataframe(merged_df[['protein_variant', 'source', 'am_pathogenicity', 'am_class']].head(20))
        
    else:
        st.error("Could not retrieve data. Please check the UniProt ID (e.g., use P63000 for Rac1).")